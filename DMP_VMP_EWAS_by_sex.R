##############################
# D Kiltschewskij, 21/08/23
#
# Case/Ctrl EWAS using methylation data from Hannon et al, 2021
# print EWASm and EWASv sumstats 
#
##############################


# libraries
library(quantro)
library(data.table)
library(ggplot2)
library(factoextra)
library(ggpubr)
library(pracma)
library(dplyr)
library(e1071)
library(R.utils)
library(optparse)
library(wateRmelon)
library(doParallel)
library(tibble)
library(s20x)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(IlluminaHumanMethylation450kmanifest)
library(IlluminaHumanMethylationEPICanno.ilm10b2.hg19)
library(IlluminaHumanMethylationEPICmanifest)


## arguments
option_list=list(
  make_option("--beta", action="store",default=NA,type="character",help="A matrix containing Beta values [required]."),
  make_option("--pheno", action="store",default=NA,type="character",help="Text file containing all phenotype data of interest [Required]."),
  make_option("--outcome", action="store",default=NA,type="character",help="Outcome variable of interest. Can be SZ or FEP [Required]."),
  make_option("--sex", action="store",default=NA,type="character",help="Sex of interest. Can be 0 (male) or 1 (female) [Required]."),
  make_option("--cores", action="store",default=1,type="numeric",help="Number of cores to use [Required]."),
  make_option("--out", action="store",default=NA,type="character",help="Output string [Required]."),
  make_option("--path", action="store", default=NA,type="character",help="Path to working directory [required].")
)

opt=parse_args(OptionParser(option_list=option_list))
setwd(paste(opt$path))
cl<-makeCluster(opt$cores)
registerDoParallel(cl)
sex_out=ifelse(opt$sex==0,"M",ifelse(opt$sex==1,"F","NA"))

# beta matrix and pheno data
betas.file<-opt$beta
pheno<-fread(paste(opt$pheno))


# call sex of interest
pheno<-pheno[which(pheno$predicted_sex_bin==opt$sex),]


# scale smoking scores and cell type proportions
pheno[,27:34]<-lapply(pheno[,27:34], function(x) c(scale(x))) 


# set Sentrix_ID to factor
pheno$Sentrix_ID<-factor(pheno$Sentrix_ID)


# set factor levels for outcome variable
pheno$Status_Kiltschewskij<-factor(pheno$Status_Kiltschewskij,levels=c("C",opt$outcome)) # ensure controls are baseline


# change "Neu" to "Gran in colnames
colnames(pheno)<-gsub("Neu","Gran",colnames(pheno))


# EWASm function
testCpG<-function(row, phe){
  model<-lm(row ~phe$Status_Kiltschewskij+phe$horvath.age+phe$smokingScore+phe$Sentrix_ID+phe$Sentrix_Position+phe$CD4T+phe$CD8T+phe$NK+phe$Mono+phe$Gran+phe$Bcell,na.action = na.exclude)
  coef<-summary(model)$coefficients
  return(coef[grepl("Status_Kiltschewskij", rownames(coef))][c(1,2,4)])
}
cat(paste("Running EWASm"))
cat("\n")


# EWASv function
varCpG<-function(row, phe){
  # residualise beta values
  resid<-residuals(lm(row ~phe$horvath.age+phe$smokingScore+phe$Sentrix_ID+phe$Sentrix_Position+phe$CD4T+phe$CD8T+phe$NK+phe$Mono+phe$Gran+phe$Bcell,na.action = na.exclude))
  # call case/ctrl means and variance
  ctrl<-resid[which(phe$Status_Kiltschewskij=="C")]
  case<-resid[which(phe$Status_Kiltschewskij==opt$outcome)]
  res<-data.frame("Mean_controls"=mean(ctrl),
                  "Variance_controls"=var(ctrl),
                  "Mean_cases"=mean(case),
                  "Variance_cases"=var(case),
                  "Mean_difference"=mean(case)-mean(ctrl),
                  "Variance_difference"=var(case)-var(ctrl))
  # variance tests on residuals
  lev<-s20x::levene.test(resid ~ phe$Status_Kiltschewskij,show.table = F,digit=200)
  bar<-bartlett.test(resid~phe$Status_Kiltschewskij)
  fli<-fligner.test(resid ~ phe$Status_Kiltschewskij)
  # add variance test results to res
  res$Lev_F<-lev$f.value
  res$Lev_P<-lev$p.value
  res$Bar_K2<-bar$statistic
  res$Bar_P<-bar$p.value
  res$FK_Chi2<-fli$statistic
  res$FK_P<-fli$p.value
  return(res)
}
cat(paste("Running EWASv"))
cat("\n")

# generate annotation
anno.cols<-c("chr","pos","strand","Islands_Name","Relation_to_Island","UCSC_RefGene_Name","UCSC_RefGene_Accession","UCSC_RefGene_Group","Regulatory_Feature_Name","Regulatory_Feature_Group")
anno.450k<-getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19) %>% as.data.frame() %>% dplyr::select(anno.cols) %>% mutate("Probe"=rownames(.))
anno.epic<-getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b2.hg19) %>% as.data.frame() %>% dplyr::select(anno.cols) %>% mutate("Probe"=rownames(.))
anno<-rbind(anno.450k,anno.epic)
anno<-anno[!duplicated(anno$Probe),]

# read in beta file in chunks and run EWAS for binary traits
cols<-fread(betas.file,nrow=2)
cols<-column_to_rownames(cols,"V1")

for (n in seq(0,850000,by=25000)){
  nrow=25000
  if (n == 0){
    skip=n
  } else {
    skip=n+1
  }
  betas<-as.matrix(fread(file=paste(betas.file),skip=skip,nrow=nrow),rownames=1)
  colnames(betas)<-colnames(cols)
  betas<-betas[,match(pheno$Sentrix,colnames(betas))]
  pheno2<-pheno[match(colnames(betas),pheno$Sentrix),]
  res<-foreach(i=1:nrow(betas), .combine=rbind) %dopar%{testCpG(betas[i,],pheno2)}
  res.var<-foreach(i=1:nrow(betas), .combine=rbind) %dopar%{varCpG(betas[i,],pheno2)}
  res<-cbind(res,res.var)
  # format EWAS results
  rownames(res)<-rownames(betas)
  colnames(res)[1:3]<-c("Status_Beta", "Status_SE", "Status_P")
  res<-as.data.frame(res)
  # Convert 0-1 proportions to percentages
  res[,"Status_Beta"]<-res[,"Status_Beta"]*100
  res[,"Status_SE"]<-res[,"Status_SE"]*100
  # Annotate results
  res<-merge(res,anno,by.x=0,by.y=0)
  # Sorting results by P.value
  res<-res[order(res$Status_P),]
  # Save results
  write.table(res, file=paste("EWAS_",opt$out,"_",sex_out,"_chunk_",skip,".txt",sep=""),quote=F,sep="\t",row.names = F)
  cat(paste("Completed analysis of ",n+25000," CpGs \n",sep=""))
}


