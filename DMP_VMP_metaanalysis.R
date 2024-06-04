##############################
# D Kiltschewskij, 08/05/24
#
# Meta analysis of Case/Ctrl EWAS using methylation data from Hannon et al, 2021
# 
#
##############################


# libraries
library(quantro)
library(stringr)
library(data.table)
library(ggplot2)
library(factoextra)
library(pracma)
library(dplyr)
library(e1071)
library(R.utils)
library(optparse)
library(wateRmelon)
library(doParallel)
library(tibble)
library(s20x)
library(metap)
library(meta)
library(bacon)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(IlluminaHumanMethylation450kmanifest)
library(IlluminaHumanMethylationEPICanno.ilm10b2.hg19)
library(IlluminaHumanMethylationEPICmanifest)


## arguments
option_list=list(
  make_option("--chr", action="store",default=NA,type="character",help="Chromosome to analyse [Required]."),
  make_option("--ewas", action="store",default=NA,type="character",help="Comma separated list of adjusted VMP EWAS files [Required]."),
  make_option("--samplesizes", action="store",default=NA,type="character",help="Comma separated list of samples sizes for each EWAS, in the same order as VMP EWAS files [Required]."),
  make_option("--cores", action="store",default=1,type="numeric",help="Number of cores to use [Required]."),
  make_option("--out", action="store",default=NA,type="character",help="Output string [Required]."),
  make_option("--path", action="store", default=NA,type="character",help="Path to working directory [required].")
)

opt=parse_args(OptionParser(option_list=option_list))

# set working directory
setwd(paste(opt$path))

# make cluster
cl<-makeCluster(opt$cores)
registerDoParallel(cl)

# read in EWAS
ewas.files<-str_split(opt$ewas,",") %>% unlist

read.ewas<-function(x,chromosome){
  y<-fread(x,data.table = F)
  z<-y[,c(1:18,(ncol(y)-6):ncol(y))]
  res<-z%>% filter(chr==chromosome)
  return(res)
}

ewas<-lapply(ewas.files,function(x){read.ewas(x,opt$chr)})
names(ewas)<-gsub("^.*\\/","",ewas.files) %>% str_extract(.,"GSE[0-9]{1,10}")

# identify probes in at least 2 studies
probes<-do.call(c,lapply(ewas,function(x){
  y<-x$Row.names
  return(y)}))

probes<-table(probes)
probes<-probes[which(probes>1)] %>% names()

# meta analyse inflation-adjusted mean effects
dmp.ivw<-function(x,prb,ewas){
  cpg<-prb[x]
  
  # collate study effect sizes and se
  means<-do.call(c,lapply(ewas,function(x){
    y<-x$Bacon_Status_Beta[x$Row.names==cpg]
    return(y)
  }))
  
  se<-do.call(c,lapply(ewas,function(x){
    y<-x$Bacon_Status_SE[x$Row.names==cpg]
    return(y)
  }))
  
  # meta analyse
  out<-metagen(means,se,control = list(stepadj=0.5,maxiter=1000))
  
  # call study directions
  directions<-sign(means) %>% gsub(-1,"-",.) %>% gsub(1,"+",.) %>% paste(.,collapse="")
  
  # call cohorts
  cohorts<-do.call(c,lapply(1:length(ewas),function(x){
    y<-ifelse(cpg %in% ewas[[x]]$Row.names,names(ewas)[x],"")
    return(y)
  }))
  cohorts.keep<-which(cohorts!="")
  cohorts2<-paste(cohorts[cohorts.keep],collapse=";")
  
  res<-data.frame(cpg,sum(!is.na(means)),cohorts2,directions,out$TE.fixed,out$seTE.fixed,out$pval.fixed, out$TE.random, out$seTE.random,out$pval.random, out$tau, out$I2, out$Q,1-pchisq(out$Q, out$df.Q))
  colnames(res)<-c("Probe","N_cohorts","Cohorts","Direction", "All_Effect_Fixed", "All_Effect_SE_Fixed", "All_P_Fixed", "All_Effect_Random", "All_Effect_SE_Random","All_P_Random", "All_tau", "All_I2", "All_Q", "All_Het P")
  return(res)
}

dmp.meta<-foreach(i=1:length(probes), .combine=rbind,.export = c("metagen"),.packages=c('dplyr')) %dopar% {dmp.ivw(i,probes,ewas)}


# meta analyse VMPs
vmp.sto<-function(x,prb,ewas,weights){
  cpg<-prb[x]
  
  # collate Levene's test, Bartlett's test and Fligner-Killeen's P values, delta variance directions and weights
  p.lev<-do.call(c,lapply(ewas,function(x){
    y<-x$Bacon_Lev_P[x$Row.names==cpg]
    return(y)
  }))
  
  p.bar<-do.call(c,lapply(ewas,function(x){
    y<-x$Bacon_Bar_P[x$Row.names==cpg]
    return(y)
  }))
  
  p.fk<-do.call(c,lapply(ewas,function(x){
    y<-x$Bacon_FK_P[x$Row.names==cpg]
    return(y)
  }))
  
  signs<-do.call(c,lapply(ewas,function(x){
    y<-sign(x$Variance_difference[x$Row.names==cpg])
    return(y)
  }))
  
  # generate signed Z-scores
  z.lev<-qnorm(p.lev/2,lower.tail = F)*signs
  z.bar<-qnorm(p.bar/2,lower.tail = F)*signs
  z.fk<-qnorm(p.fk/2,lower.tail = F)*signs
  
  # call cohorts
  cohorts<-do.call(c,lapply(1:length(ewas),function(x){
    y<-ifelse(cpg %in% ewas[[x]]$Row.names,names(ewas)[x],"")
    return(y)
  }))
  cohorts.keep<-which(cohorts!="")
  cohorts2<-paste(cohorts[cohorts.keep],collapse=";")
  
  # filter weights for missingness
  weights<-weights[cohorts.keep]
  
  # meta analyse
  z.out.lev<-sum(z.lev*weights) / sqrt(sum(weights^2))
  z.out.bar<-sum(z.bar*weights) / sqrt(sum(weights^2))
  z.out.fk<-sum(z.fk*weights) / sqrt(sum(weights^2))
  
  p.out.lev<-pnorm(abs(z.out.lev),lower.tail=F)*2
  p.out.bar<-pnorm(abs(z.out.bar),lower.tail=F)*2
  p.out.fk<-pnorm(abs(z.out.fk),lower.tail=F)*2
  
  
  # call study directions
  directions<-signs%>% gsub(-1,"-",.) %>% gsub(1,"+",.) %>% paste(.,collapse="")
  
  res<-data.frame(cpg,sum(!is.na(z.lev)),cohorts2,directions,z.out.lev,p.out.lev,z.out.bar,p.out.bar,z.out.fk,p.out.fk)
  colnames(res)<-c("Probe","N_cohorts","Cohorts","Direction", "All_Lev_Z","All_Lev_P","All_Bar_Z","All_Bar_P","All_FK_Z","All_FK_P")
  return(res)
}

weights<-str_split(opt$samplesizes,",") %>% unlist %>% as.numeric
weights<-sqrt(weights)
vmp.meta<-foreach(i=1:length(probes), .combine=rbind, .packages=c('dplyr')) %dopar% {vmp.sto(i,probes,ewas,weights)}

# merge results
dmp.vmp.meta<-merge(dmp.meta,vmp.meta,by="Probe")
colnames(dmp.vmp.meta)<-gsub("\\.x","_DMP",colnames(dmp.vmp.meta))
colnames(dmp.vmp.meta)<-gsub("\\.y","_VMP",colnames(dmp.vmp.meta))

# annotate results
anno.cols<-c("chr","pos","strand","Islands_Name","Relation_to_Island","UCSC_RefGene_Name","UCSC_RefGene_Accession","UCSC_RefGene_Group","Regulatory_Feature_Name","Regulatory_Feature_Group")
anno.450k<-getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19) %>% as.data.frame() %>% dplyr::select(anno.cols) %>% mutate("Probe"=rownames(.))
anno.epic<-getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b2.hg19) %>% as.data.frame() %>% dplyr::select(anno.cols) %>% mutate("Probe"=rownames(.))
anno<-rbind(anno.450k,anno.epic)
anno<-anno[!duplicated(anno$Probe),]
dmp.vmp.meta<-merge(dmp.vmp.meta, anno, by="Probe")
cat("Finished annotating results")
cat("\n")

# save results
write.table(dmp.vmp.meta,file=paste(opt$out,"_",opt$chr,"_DMP_VMP_meta.txt",sep=""),sep="\t",quote=F,row.names = T)
