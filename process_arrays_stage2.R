##############################
# D Kiltschewskij, 09/08/22
#
# Process raw signals from Hannon et al, 2021 raw methylation data.
#
# Stage 2: predict age, cell type proportions and smoking status using filtered/normalised betas/intensities
# Note: Run on neuromol
#
##############################


# libraries
library(quantro)
library(data.table)
library(ggplot2)
library(factoextra)
library(pracma)
library(dplyr)
library(wateRmelon)
library(tibble)
library(optparse)
library(R.utils)
library(tidyr)
library(stringr)
library(ENmix)
library(meffil)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(IlluminaHumanMethylation450kmanifest)
library(IlluminaHumanMethylationEPICanno.ilm10b2.hg19)
library(IlluminaHumanMethylationEPICmanifest)
library(EpiSmokEr)


## arguments
option_list=list(
  make_option("--m_int", action="store",default=NA,type="character",help="Filtered/normalised methylated intensities [required]."),
  make_option("--u_int", action="store",default=NA,type="character",help="Filtered/normalised unmethylated intensities [required]."),
  make_option("--pheno", action="store",default=NA,type="character",help="Text file phenotype data [Required]."),
  make_option("--array", action="store", default=NA,type="character",help="Array type (450K or EPIC) [required]."),
  make_option("--path", action="store", default=NA,type="character",help="Path to working directory [required]."),
  make_option("--out", action="store", default=NA,type="character",help="Output prefix [required].")
)

opt=parse_args(OptionParser(option_list=option_list))
setwd(paste(opt$path))


# read in phenotype data
pheno<-fread(opt$pheno)


# read in intensities
m.int<-fread(opt$m_int)
m.int<-column_to_rownames(m.int,"V1")
u.int<-fread(opt$u_int)
u.int<-column_to_rownames(u.int,"V1")
cat("Read in data")
cat("\n")


# generate methylset
if(opt$array == "450K"){
  mset<-MethylSet(Meth=m.int,Unmeth=u.int,annotation=c(array="IlluminaHumanMethylation450k",annotation="ilmn12.hg19"))  
} else if (opt$array == "EPIC"){
  mset<-MethylSet(Meth=m.int,Unmeth=u.int,annotation=c(array="IlluminaHumanMethylationEPIC",annotation="ilm10b2.hg19"))
}
cat("Generated methylset")
cat("\n")
rm(list=c("m.int","u.int"))


## predict age, cell composition and smoking
# age
age<-wateRmelon::agep(mset)


# cells
if(opt$array == "450K"){
  cells<-wateRmelon::estimateCellCounts.wateRmelon(mset,referencePlatform = "IlluminaHumanMethylation450k",compositeCellType = "Blood")  
} else if (opt$array == "EPIC"){
  cells<-wateRmelon::estimateCellCounts.wateRmelon(mset,referencePlatform = "IlluminaHumanMethylationEPIC",compositeCellType = "Blood")
}


# smoking scores, using a custom sample sheet compatible with epismoker
epi.ss<-pheno
rownames(epi.ss)<-epi.ss$Sentrix
colnames(epi.ss)<-gsub("predicted_sex_bin","sex",colnames(epi.ss))
epi.ss$sex<-ifelse(epi.ss$sex==1,2,ifelse(epi.ss$sex==0,1,NA))
smoking<-epismoker(getBeta(mset,offset=100),samplesheet = epi.ss,method = "all")
colnames(smoking)[1]<-"Sentrix"

# merge
age$Sentrix<-rownames(age)
cells<-as.data.frame(cells)
cells$Sentrix<-rownames(cells)

pheno2<-merge(pheno,age,by="Sentrix")
pheno3<-merge(pheno2,cells,by="Sentrix")
pheno4<-merge(pheno3,smoking,by="Sentrix")


# print results
write.table(pheno4,file=paste(opt$out,"_pheno.txt",sep=""),sep="\t",quote=F, row.names=F)


# print beta distributions
png(paste(opt$out,"_norm_beta_distributions.png",sep=""),units="in",height=5,width=7,res = 600,type="cairo")
matdensity(as.matrix(getBeta(mset,offset=100)),xlab="Beta",ylab="Density")
dev.off()

