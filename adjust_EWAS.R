##############################
# D Kiltschewskij, 29/08/23
#
# Adjust EWAS mean effects and variance via Bacon
# Run Prior to meta-analysis
#
##############################


# libraries
library(quantro)
library(data.table)
library(dplyr)
library(R.utils)
library(optparse)
library(tibble)
library(bacon)


## arguments
option_list=list(
  make_option("--string", action="store",default=NA,type="character",help="String common to VMP EWAS files [Required]."),
  make_option("--out", action="store",default=NA,type="character",help="Output string [Required]."),
  make_option("--path", action="store", default=NA,type="character",help="Path to working directory containing EWAS files [required].")
)

opt=parse_args(OptionParser(option_list=option_list))


# set working directory
setwd(paste(opt$path))


# read in EWAS
files<-list.files(".",pattern=paste(opt$string),full.names = T)

read.ewas<-function(x){
  res<-lapply(x, function(y){
    z<-fread(y,data.table = F)
    return(z)
  })
  res<-do.call(rbind,res)
  return(res)
}

ewas<-read.ewas(files)


# adjust for P value inflation via bacon
bc.adj<-function(ew){
  bc.m<-bacon(NULL,ew$Status_Beta,ew$Status_SE)
  ew$Bacon_Status_Beta<-as.numeric(es(bc.m))
  ew$Bacon_Status_SE<-as.numeric(se(bc.m))
  ew$Bacon_Status_P<-as.numeric(pval(bc.m))
  
  bc.lev<-bacon(teststatistics = qnorm(ew$Lev_P/2,lower.tail = F)*sign(ew$Variance_difference))
  bc.bar<-bacon(teststatistics = qnorm(ew$Bar_P/2,lower.tail = F)*sign(ew$Variance_difference))
  bc.FK<-bacon(teststatistics = qnorm(ew$FK_P/2,lower.tail = F)*sign(ew$Variance_difference))
  ew$Bacon_Lev_P<-as.numeric(pval(bc.lev))
  ew$Bacon_Bar_P<-as.numeric(pval(bc.bar))
  ew$Bacon_FK_P<-as.numeric(pval(bc.FK))
  
  return(ew)
}

ewas<-bc.adj(ewas)


# save results
write.table(ewas,file=paste(opt$out),quote=F,sep="\t",row.names = F)
