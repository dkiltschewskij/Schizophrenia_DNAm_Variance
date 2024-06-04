##############################
# D Kiltschewskij
#
# Call phenotypes from GEO series matrix files
# note that each data frame is ordered differently, thus full processing cannot be automated in one function
# all directories have been censored for security purposes
#
##############################


# libraries
library(data.table)
library(dplyr)
library(janitor)
library(tidyr)
library(tibble)

# setwd


# all fields required in the final phenotype files
fields<-c("GSE_ID","GEO_ID","Species","Platform_ID","Sentrix","Sentrix_ID","Sentrix_Position","Tissue","Age","Sex","Diagnosis","Status_Hannon","Exclude_Hannon","Status_Kiltschewskij","Exclude_Kiltschewskij")


# list all GEO series matrix files
series.mat.files<-list.files(pattern="*series_matrix.txt.gz",recursive = T)


# function to read in files
read.geo<-function(x){
  geo<-fread(paste("zcat < ",x," | sed -n '/!Sample_title/,$p'"),fill=T,sep = "\t",header=F)
  geo2<-as.data.frame(t(geo))
  geo2<-row_to_names(geo2,row_number = 1) # set first row to column names
  return(geo2)
}


## GSE147221
# read in file
GSE147221<-read.geo(paste(series.mat.files[[1]]))

# keep columns of interest
GSE147221<-GSE147221[,c(1,2,8,9,10,11,12,13,22,24)]

# delimit sentrix data
GSE147221<-separate(data = GSE147221, col = "!Sample_title", into = c("Sentrix"), sep = ":")

# generate columns for sentrix_id and sentrix_pos
GSE147221$Sentrix_ID<-sapply(strsplit(GSE147221$Sentrix,split = "_"),'[',1)
GSE147221$Sentrix_Position<-sapply(strsplit(GSE147221$Sentrix,split = "_"),'[',2)

# generate column for GEO_ID
GSE147221$GEO_ID<-rep("GSE147221",nrow(GSE147221))

# rename and reorder columns
colnames(GSE147221)<-c("Sentrix","GSE_ID","Tissue","Species","Sex","Diagnosis","Status_Hannon","Age","Exclude_Hannon","Platform_ID","Sentrix_ID","Sentrix_Position","GEO_ID")
GSE147221<-GSE147221[,match(fields[1:13],colnames(GSE147221))]

# code sex to 0,1,NA
GSE147221$Sex<-ifelse(GSE147221$Sex=="Sex: M", 0, 
                 ifelse(GSE147221$Sex == "Sex: F", 1, "NA"))

# remove superfluous prefixes
GSE147221$Diagnosis<-gsub("consensus_diagnosis: ","",GSE147221$Diagnosis)
GSE147221$Status_Hannon<-gsub("status: ","",GSE147221$Status_Hannon)
GSE147221$Age<-gsub("age: ","",GSE147221$Age)

# define case-control status
GSE147221$Status_Kiltschewskij<-ifelse(grepl("[Ss]chizophrenia|[Ss]chizo",GSE147221$Diagnosis),"SZ",
                                  ifelse(grepl("[Cc]ontrol",GSE147221$Diagnosis),"C","NA"))

# flag samples with no case/ctrl status and those excluded in the Hannon analysis
GSE147221$Exclude_Kiltschewskij<-ifelse(grepl("NA",GSE147221$Status_Kiltschewskij)|grepl("excluded",GSE147221$Exclude_Hannon),"Exclude","Keep")

# change Exclude_Hannon to Keep/Exclude
GSE147221$Exclude_Hannon<-ifelse(grepl("excluded",GSE147221$Exclude_Hannon),"Exclude","Keep")

# write to phenotype file with all samples
write.table(GSE147221,file="GSE147221_pheno_all.txt",sep="\t",row.names=F,quote=F)


# write phenotype file with cases and controls only
GSE147221.keep<-GSE147221[GSE147221$Exclude_Kiltschewskij=="Keep",]
write.table(GSE147221.keep,file="GSE147221_pheno_szcasectrl.txt",sep="\t",row.names=F,quote=F)
rm(list=ls(pattern = "GSE147221"))



## GSE152026
# read in file
GSE152026<-read.geo(paste(series.mat.files[[2]]))

# keep columns of interest
GSE152026<-GSE152026[,c(1,2,8,9,10,11,12,24)]

# delimit sentrix data
GSE152026<-separate(data = GSE152026, col = "!Sample_title", into = c("Sentrix"), sep = " ")

# generate columns for sentrix_id and sentrix_pos
GSE152026$Sentrix_ID<-sapply(strsplit(GSE152026$Sentrix,split = "_"),'[',1)
GSE152026$Sentrix_Position<-sapply(strsplit(GSE152026$Sentrix,split = "_"),'[',2)

# generate column for GEO_ID
GSE152026$GEO_ID<-rep("GSE152026",nrow(GSE152026))

# dummy column for Status_Hannon
GSE152026$`!Sample_characteristics_ch1`<-gsub("phenotype: ","",GSE152026$`!Sample_characteristics_ch1`)
GSE152026$Status_Hannon<-GSE152026$`!Sample_characteristics_ch1`

# rename and reorder columns
colnames(GSE152026)<-c("Sentrix","GSE_ID","Tissue","Species","Diagnosis","Sex","Age","Platform_ID","Sentrix_ID","Sentrix_Position","GEO_ID","Status_Hannon")
GSE152026<-GSE152026[,match(fields[1:12],colnames(GSE152026))]

# dummy Exclude_Hannon column
GSE152026$Exclude_Hannon<-"Keep"

# code sex to 0,1,NA
GSE152026$Sex<-ifelse(GSE152026$Sex=="Sex: Male", 0, 
                      ifelse(GSE152026$Sex == "Sex: Female", 1, "NA"))

# remove superfluous prefixes
GSE152026$Age<-gsub("age: ","",GSE152026$Age)

# define case-control status
GSE152026$Status_Kiltschewskij<-ifelse(grepl("[Cc]ase",GSE152026$Diagnosis),"FEP",
                                       ifelse(grepl("[Cc]ontrol",GSE152026$Diagnosis),"C","NA"))

# flag samples with no case/ctrl status and those excluded in the Hannon analysis
GSE152026$Exclude_Kiltschewskij<-ifelse(grepl("NA",GSE152026$Status_Kiltschewskij)|grepl("excluded",GSE152026$Exclude_Hannon),"Exclude","Keep")

# write to phenotype file with all samples
write.table(GSE152026,file="GSE152026_pheno_all.txt",sep="\t",row.names=F,quote=F)


# write phenotype file with cases and controls only
GSE152026.keep<-GSE152026[GSE152026$Exclude_Kiltschewskij=="Keep",]
write.table(GSE152026.keep,file="GSE152026_pheno_fepcasectrl.txt",sep="\t",row.names=F,quote=F)
rm(list=ls(pattern = "GSE152026"))



## GSE152027
# read in file
GSE152027<-read.geo(paste(series.mat.files[[3]]))

# keep columns of interest
GSE152027<-GSE152027[,c(1,2,8,9,10,11,12,24)]

# delimit sentrix data
GSE152027<-separate(data = GSE152027, col = "!Sample_title", into = c("Sentrix"), sep = " ")

# generate columns for sentrix_id and sentrix_pos
GSE152027$Sentrix_ID<-sapply(strsplit(GSE152027$Sentrix,split = "_"),'[',1)
GSE152027$Sentrix_Position<-sapply(strsplit(GSE152027$Sentrix,split = "_"),'[',2)

# generate column for GEO_ID
GSE152027$GEO_ID<-rep("GSE152027",nrow(GSE152027))

# dummy column for Status_Hannon
GSE152027$`!Sample_characteristics_ch1`<-gsub("status: ","",GSE152027$`!Sample_characteristics_ch1`)
GSE152027$Status_Hannon<-GSE152027$`!Sample_characteristics_ch1`

# rename and reorder columns
colnames(GSE152027)<-c("Sentrix","GSE_ID","Tissue","Species","Diagnosis","Sex","Age","Platform_ID","Sentrix_ID","Sentrix_Position","GEO_ID","Status_Hannon")
GSE152027<-GSE152027[,match(fields[1:12],colnames(GSE152027))]

# fill missing age value in GSE152027, which appears to have been caused by a missing cell for this sample's sex
GSE152027[58,9]<-GSE152027[58,10]

# dummy Exclude_Hannon column
GSE152027$Exclude_Hannon<-"Keep"

# code sex to 0,1,NA
GSE152027$Sex<-ifelse(GSE152027$Sex=="gender: M", 0, 
                      ifelse(GSE152027$Sex == "gender: F", 1, "NA"))

# remove superfluous prefixes
GSE152027$Age<-gsub("ageatbloodcollection: ","",GSE152027$Age)

# define case-control status
GSE152027$Status_Kiltschewskij<-ifelse(GSE152027$Status_Hannon=="SCZ", "SZ", 
                                       ifelse(GSE152027$Status_Hannon == "FEP", "FEP", 
                                              ifelse(GSE152027$Status_Hannon=="CON","C","NA")))

# flag samples with no case/ctrl status and those excluded in the Hannon analysis
GSE152027$Exclude_Kiltschewskij<-ifelse(grepl("NA",GSE152027$Status_Kiltschewskij)|grepl("excluded",GSE152027$Exclude_Hannon),"Exclude","Keep")

# write to phenotype file with all samples
write.table(GSE152027,file="GSE152027_pheno_all.txt",sep="\t",row.names=F,quote=F)

# write phenotype file with cases and controls only
GSE152027.keep<-GSE152027[GSE152027$Exclude_Kiltschewskij=="Keep",]
GSE152027.keep.FEP<-GSE152027.keep[GSE152027.keep$Status_Kiltschewskij!="SZ",]
GSE152027.keep.SZ<-GSE152027.keep[GSE152027.keep$Status_Kiltschewskij!="FEP",]

write.table(GSE152027.keep,file="GSE152027_pheno_fepszcasectrl.txt",sep="\t",row.names=F,quote=F)
write.table(GSE152027.keep.FEP,file="GSE152027_pheno_fepcasectrl.txt",sep="\t",row.names=F,quote=F)
write.table(GSE152027.keep.SZ,file="GSE152027_pheno_szcasectrl.txt",sep="\t",row.names=F,quote=F)
rm(list=ls(pattern = "GSE152027"))


## GSE84727
# read in file
GSE84727<-read.geo(paste(series.mat.files[[5]]))

# keep columns of interest
GSE84727<-GSE84727[,c(2,8,9,10,11,12,13,23)]

# delimit sentrix data
GSE84727$`!Sample_characteristics_ch1`<-gsub("sentrixids: ","",GSE84727$`!Sample_characteristics_ch1`)
colnames(GSE84727)[4]<-"Sentrix"

# generate columns for sentrix_id and sentrix_pos
GSE84727$Sentrix_ID<-sapply(strsplit(GSE84727$Sentrix,split = "_"),'[',1)
GSE84727$Sentrix_Position<-sapply(strsplit(GSE84727$Sentrix,split = "_"),'[',2)

# generate column for GEO_ID
GSE84727$GEO_ID<-rep("GSE84727",nrow(GSE84727))

# dummy column for Status_Hannon
GSE84727$`!Sample_characteristics_ch1.3`<-gsub("disease_status: ","",GSE84727$`!Sample_characteristics_ch1.3`)
GSE84727$Status_Hannon<-ifelse(GSE84727$`!Sample_characteristics_ch1.3`==1,"C",
                               ifelse(GSE84727$`!Sample_characteristics_ch1.3`==2,"SZ","NA"))

# rename and reorder columns
colnames(GSE84727)<-c("GSE_ID","Tissue","Species","Sentrix","Sex","Age","Diagnosis","Platform_ID","Sentrix_ID","Sentrix_Position","GEO_ID","Status_Hannon")
GSE84727<-GSE84727[,match(fields[1:12],colnames(GSE84727))]

# dummy Exclude_Hannon column
GSE84727$Exclude_Hannon<-"Keep"

# code sex to 0,1,NA
GSE84727$Sex<-ifelse(GSE84727$Sex=="Sex: M", 0, 
                     ifelse(GSE84727$Sex == "Sex: F", 1, "NA"))

# remove superfluous prefixes
GSE84727$Age<-gsub("age: ","",GSE84727$Age)

# define case-control status
GSE84727$Status_Kiltschewskij<-GSE84727$Status_Hannon

# flag samples with no case/ctrl status and those excluded in the Hannon analysis
GSE84727$Exclude_Kiltschewskij<-ifelse(grepl("NA",GSE84727$Status_Kiltschewskij)|grepl("excluded",GSE84727$Exclude_Hannon),"Exclude","Keep")

# write to phenotype file with all samples
write.table(GSE84727,file="GSE84727_pheno_all.txt",sep="\t",row.names=F,quote=F)

# write phenotype file with cases and controls only
GSE84727.keep<-GSE84727[GSE84727$Exclude_Kiltschewskij=="Keep",]
write.table(GSE84727.keep,file="GSE84727_pheno_szcasectrl.txt",sep="\t",row.names=F,quote=F)
rm(list=ls(pattern = "GSE84727"))






