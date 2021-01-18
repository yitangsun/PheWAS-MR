'%ni%' <- Negate('%in%')
library(plyr)
library(dplyr)
library(tidyverse)
library(PheWAS)
library(readr)

setwd("/Users/yelab/UKBiobank/")
bd <- read.table("/Users/yelab/UKBiobank/ukb34137.tab", header=TRUE, sep="\t")
bd <- read.table("ukb34137.tab", header=TRUE, sep="\t")
source('/Users/yelab/UKBiobank/ukb34137.r') #generates bd (UKB dataset)

bd_kinship <- read.table("/Users/yelab/UKBiobank/ukb48818_rel_s488282_output.dat", header=F, sep=" ")



bd_pheno<- bd %>% select(f.eid, f.31.0.0, f.21000.0.0, f.40007.0.0, f.21003.0.0,f.53.0.0,
                         f.30690.0.0, f.30780.0.0,f.30760.0.0, f.30870.0.0, f.30710.0.0,f.4080.0.0,f.4079.0.0, f.2443.0.0,f.6150.0.0,f.22020.0.0,
                         f.22027.0.0, f.22019.0.0, 
                         f.22021.0.0, f.22006.0.0,f.22001.0.0 ,f.189.0.0, f.22000.0.0, f.54.0.0)
colnames(bd_pheno)<-c("FID", "Sex", "Race", "Death_Age",  "Age", "Date_of_attending",
                      "Tot_Chol", "LDL", "HDL","TAGs","CRP","SBP","DBP", "Diabetes","Heart_Attack", "Used_in_PCA",
                      "Outliers_for_het_or_missing", "SexchrAneuploidy",
                      "Genetic_kinship", "Genetic_ethnic_grouping","Genetic_Sex","Townsend", "Array", "Assessment_centres")
bd_pheno<-as_tibble(bd_pheno)

#####################kinship $ race
#bd_pheno$kinship=bd$f.22021.0.0
#bd_pheno$race=bd$f.22006.0.0

#Sex
bd_pheno$Sex<-mapvalues(as.character(bd_pheno$Sex), 
                        c("Male", "Female"), c(1,0))
bd_pheno$Genetic_Sex<-mapvalues(as.character(bd_pheno$Genetic_Sex), 
                                c("Male", "Female"), c(1,0))

#Death_Age
bd_pheno$Death_Age<-replace_na(bd_pheno$Death_Age, "-9")

# Finalize/write  pheno file -----------------------------------------

bd_pheno$IID<-bd_pheno$FID


bd_new <- read.table("ukb37330.tab",  header=TRUE, sep="\t")
bd_new1 <- read.table("ukb40720.tab", header=TRUE, sep="\t")
bd_new2 <- read.table("ukb43083.tab", header=TRUE, sep="\t")

bd_pheno3<- bd_new2 %>% select(f.eid,f.30080.0.0, f.30100.0.0, f.30110.0.0, f.30090.0.0, f.30010.0.0, f.30040.0.0, f.30270.0.0,f.30030.0.0,
                          f.30050.0.0, f.30060.0.0, f.30020.0.0,f.30070.0.0, f.30250.0.0, f.30260.0.0, f.30240.0.0,
                          f.30280.0.0, f.30300.0.0, f.30290.0.0, f.30170.0.0, f.30230.0.0, f.30130.0.0, f.30140.0.0,
                          f.30150.0.0,f.30160.0.0, f.30120.0.0, f.30000.0.0, f.30190.0.0, f.30200.0.0, f.30210.0.0,
                          f.30220.0.0, f.30180.0.0)
colnames(bd_pheno3)<-c("FID", "Platelet_count", "Mean_platelet_(thrombocyte)_volume", "Platelet_distribution_width", 
"Platelet_crit", "Red_blood_cell_(erythrocyte)_count", "Mean_corpuscular_volume", "Mean_sphered_cell_volume", 
"Haematocrit_percentage", "Mean_corpuscular_haemoglobin", "Mean_corpuscular_haemoglobin_concentration", 
"Haemoglobin_concentration", "Red_blood_cell_(erythrocyte)_distribution_width", "Reticulocyte_count", 
"Mean_reticulocyte_volume", "Reticulocyte_percentage", "Immature_reticulocyte_fraction", 
"High_light_scatter_reticulocyte_count", "High_light_scatter_reticulocyte_percentage", 
"Nucleated_red_blood_cell_count", "Nucleated_red_blood_cell_percentage", "Monocyte_count", "Neutrophill_count", 
"Eosinophill_count", "Basophill_count", "Lymphocyte_count", "White_blood_cell_(leukocyte)_count", 
"Monocyte_percentage", "Neutrophill_percentage", "Eosinophill_percentage", "Basophill_percentage", 
"Lymphocyte_percentage")
bd_pheno3<-as_tibble(bd_pheno3)
bd_pheno<-bd_pheno %>% left_join(bd_pheno3, by= "FID")


bd_pheno1<- bd %>% select(f.eid,f.30510.0.0, f.30515.0.0, f.30500.0.0, f.30505.0.0, f.30520.0.0, f.30525.0.0, f.30530.0.0, f.30535.0.0)
colnames(bd_pheno1)<-c("FID", "Creatinine", "Creatinine_result_flag","Microalbumin","Microalbumin_result_flag","Potassium","Potassium_result_flag","Sodium","Sodium_result_flag")
bd_pheno1<-as_tibble(bd_pheno1)
bd_pheno<-bd_pheno %>% left_join(bd_pheno1, by= "FID")



today <- Sys.Date()
for (data_nu in 1:nrow(bd_pheno)) {
  bd_pheno$Age_now[data_nu]=bd_pheno$Age[data_nu]+(as.numeric(as.Date(today)-as.Date(bd_pheno$Date_of_attending[data_nu])))%/%365 
}
# bd_pheno3<- bd_pheno %>% select(FID,Age, Date_of_attending, Death_Age, Age_now)
# save(bd_pheno3,file = "Age_09_06_unQC_bd_pheno.RData")
# load('Age_09_06_unQC_bd_pheno.RData')

#input_df1 <- read.table('/Users/yelab/Covid-19/GWAS summary/ukb-covid19.20200523.txt',header=T)
input_df2 <- read.table('/Users/yelab/Covid-19/GWAS summary/ukb-covid19.20200607.txt',header=T)
input_df1 <- read.table('covid19_result.20200830.txt',header=T)


## Total 4510 Result (COVID19) 1326 Inpatient COVID19 932 # Inpatient 3186
eid=unique(input_df1$eid)

#####Total: 6117 Positive: 1474 origin: 4448 serious: 991

xyy2=input_df1[input_df1$result==1 ,]
xyy3=input_df1[input_df1$origin==1 ,]



xy=unique(xyy2$eid)
xy2=unique(xyy3$eid)

input_df1_origin=data.frame(eid)
input_df1_origin$origin=ifelse(input_df1_origin$eid %in%xy2,1,0)
input_df1_origin$result=ifelse(input_df1_origin$eid %in%xy,1,0)
input_df1_origin$In_patient_positive=ifelse(input_df1_origin$result==1 & input_df1_origin$origin==1,1,0)

xy1=input_df1_origin[input_df1_origin$In_patient_positive==1 ,]
xy2=input_df1_origin[input_df1_origin$origin==1 ,]
xy3=input_df1_origin[input_df1_origin$result==1 ,]

bd_pheno$In_patient_positive=ifelse(bd_pheno$FID %in%xy1$eid,1,0)
bd_pheno$origin=ifelse(bd_pheno$FID %in%xy2$eid,1,0)
bd_pheno$result=ifelse(bd_pheno$FID %in%xy3$eid,1,0)
bd_pheno$attend=ifelse(bd_pheno$FID %in%eid,1,0)


bd_pheno$result <- as.integer(bd_pheno$result)
bd_pheno$origin <- as.integer(bd_pheno$origin)
bd_pheno$In_patient_positive <- as.integer(bd_pheno$In_patient_positive)


# Generate pheno 4882le -------------------------------------------
#######BMI
bd_pheno1<- bd_new %>% select(f.eid,f.21001.0.0)
colnames(bd_pheno1)<-c("FID", "BMI")
bd_pheno1<-as_tibble(bd_pheno1)
bd_pheno<-bd_pheno %>% left_join(bd_pheno1, by= "FID")


# Generate Alcohol_status Smoking_status

bd_pheno2<- bd %>% select(f.eid, f.20117.0.0)
colnames(bd_pheno2)<-c("FID", "Alcohol_status")

bd_pheno2<-as_tibble(bd_pheno2)

bd_pheno<-bd_pheno %>% left_join(bd_pheno2, by= "FID")

bd_pheno$Alcohol_status<-mapvalues(as.character(bd_pheno$Alcohol_status), 
                                   c("Prefer not to answer", "Never", "Previous", "Current"), c(-3,0,1,2))


bd_pheno1<- bd_new1 %>% select(f.eid,f.20116.0.0)
colnames(bd_pheno1)<-c("FID", "Smoking_status")

bd_pheno1<-as_tibble(bd_pheno1)


bd_pheno<-bd_pheno %>% left_join(bd_pheno1, by= "FID")




########blood cell
bd_pheno3<- bd_new2 %>% select(f.eid,f.30080.0.0, f.30100.0.0, f.30110.0.0, f.30090.0.0, f.30010.0.0, f.30040.0.0, f.30270.0.0,f.30030.0.0,
                               f.30050.0.0, f.30060.0.0, f.30020.0.0,f.30070.0.0, f.30250.0.0, f.30260.0.0, f.30240.0.0,
                               f.30280.0.0, f.30300.0.0, f.30290.0.0, f.30170.0.0, f.30230.0.0, f.30130.0.0, f.30140.0.0,
                               f.30150.0.0,f.30160.0.0, f.30120.0.0, f.30000.0.0, f.30190.0.0, f.30200.0.0, f.30210.0.0,
                               f.30220.0.0, f.30180.0.0)
colnames(bd_pheno3)<-c("FID", "Platelet_count", "Mean_platelet_(thrombocyte)_volume", "Platelet_distribution_width", 
                       "Platelet_crit", "Red_blood_cell_(erythrocyte)_count", "Mean_corpuscular_volume", "Mean_sphered_cell_volume", 
                       "Haematocrit_percentage", "Mean_corpuscular_haemoglobin", "Mean_corpuscular_haemoglobin_concentration", 
                       "Haemoglobin_concentration", "Red_blood_cell_(erythrocyte)_distribution_width", "Reticulocyte_count", 
                       "Mean_reticulocyte_volume", "Reticulocyte_percentage", "Immature_reticulocyte_fraction", 
                       "High_light_scatter_reticulocyte_count", "High_light_scatter_reticulocyte_percentage", 
                       "Nucleated_red_blood_cell_count", "Nucleated_red_blood_cell_percentage", "Monocyte_count", "Neutrophill_count", 
                       "Eosinophill_count", "Basophill_count", "Lymphocyte_count", "White_blood_cell_(leukocyte)_count", 
                       "Monocyte_percentage", "Neutrophill_percentage", "Eosinophill_percentage", "Basophill_percentage", 
                       "Lymphocyte_percentage")
bd_pheno3<-as_tibble(bd_pheno3)
bd_pheno<-bd_pheno %>% left_join(bd_pheno3, by= "FID")






# Add PCs -------------------------------------------------------

cols_PCA<-c("f.eid", sprintf("f.22009.0.%s", 1:20))
PCA<-bd[,cols_PCA]
colnames(PCA)<-c("FID", sprintf("PCA%s", 1:20))

bd_pheno<-bd_pheno %>% inner_join(PCA, by= "FID")



#function 3
myfunction <- function(x){
  ifelse(nchar(x)>3,
         sub("(.{3})(.*)", "\\1.\\2", x),x)
}

# ### Date
# for (i in 0:212){
#   temp <- select(bd,c(paste0("f.41280.0.",i)))
#   names(temp)[1] <- "V1"
#   temp_dat=max(temp$V1, na.rm = TRUE)
#   write.table(temp_dat, file= "/Users/yelab/Downloads/date.txt", col.names = FALSE, append = TRUE, row.names = FALSE, quote = FALSE, sep='\t')
# }
# temp_dat_max <- read.csv("/Users/yelab/Downloads/date.txt",as.is=T, header=F, sep="\t")
# max(temp_dat_max$V1, na.rm = TRUE)
# for (i in 0:212){
#   temp <- select(bd,c(paste0("f.41280.0.",i)))
#   names(temp)[1] <- "V1"
#   temp_dat=min(temp$V1, na.rm = TRUE)
#   write.table(temp_dat, file= "/Users/yelab/Downloads/date2.txt", col.names = FALSE, append = TRUE, row.names = FALSE, quote = FALSE, sep='\t')
# }
# temp_dat_max <- read.csv("/Users/yelab/Downloads/date2.txt",as.is=T, header=F, sep="\t")
# min(temp_dat_max$V1, na.rm = TRUE)
# 
# for (i in 0:46){
#   temp <- select(bd,c(paste0("f.41281.0.",i)))
#   names(temp)[1] <- "V1"
#   temp_dat=max(temp$V1, na.rm = TRUE)
#   write.table(temp_dat, file= "/Users/yelab/Downloads/date1.txt", col.names = FALSE, append = TRUE, row.names = FALSE, quote = FALSE, sep='\t')
# }
# temp_dat_max <- read.csv("/Users/yelab/Downloads/date1.txt",as.is=T, header=F, sep="\t")
# max(temp_dat_max$V1, na.rm = TRUE)
# for (i in 0:46){
#   temp <- select(bd,c(paste0("f.41281.0.",i)))
#   names(temp)[1] <- "V1"
#   temp_dat=min(temp$V1, na.rm = TRUE)
#   write.table(temp_dat, file= "/Users/yelab/Downloads/date11.txt", col.names = FALSE, append = TRUE, row.names = FALSE, quote = FALSE, sep='\t')
# }
# temp_dat_max <- read.csv("/Users/yelab/Downloads/date11.txt",as.is=T, header=F, sep="\t")
# min(temp_dat_max$V1, na.rm = TRUE)


for (i in 0:212){
  temp <- select(bd,c("f.eid",paste0("f.41270.0.",i)))
  colnames(temp)[2] <- 'icd10'
  temp$icd10 <- as.character(temp$icd10)
  temp <- na.omit(temp)
  if(i == 0) { icd10_code = temp}
  else { icd10_code = rbind(icd10_code,temp)}
}
icd10_code <- ddply(icd10_code,.(f.eid,icd10),nrow)
colnames(icd10_code)[3] <- "count"
colnames(icd10_code)[2] <- "code"
colnames(icd10_code)[1] <- "id"
icd10_code$vocabulary_id <- "ICD10CM"
icd10_code <- icd10_code[,c(1,4,2,3)]
id_sex <- select(bd,c("f.eid","f.31.0.0"))
colnames(id_sex) <- c("id","sex")
icd10_code$code <- myfunction(icd10_code$code)
phenotype1 <- createPhenotypes(icd10_code,aggregate.fun=sum, id.sex=bd_pheno$Sex,min.code.count = 1)
#(because for CreatPhenoypte function, the order of colnames is essential)




for (i in 0:46){
  temp <- select(bd,c("f.eid",paste0("f.41271.0.",i)))
  colnames(temp)[2] <- 'icd9'
  temp$icd9 <- as.character(temp$icd9)
  temp <- na.omit(temp)
  if(i == 0) { icd9_code = temp}
  else { icd9_code = rbind(icd9_code,temp)}
}
icd9_code <- ddply(icd9_code,.(f.eid,icd9),nrow)
colnames(icd9_code)[3] <- "count"
colnames(icd9_code)[2] <- "code"
colnames(icd9_code)[1] <- "id"
icd9_code$vocabulary_id <- "ICD9CM"
icd9_code <- icd9_code[,c(1,4,2,3)]


id_sex <- select(bd,c("f.eid","f.31.0.0"))
colnames(id_sex) <- c("id","sex")
icd9_code$code <- myfunction(icd9_code$code)
phenotype2 <- createPhenotypes(icd9_code,aggregate.fun=sum, id.sex=bd_pheno$Sex,min.code.count = 1)
#(Parameter "min.code.count" means the minimum time of one disease report, default is 2)

phenotype3=merge(phenotype1,phenotype2,  all=TRUE)

remove_duplicates <- function(input_df){
  output_df <- input_df[!duplicated(input_df$id),]
  dup_ids=unique(input_df[duplicated(input_df$id),"id"])
  for (i in dup_ids){
    rows_with_dup_ids <- input_df[input_df$id == i,]
    combined_row=apply(rows_with_dup_ids[,names(input_df)[2:ncol(input_df)]],2,any)
    combined_row[is.na(combined_row)]=FALSE
    output_df[output_df$id==i, names(input_df)[2:ncol(input_df)]]<- combined_row
  }
  output_df
}

phenotype4=remove_duplicates(phenotype3)

phenotype4$id <-as.numeric(as.character(phenotype4$id))

bd_pheno<-bd_pheno%>%inner_join(phenotype4, c("FID" = "id"))



########bd_pheno1=bd_pheno

#################################         !!!!!!!!!!!!!!!! Need to change the number!!!!!!!!!
#Map values 0 and 1 to TRUE/FALSE
bd_pheno[,15:ncol(bd_pheno)]<-
  sapply(bd_pheno[,15:ncol(bd_pheno)], mapvalues, c(TRUE, FALSE), c(1, 0))



#QC filtering
bd_pheno<-bd_pheno%>%filter(Genetic_ethnic_grouping == "Caucasian")

#409628
bd_pheno<-bd_pheno%>%filter(Used_in_PCA == "Yes")

#337483
bd_pheno<-bd_pheno%>%filter(is.na(Outliers_for_het_or_missing) | Outliers_for_het_or_missing !="Yes")
#337483
bd_pheno<-bd_pheno%>%filter(is.na(SexchrAneuploidy) | SexchrAneuploidy != "Yes")
#337146
bd_pheno<- bd_pheno%>%filter(is.na(Genetic_kinship) |Genetic_kinship != "Ten or more third-degree relatives identified")
#337146 If Sex does not equal genetic sex, exclude participant
bd_pheno<-bd_pheno[bd_pheno$Sex == bd_pheno$Genetic_Sex,] #remove 378 individuals


#From maximum_set_of_unrelated_individuals.MF.pl output
max_unrelated<-read.table("ukb48818_rel_s488282_output.dat")
max_unrelated<-as.integer(unlist(max_unrelated))
bd_pheno<-bd_pheno%>%filter(!FID %in% max_unrelated)

#298323
bd_pheno<-bd_pheno%>%filter(bd_pheno$Death_Age == -9)

######263231
#Assessment_Non_England=c(11004,11005,11003,11022,11023)
#bd_pheno<-bd_pheno%>%filter(Assessment_centres %ni% Assessment_Non_England )




# bd_pheno$Assessment_centres<-mapvalues(as.character(bd_pheno$Assessment_centres), 
#                                        c(11012, 11021,11011,11008,11003,11024,11020,11005,11004,11018,11010,11016,
#                                          11001, 11017,11009,11013,11002,11007,11014,10003,11006,11022,11023,11025,
#                                          11026,11027,11028), 
#                                        c("Barts","Birmingham","Bristol","Bury","Cardiff","Cheadle_re","Croydon","Edinburgh","Glasgow",
#                                          "Hounslow","Leeds","Liverpool","Manchester","Middlesborough","Newcastle","Nottingham","Oxford",
#                                          "Reading","Sheffield","Stockport_pi","Stoke","Swansea","Wrexham","Cheadle_im",
#                                          "Reading_im","Newcastle_im","Bristol_im"))
# 
#The following `from` values were not present in `x`: 11024, 11025, 11026, 11027, 11028
bd_pheno$Assessment_centres<-mapvalues(as.character(bd_pheno$Assessment_centres),
                                       c(11012, 11021,11011,11008,11003,11020,11005,11004,11018,11010,11016,
                                         11001, 11017,11009,11013,11002,11007,11014,10003,11006,11022,11023),                                      
                                       c('a11012',"a11021","a11011","a11008","a11003","a11020",'a11005',"a11004","a11018","a11010","a11016",
                                         "a11001","a11017","a11009","a11013","a11002","a11007","a11014","a10003","a11006","a11022","a11023"))


save(bd_pheno,file = "Blood_08_05_QC_bd_pheno.RData")
load("Blood_08_05_unQC_bd_pheno.RData")
load("Urine_07_06_QC_bd_pheno.RData")
load("phecode_05_25_bd_pheno.RData")

