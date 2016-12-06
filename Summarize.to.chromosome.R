###########################################
# Summarizes results of PM2.5 GxE GWIS
# Read in segmented rds files
###########################################
library(qqman) 
library(data.table)
#Local
# output_path <- "D:/PM25_Imputed_Results/Imputed/Summary/"
# info_path <- "C:/Users/Lilith Moss/Documents/MATH/RESEARCH/BMA_Project/#Real_Data/Imputed_Info_Scores/"

#HPC
# output_path <- "/home/pmd-01/chemenya/CHS/PM25/Imputed/Summary/BMA_DF2/"
output_path <- "/home/pmd-01/chemenya/CHS/asthma/Imputed/Summary/"
info_path <- "/auto/pmd-01/chemenya/CHS/Imputed_Info_Scores/"
segmented_path <- "/home/pmd-01/chemenya/CHS/asthma/Imputed/"
map_path <- "/auto/pmd-01/chemenya/CHS/txtDosage/"


for(i in 1:6)
#Define all methods for summary
methods <- c("BMA_DF2","CC","CO","GENOTYPE.Y","DF2","COgetw_DF2")
i=methods[1]

#Read results from former save
# system.time( results <- readRDS(paste0(output_path,"PM25_all_chr",i,".rds")) )
# test <- readRDS(paste0(output_path,"chr.22.BMA_DF2.rds"))
# results <- test

all.files <- list.files(segmented_path,pattern=i)
#Composite.chr will be a composite of all chromosomes for method i
composite.chr <- do.call(rbind,lapply(1:22,function(chr){ #For each chromosome
  #chr.files is a list of all files in chromosome chr and method i
  chr.files <- list.files(segmented_path,pattern=paste0(i,".chr",chr,".0"))
  #composite.file will be a composite of one chromosome for method i
  composite.file <- do.call(rbind,lapply(1:length(chr.files),function(y){ #For each chromosome file
    rds.file <- readRDS(paste0(segmented_path,chr.files[y]))  #Read file
  }))
  map <- read.table(paste0(map_path,"chr",chr,".map"))
  final.chr <- cbind(map,composite.file)
  colnames(final.chr)[1:4] <- c("CHR","SNP.map","V","BP")
  saveRDS(final.chr,paste0(output_path,"chr.",chr,".",i,".rds"))
  final.chr
}))
