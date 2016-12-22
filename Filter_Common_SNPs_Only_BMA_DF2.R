library(qqman) 
library(data.table)
#Local
# output_path <- "D:/PM25_Imputed_Results/Imputed/Summary/"
# info_path <- "C:/Users/Lilith Moss/Documents/MATH/RESEARCH/BMA_Project/#Real_Data/Imputed_Info_Scores/"

#HPC
# output_path <- "/home/pmd-01/chemenya/CHS/PM25/Imputed/Summary/BMA_DF2/"
output_path <- "/home/pmd-01/chemenya/CHS/PM25/Imputed/Summary/Calls_0.5/"
info_path <- "/auto/pmd-01/chemenya/CHS/Imputed_Info_Scores/"
segmented_path <- "/home/pmd-01/chemenya/CHS/PM25/Imputed/Calls_0.5/"
map_path <- "/auto/pmd-01/chemenya/CHS/txtDosage/"
common_snp_path <- "/auto/pmd-01/chemenya/CHS/Split_Imputed_Scripts/posthoc.code/"


#Define all methods for summary
methods <- c("BMA_DF2","CC","CO","GENOTYPE.Y","DF2","CO_DF2")
i=methods[1]

#Read back SNP list if saved
mdat.2 <- readRDS(paste0(output_path,"PM25_all_Info_Freq_",i,".rds"))

mdat.filtered3 <- readRDS(paste0(output_path,"PM25_all_Filtered_",i,".rds"))
mdat.filtered3 <- readRDS(pasteo(output_path,"PM25_all_Filtered_BMA_DF2.rds"))