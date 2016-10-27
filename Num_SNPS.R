#Count the number of scripts per chromosome for proliferation purposes
#Cluster Paths
dimension_path <- "/home/pmd-01/chemenya/CHS/txtDosage/dimensions/"
output_path <- "/home/pmd-01/chemenya/CHS/Split_Imputed_Results/"
dosage_path <- "/home/pmd-01/chemenya/CHS/Split_Imputed/"

count <- do.call(rbind,lapply(1:22,function(x){
  read.table(paste0(dimension_path,"chr",x,".row"))
}))

number.of.snps <- sum(count[,1])

#Number of SNPs = 23,207,726