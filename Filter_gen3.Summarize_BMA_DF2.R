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

#Define all methods for summary
methods <- c("BMA_DF2","BMA2","CC","CO","GENOTYPE.Y","GENOTYPE.E","DF2","CO_DF2")
i=methods[1]

#Read results from former save
# system.time( results <- readRDS(paste0(output_path,"PM25_all_chr",i,".rds")) )
# test <- readRDS(paste0(output_path,"chr.22.BMA_DF2.rds"))
# results <- test

# all.files <- list.files(segmented_path,pattern=i)
# #Composite.chr will be a composite of all chromosomes for method i
# composite.chr <- do.call(rbind,lapply(1:22,function(chr){ #For each chromosome
#   #chr.files is a list of all files in chromosome chr and method i
#   chr.files <- list.files(segmented_path,pattern=paste0(i,".chr",chr,".0"))
#   #composite.file will be a composite of one chromosome for method i
#   composite.file <- do.call(rbind,lapply(1:length(chr.files),function(y){ #For each chromosome file
#     rds.file <- readRDS(paste0(segmented_path,chr.files[y]))  #Read file
#   }))
#   map <- read.table(paste0(map_path,"chr",chr,".map"))
#   final.chr <- cbind(map,composite.file)
#   colnames(final.chr)[1:4] <- c("CHR","SNP.map","V","BP")
#   saveRDS(final.chr,paste0(output_path,"chr.",chr,".",i,".rds"))
#   final.chr
# }))

results <- do.call(rbind,lapply(1:22,function(chr){
  dat <- readRDS(paste0(output_path,"chr.",chr,".",i,".rds"))
}))

#Read in the info files
#info.hw <- as.data.frame(fread(paste0(info_path,"hw.info.freq")))[,-1]
#info.nhw <- as.data.frame(fread(paste0(info_path,"nhw.info.freq")))[,-1]

#names(info.hw) <- c("SNP","Freq.hw","Info.hw")
#names(info.nhw) <- c("SNP","Freq.nhw","Info.nhw")

#info <- as.data.frame(fread(paste0(info_path,"combined.info.freq")))
info <- as.data.frame(readRDS(paste0(info_path,"combined.info.freq.rds")))
names(info) <- c("SNP","Freq.nhw","Info.nhw","Freq.hw","Info.hw","Freq")

#mdat.1 <- merge(results,info.nhw,by="SNP",all.x=T)
#mdat.2 <- merge(mdat.1, info.hw,by="SNP",all.x=T)
mdat.2 <- merge(results,info,by="SNP",all.x=T)

#Write info matched data
saveRDS("mdat.2",paste0(output_path,"PM25_all_Info_Freq_",i,".rds"))

#Filter on info >=0.7
mdat.filtered <- mdat.2[mdat.2$Info.hw>=0.7 & mdat.2$Info.nhw>=0.7,]
low <- mdat.2[mdat.2$Info.hw<0.7 | mdat.2$Info.nhw<0.7,]   
range(low$Info.hw)
range(low$Info.nhw)
head(low,n=50)

#Filter on MAF <0.5
mdat.filtered.final <- mdat.filtered[mdat.filtered$Freq>0.05,]

#Structure data types for plotting and lambda calculation
mdat.filtered.final$P <- mdat.filtered.final$BMA_P
which.name <- "BMA_P" #Only for BMA_DF2
which.col <- which(names(mdat.filtered.final)==which.name)
mdat.filtered.final$P <- as.numeric(levels(mdat.filtered.final[,which.col]))[mdat.filtered.final[,which.col]]
mdat.filtered.final$CHR <- as.numeric(mdat.filtered.final$CHR)
mdat.filtered.final$BP <- as.numeric(mdat.filtered.final$BP)
mdat.filtered2 <- mdat.filtered.final[mdat.filtered.final$P>-99,]

#Calculate Lambda (Genetic Inflation Factor)
chisq <- qchisq(1-mdat.filtered2$P,1)
lambda = median(chisq)/qchisq(0.5,1) #CC = 0.3901901 #BMA_DF2 = 0.6731412
lambda2 <- round(lambda,3)

   #Manhattan only
   jpeg(paste0(output_path,"PM25_MAN_",i,"_Filtered.jpg"),width=800)
   manhattan(mdat.filtered2,main=paste0(i,", ",nrow(mdat.filtered2)," SNPs, E = PM2.5"))
   dev.off()
   
   #QQ Only
   jpeg(paste0(output_path,"PM25_QQ_",i,"_Filtered.jpg"),width=800)
   qq(mdat.filtered2$P,main=paste0(i,", ",nrow(mdat.filtered2)," SNPs, E = PM2.5"))
   mtext(lambda,side=1,line=-3.5,cex=1) #Add lamda value to plot
   dev.off()
   
   #Order by significance and save results
   mdat.filtered3 <- mdat.filtered2[order(mdat.filtered2$P),]
   saveRDS(mdat.filtered3,paste0(output_path,"PM25_all_Filtered_",i,".rds"))
   top1000 <- mdat.filtered3[1:1000,]
   write.csv(top1000,paste0(output_path,"PM25_1000_Filtered_",i,".csv"),row.names=F)
   
   # #Zoom in on chromosome 22
   # chr22 <- mdat.filtered2[mdat.filtered2$CHR==22,]
   # jpeg(paste0(output_path,"PM25_MAN_Chr22_",i,"_Filtered.jpg"),width=800)
   # manhattan(chr22,main=paste0(i,", Chr22, ",nrow(mdat.filtered2)," SNPs, E = PM2.5"))
   # dev.off()
   # 
   # #Zoom in to significant region
   # jpeg(paste0(output_path,"PM25_MAN__rs11703393_",i,"_Filtered.jpg"),width=800)
   # manhattan(chr22,xlim=c(44.48,44.54),highlight=c("rs11703393"),main = "Chr22 rs11703393 Region")
   # dev.off()
 
   
   
   
   