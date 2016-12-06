###########################################
# Summarizes results of PM2.5 GxE GWIS
# Read in segmented rds files
###########################################
library(qqman) 
library(data.table)
#Local
 output_path <- "D:/Hispanicity_Imputed_Results/Imputed/Summary/BMA_DF2/"
 info_path <- "C:/Users/Lilith Moss/Documents/MATH/RESEARCH/BMA_Project/#Real_Data/Imputed_Info_Scores/"

#HPC
output_path <- "/home/pmd-01/chemenya/CHS/PM25/Imputed/Summary/BMA_DF2/"
info_path <- "/auto/pmd-01/chemenya/CHS/Imputed_Info_Scores/"

#Define all methods for summary
methods <- c("BMA_DF2","BMA2","CC","CO","GENOTYPE.Y","GENOTYPE.E","DF2","CO_DF2")
i=methods[1]

#Read results from former save
#system.time( results <- readRDS(paste0(output_path,"PM25_all_chr",i,".rds")) )
system.time( results <- readRDS(paste0(output_path,"HISP_all_Filtered_",i,".rds")) )

# #Info filter - Don't have to do this for HISPANICITY Files,
# #HISPANICITY results are already filtered.
# info <- as.data.frame(readRDS(paste0(info_path,"combined.info.freq.rds")))
# names(info) <- c("SNP","Freq.nhw","Info.nhw","Freq.hw","Info.hw","Freq")
# 
# # mdat.1 <- merge(results,info.nhw,by="SNP",all.x=T)
# # mdat.2 <- merge(mdat.1, info.hw,by="SNP",all.x=T)
# 
# mdat.2 <- merge(results,info,by="SNP",all.x=T)
# 
# #Write info matched data
# saveRDS("mdat.2",paste0(output_path,"HISP_all_Info_Freq_",i,".rds"))
# 
# #Filter on >=0.7
# mdat.filtered <- mdat.2[mdat.2$Info.hw>=0.7 & mdat.2$Info.nhw>=0.7,]
# low <- mdat.2[mdat.2$Info.hw<0.7 | mdat.2$Info.nhw<0.7,]   
# range(low$Info.hw)
# range(low$Info.nhw)
# head(low,n=50)
# 
# #Filter on MAF <0.5
# mdat.filtered.final <- mdat.filtered[mdat.filtered$Freq>0.05,]
# 
# #Structure data types for plotting and lambda calculation
# #mdat.filtered$P <- mdat.filtered$BMA_P
# which.name <- "GENY_P" #Only for BMA_DF2
# which.col <- which(names(mdat.filtered.final)==which.name)
# mdat.filtered.final$P <- as.numeric(levels(mdat.filtered.final[,which.col]))[mdat.filtered.final[,which.col]]
# mdat.filtered.final$CHR <- as.numeric(mdat.filtered.final$CHR)
# mdat.filtered.final$BP <- as.numeric(mdat.filtered.final$BP)
# mdat.filtered2 <- mdat.filtered.final[mdat.filtered.final$P>-99,]

#Structure data for plotting
mdat.filtered2 <- results



#Calculate Lambda (Genetic Inflation Factor)
chisq <- qchisq(1-mdat.filtered2$P,1)
lambda = median(chisq)/qchisq(0.5,1) #CC = 0.3901901 #BMA_DF2 = 0.6731412
lambda2 <- round(lambda,3)

   #Manhattan only
   jpeg(paste0(output_path,"HISP_MAN_",i,"_Filtered2.jpg"),width=800)
   manhattan(mdat.filtered2,main=paste0(i,", ",nrow(mdat.filtered2)," SNPs, E = HISP"))
   dev.off()
   
   #QQ Only
   jpeg(paste0(output_path,"HISP_QQ_",i,"_Filtered2.jpg"),width=800)
   qq(mdat.filtered2$P,main=paste0(i,", ",nrow(mdat.filtered2)," SNPs, E = HISP"))
   mtext(lambda2,side=1,line=-3.5,cex=1) #Add lamda value to plot
   dev.off()
   
   #Order by significance and save results
   mdat.filtered3 <- mdat.filtered2[order(mdat.filtered2$P),]
   top1000 <- mdat.filtered3[1:1000,]
   write.csv(top1000,paste0(output_path,"HISP_1000_Filtered2_",i,".csv"),row.names=F)
   saveRDS(mdat.filtered3,paste0(output_path,"HISP_all_Filtered2_",i,".rds"))
   
   #Zoom in on chromosome 8
   chr8 <- mdat.filtered2[mdat.filtered2$CHR==8,]
   jpeg(paste0(output_path,"HISP_MAN_Chr8_",i,"_Filtered2.jpg"),width=800)
   manhattan(chr8,main=paste0(i,", Chr8, ",nrow(mdat.filtered2)," SNPs, E = HISP"))
   dev.off()
   
   #Zoom in to significant region
   jpeg(paste0(output_path,"HISP_MAN__rs4738376_",i,"_Filtered2.jpg"),width=800)
   manhattan(chr8,xlim=c(74,75),highlight=c("rs4738376"),main = "Chr8 rs4738376 Region")
   dev.off()
 
   
   