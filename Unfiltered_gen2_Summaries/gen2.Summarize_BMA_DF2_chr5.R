###########################################
# Summarizes results of PM2.5 GxE GWIS
# Read in segmented rds files
###########################################
library(qqman) 
# segmented_path <- "/home/pmd-01/chemenya/CHS/PM25/Imputed/"
# map_path <- "/auto/pmd-01/chemenya/CHS/txtDosage/"
# output_path <- "/home/pmd-01/chemenya/CHS/PM25/Imputed/Summary/"

#Local
output_path <- "D:/PM25_Imputed_Results/Imputed/Summary/"

#Define all methods for summary
methods <- c("BMA_DF2","BMA2","CC","CO","GENOTYPE.Y","GENOTYPE.E","DF2","CO_DF2")
i=methods[1]

   # composite.chr <- do.call(rbind,lapply(1:22,function(chr){
   #   dat <- readRDS(paste0(output_path,"chr.",chr,".",i,".rds"))
   # }))

if(i=="BMA_DF2"){
  tot <- matrix(,ncol = 16, nrow = 0) #If BMA_DF2
} else {
  tot <- matrix(,ncol = 13 ,nrow = 0) #If CC,CO, GenY
} 
   system.time(for(chr in 5:5){
     #assign(paste("chr",chr,sep=""), readRDS(paste0(output_path,"chr.",chr,".",i,".rds")))
     chrx <- readRDS(paste0(output_path,"chr.",chr,".",i,".rds"))
     tot <- rbind(tot,chrx)
     rm(chrx)
     print(object.size(tot))
     print(chr)
   })
   
   #Bind all chromosomes together
   # composite.chr <- rbind(chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,
   #                        chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22) 
   
   #QQ and Manhattan Plot Set-Up
   #which.name <- paste0(i,"_P")
   which.name <- "BMA_P" #Only for BMA_DF2
   #which.name <- "GENY_P" #Only for GENOTYPE.Y
   composite.chr <- tot
   which.col <- which(names(composite.chr)==which.name)
   #composite.chr$P <- as.numeric(levels(composite.chr[,10]))[composite.chr[,10]]
   composite.chr$P <- as.numeric(levels(composite.chr[,which.col]))[composite.chr[,which.col]]
   composite.chr$CHR <- as.numeric(composite.chr$CHR)
   composite.chr$BP <- as.numeric(composite.chr$BP)
   composite.chr2 <- composite.chr[composite.chr$P>-99,]
   
   #Calculate Lambda ########################################
   # library(gap)
   # r <- gcontrol2(composite.chr$P)
   # lam <- round(r$lambda,3)
   # lam.lab <- bquote(lambda ~ "=" ~ .(lam))
   
   composite.chr <- readRDS("D:/PM25_Imputed_Results/Imputed/Summary/BMA_DF2/PM25_all_chrBMA_DF2.rds")
   library(gap)
   r <- gcontrol2(composite.chr$P)
   lambda <- r$lambda
   #Lambda manually
   chisq <- qchisq(1-composite.chr$P,1)
   lambda = median(chisq)/qchisq(0.5,1) #CC = 0.3901901 #BMA_DF2 = 0.6731412
   lambda
   ########################################################
   
   
   #Manhattan only
   jpeg(paste0(output_path,"PM25_MAN_",i,".jpg"),width=800)
   manhattan(composite.chr2,main=paste0(i,", ",nrow(composite.chr2)," SNPs, E = PM2.5"))
   dev.off()
   
   #QQ Only
   jpeg(paste0(output_path,"PM25_QQ_",i,".jpg"),width=800)
   qq(composite.chr2$P,main=paste0(i,", ",nrow(composite.chr2)," SNPs, E = PM2.5"))
   #mtext(lam.lab,side=1,line=-3.5,cex=1) #Add lamda value to plot
   dev.off()
   
   #Both QQ and Manhattan
   # jpeg(paste0(output_path,"PM25_Both_",i,".jpg"),width=1000)
   # par(mfcol=c(1,2))
   # qq(composite.chr2$P,ylim=0.8,main=paste0(i,", ",nrow(composite.chr2)," SNPs, E = PM2.5"))
   # mtext(lam.lab,side=1,line=-3.5,cex=1) #Add lamda value to plot
   # manhattan(composite.chr2,main=paste0(i,", ",nrow(composite.chr2)," SNPs, E = PM2.5"))
   # dev.off()
   
   #Order by significance and save results
   composite.chr3 <- composite.chr2[order(composite.chr2$P),]
   saveRDS(composite.chr3,paste0(output_path,"PM25_all_chr",i,".rds"))
   top1000 <- composite.chr3[1:1000,]
   write.csv(top1000,paste0(output_path,"PM25_1000_chr",i,".csv"),row.names=F)
   
   #Zoom in on chromosome 22
   chr22 <- composite.chr2[composite.chr2$CHR==22,]
   jpeg(paste0(output_path,"PM25_MAN_Chr22_",i,".jpg"),width=800)
   manhattan(chr22,main=paste0(i,", Chr22, ",nrow(composite.chr2)," SNPs, E = PM2.5"))
   dev.off()
   
   #Zoom in to significant region
   jpeg(paste0(output_path,"PM25_MAN__rs11703393_",i,".jpg"),width=800)
   manhattan(chr22,xlim=c(44.48,44.54),highlight=c("rs11703393"),main = "Chr22 rs11703393 Region")
   dev.off()
    
   #Looking at chr5
   #Zoom in on chromosome 5
   chr5 <- composite.chr2
   chr5 <- chr5[order(chr5$P),]
   jpeg(paste0(output_path,"PM25_MAN_Chr5_",i,".jpg"),width=800)
   manhattan(chr5,main=paste0(i,", Chr5, ",nrow(composite.chr2)," SNPs, E = PM2.5"))
   dev.off()
   
   #Zoom in to significant region
   jpeg(paste0(output_path,"PM25_MAN_chr5_Zoom_",i,".jpg"),width=800)
   manhattan(chr5,xlim=c(148.5,149.5),main = "Chr5 Significant Region")
   dev.off()
   