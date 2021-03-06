###########################################
# Summarizes results of PM2.5 GxE GWIS
# Read in segmented rds files
###########################################
library(qqman) 
segmented_path <- "/home/pmd-01/chemenya/CHS/PM25/Imputed/"
map_path <- "/auto/pmd-01/chemenya/CHS/txtDosage/"
output_path <- "/home/pmd-01/chemenya/CHS/PM25/Imputed/Summary/"

#Define all methods for summary
methods <- c("BMA_DF2","BMA2","CC","CO","GENOTYPE.Y","GENOTYPE.E","DF2","CO_DF2")
i=methods[8]
#Loop through each method
#for(i in methods){
  all.files <- list.files(segmented_path,pattern=i)
   #Composite.chr will be a composite of all chromosomes for method i 
   #composite.chr <- do.call(rbind,lapply(1:22,function(chr){ #For each chromosome
   composite.chr <- do.call(rbind,lapply(21:22,function(chr){ #For each chromosome
     #chr.files is a list of all files in chromosome chr and method i 
     chr.files <- list.files(segmented_path,pattern=paste0(i,".chr",chr,".0"))
        #composite.file will be a composite of one chromosome for method i
        #composite.file <- do.call(rbind,lapply(1:length(chr.files),function(y){ #For each chromosome file
         composite.file <- do.call(rbind,lapply(1:length(chr.files),function(y){ #For each chromosome file    
             rds.file <- readRDS(paste0(segmented_path,chr.files[y]))  #Read file
         }))
      map <- read.table(paste0(map_path,"chr",chr,".map"))
      final.chr <- cbind(map,composite.file)
      colnames(final.chr)[1:4] <- c("CHR","SNP.map","V","BP")
      saveRDS(final.chr,paste0(output_path,"chr.",chr,".",i,".rds"))
      final.chr
    }))
   #QQ and Manhattan Plots
   composite.chr$P <- as.numeric(levels(composite.chr[,10]))[composite.chr[,10]]
   composite.chr$CHR <- as.numeric(composite.chr$CHR)
   composite.chr$BP <- as.numeric(composite.chr$BP)
   composite.chr$SNP <- as.character(composite.chr$SNP)
   composite.chr2 <- composite.chr[composite.chr$P>-99,]
   
   #Plot and save as jpeg
   jpeg(paste0(output_path,"PM25_",i,".jpg"),width=1000)
   par(mfcol=c(1,2))
   qq(composite.chr2$P,main=paste0(i,", ",nrow(composite.chr2)," SNPs, E = PM2.5"))
   manhattan(composite.chr2,main=paste0(i,", ",nrow(composite.chr2)," SNPs, E = PM2.5"))
   dev.off()
   
   #Order by significance
   composite.chr2 <- composite.chr2[order(composite.chr2$P),]
   saveRDS(composite.chr2,paste0(output_path,"PM25_all_chr",i,".rds"))
   write.csv(composite.chr2,paste0(output_path,"PM25_all_chr",i,".csv"),row.names=F)
#} For-loop

   
   
   
   
   
#Test

rds.dir <- "/home/pmd-01/chemenya/CHS/PM25/Imputed/Summary/"
composite.chr <- do.call(rbind,lapply(21:22,function(chr){
  chr.file <- readRDS(paste0(rds.dir,"chr.",chr,".CO_DF2.rds"))
}))
jpeg(paste0("PM25_",i,".jpg"),width=1000)
par(mfcol=c(1,2))
qq(composite.chr2$P,main=paste0(i,", ",nrow(composite.chr2)," SNPs, E = PM2.5"))
manhattan(composite.chr2,main=paste0(i,", ",nrow(composite.chr2)," SNPs, E = PM2.5"))
dev.off()

test <- composite.chr2[1:100,]
jpeg(paste0("PM25_test.jpg"),width=1000)
#par(mfcol=c(1,2))
#qq(composite.chr2$P,main=paste0(i,", ",nrow(composite.chr2)," SNPs, E = PM2.5"))
manhattan(test)
dev.off()

test <- readRDS("test.rds")
manhattan(test)

test$CHR[test$CHR==21] <- 1
test$CHR[test$CHR==22] <- 2
