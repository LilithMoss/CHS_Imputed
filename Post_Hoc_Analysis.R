snp <- "rs6866110" #SNP of post-hoc interest

#Path to extract data from
output_path <- "/home/pmd-01/chemenya/CHS/PM25/Imputed/Summary/Calls_0.5/"

#Check results for particular method
#Marginal Method
marg <- readRDS(paste0(output_path,"PM25_all_Filtered_GENOTYPE.Y.rds"))
head(marg)
#row.marg <- which(marg$SNP==snp)
res.marg <- marg[marg$SNP==snp,]
#This is the file where chr5 significant SNP is
dat5 <- readRDS(paste0(G_path,"G.chr5.0257"))
j=257
#Extract the SNP
G.snp <- dat5[dat5$V1==snp,]

#Original file
chr5 <- read.table("/home/pmd-01/chemenya/CHS/Split_Imputed/chr5.0257")
rs <- chr5[chr5$V1=="rs6866110",]

rsm <- matrix(rs,ncol=3,byrow=T)
df.rsm <- as.data.frame(rsm)
rsm.small <- df.rsm[df.rsm$V1<0.9 & df.rsm$V2<0.9 & df.rsm$V3<0.9,]
rsm.small0 <- rsm.small[rsm.small$V1>=0.5,]
rsm.small1 <- rsm.small[rsm.small$V1<0.5,]
close <- rsm.small[rsm.small$V1>0.45 & rsm.small$V1<0.55 ,]


t <- t(c(1,2,3,4,5,6))
matrix(t,ncol=3,byrow=T)
exclude <- c(154,556,637,853,869,896,1142,1157,1253,1453,1535,1594,1625,1942,
             1985,1993,2015,2110,2291,2384,2408,2602,2728,2819,2828)

all_path <- "D:/Hispanicity_Imputed_Results/Imputed/Weighted1_1/Marg/"
setwd(all_path)
m <- readRDS("HISP_all_Filtered_GENOTYPE.Y.rds")


