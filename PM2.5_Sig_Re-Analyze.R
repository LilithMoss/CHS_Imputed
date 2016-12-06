###################################################################
# Re-analyze the SNPs from the significant chromosome region
# from BMA_DF2 GxHispanicity interaction findin (most sig: rs4738376)
###################################################################
#PM2.5 Significant SNPs
# Chr	SNP
# 22	rs11705133
# 22	rs67355122
# 22	rs72619560
# 22	rs11703393
# 5	  rs6866110
# 22	rs62227671
# 22	rs2267611
# 19	rs12980282
# 19	rs2335626
# 19	rs2288890
pm2.5 <- c("rs11705133","rs11703393","rs6866110","rs12980282") #chr: 22,22,5,19 
hisp <- c("rs4672623","rs73351224","chr8:3242613:D","rs75265348") #chr: 2,14,8,1
search_path <- "/home/pmd-01/chemenya/CHS/Split_Imputed_Results/"
G_path <- search_path
I_ge = 0.5
I_g = 0.5 
I_gxe = 0.5



#PM2.5
#Chr22: rs11705133
#Chr5: rs6866110
#Chr19: rs12980282

# chr=19
# snp=pm2.5[4]
# #Look for the gen file in RDS format
# for(i in 0:9){
#   #t <- readRDS(paste0("G.chr2.000",i))
#   #a <- "rs4672623" %in% t[,1]
#   t <- readRDS(paste0(search_path,"G.chr",chr,".000",i))
#   a <- snp %in% t[,1]
#   print(a);print(i)
# }
# 
# for(i in 10:99){
#   #t <- readRDS(paste0("G.chr2.00",i))
#   t <- readRDS(paste0(search_path,"G.chr",chr,".00",i))
#   a <- snp %in% t[,1]
#   print(a);print(i)
# }
# 
# for(i in 100:350){
#   #t <- readRDS(paste0("G.chr2.0",i))
#   t <- readRDS(paste0(search_path,"G.chr",chr,".0",i))
#   a <- snp %in% t[,1]
#   print(a);print(i)
# } 
# 
# G_path <- search_path
# dat22 <- readRDS(paste0(G_path,"G.chr22.0045"))
# sig <- c(pm2.5[1],pm2.5[2])
# sig %in% dat22[,1] #All in this file from this region?
# dat5 <- readRDS(paste0(G_path,"G.chr5.0257"))
# sig <- c(pm2.5[3])
# sig %in% dat5[,1] #All in this file from this region?
# dat19 <- readRDS(paste0(G_path,"G.chr19.0009"))
# sig <- c(pm2.5[4])
# sig %in% dat19[,1] #All in this file from this region?
# 
# #HISPANICITY
# #Chr2: rs4672623
# #Chr14: rs73351224
# #Chr8: chr8:3242613:D
# #Chr1: rs75265348
# 
# chr=1
# snp=hisp[4]
# #Look for the gen file in RDS format
# rd2 <- do.call(rbind,lapply(0:9,function(i){
#   t <- readRDS(paste0(search_path,"G.chr",chr,".000",i))
#   a <- snp %in% t[,1]
#   print(a);print(i)
#   cbind(a,i)
# }))
# 
# rd2 <- do.call(rbind,lapply(10:99,function(i){
#   t <- readRDS(paste0(search_path,"G.chr",chr,".00",i))
#   a <- snp %in% t[,1]
#   print(a);print(i)
#   cbind(a,i)
# }))
# 
# rd2 <- do.call(rbind,lapply(100:400,function(i){
#   t <- readRDS(paste0(search_path,"G.chr",chr,".0",i))
#   a <- snp %in% t[,1]
#   print(a);print(i)
#   cbind(a,i)
# }))
# 
# dat2 <- readRDS(paste0(G_path,"G.chr2.0351"))
# sig <- c(hisp[1])
# sig %in% dat2[,1] #All in this file from this region?
# 
# dat14 <- readRDS(paste0(G_path,"G.chr14.0132"))
# sig <- c(hisp[2])
# sig %in% dat14[,1] #All in this file from this region?
# 
# dat8 <- readRDS(paste0(G_path,"G.chr8.0009"))
# sig <- c(hisp[3])
# sig %in% dat8[,1] #All in this file from this region?
# 
# dat1 <- readRDS(paste0(G_path,"G.chr1.0269"))
# sig <- c(hisp[4])
# sig %in% dat1[,1] #All in this file from this region?

#Significant SNPs
pm2.5.chr <- c(22,22,5,19)
pm2.5.file <- c("0045","0045","0257","0009")
hisp.file <- c("0351","0132","0009","0269")
hisp.chr <- c(2,14,8,1)

pm2.5.sig <- data.frame( cbind(pm2.5.chr,pm2.5,pm2.5.file) )
hisp.sig <- data.frame( cbind(hisp.chr,hisp,hisp.file) )

########################
#Re-Analyze for PM2.5
########################
#For each of the SNPs, do an analysis to make sure I get the same thing
#Load in phen and cov
mainDir <- "/auto/pmd-01/chemenya/CHS/"
#source(paste0(mainDir,"Real.functions_CapPP.R"))
epi_path <- "/auto/pmd-01/chemenya/CHS/CHS_GxEScan/Data/" #Path from which to read cov/phen files
fam_path <- "/home/pmd-01/chemenya/CHS/txtDosage/" #Path from which to read whole chromosome data, fam, map
dimension_path <- "/home/pmd-01/chemenya/CHS/txtDosage/dimensions/"

#Read in Epi files
phen <- read.table(paste0(epi_path,"chs.pheno"),header=T) #Read in phenotype #3000 samples
fam <- read.table(paste0(fam_path,"chs3000.fam"),header=F) #Read in fam file #3000 Samples
names(fam) <- c("FID","IID","PID","MID","SEX","PHEN")
cov <- read.table(paste0(epi_path,"chs.cov"),header=T) #Read in Covariates #3000 Samples

#Match all data up to fam order
cov.matched <- cov[match(fam$IID,cov$IID),]
phen.matched <- phen[match(fam$IID,phen$IID),]

#Create useful Covariate file with only relevant variables
#With Hispanicity
Cov <- data.frame( cbind(cov.matched$male,cov.matched$afr,cov.matched$natam,
                         cov.matched$asian,cov.matched$ses,cov.matched$hw) )
names(Cov) <- c("male","afr","natam","asian","ses","hw")

#Assign Vectors
# E <- cov.matched$hw #Environmental Factor = hispanic whites
# E <- cov$hw
E <- ifelse(cov.matched$pm25 <= median(cov.matched$pm25),0,1) #Environmental Factor = dichotomized pm25
Y <- phen.matched$asthma #Phenotype = asthma  

#Categorize continuous variables into dichotomies
cov1 <- Cov$male
Cov$natam1 <- ifelse(Cov$natam<=0.05,1,0) #European
Cov$natam2 <- ifelse(Cov$natam<=0.5 & Cov$natam>0.05,1,0) #HW
Cov$natam3 <- ifelse(Cov$natam>0.5,1,0) 
cov2 <- Cov$natam2
cov3 <- Cov$natam3
cov4 <- Cov$hw

tab.pm2.5 <- lapply(1:4,function(i){
  dat <- readRDS(paste0(G_path,"G.chr",pm2.5.sig[i,1],".",pm2.5.sig[i,3]))
  g.rs <- as.character( pm2.5.sig[i,2] )
  r <- which(dat[,1]==g.rs)
  G <- t( dat[r,4:3003] )
  G[G==-9] <- NA
  #dat.1 <- as.data.frame( cbind(Y,G,E,cov1,cov2,cov3,cov4,ses2,ses3,ses4,ses5) )
  dat.1 <- as.data.frame( cbind(Y,G,E,cov1,cov2,cov3,cov4) )
  colnames(dat.1)[2] <- "G" 
  Dat <- dat.1[complete.cases(dat.1),]
  
  #Get Cell Counts
  cell.counts <- as.data.frame(ftable(Dat$Y,Dat$G,Dat$E,Dat$cov1,Dat$cov2,Dat$cov3,Dat$cov4))
  names(cell.counts) <- c("Y","G","E","Male","NA1","NA2","HW","Freq")
  write.csv(cell.counts,paste0("cell.counts_PM2.5_",g.rs,"csv"))
  #Analysis
  #CC
  CC <- run.CC(Dat)
  #CO
  CO <- run.CO(Dat)
  #MARG
  Marg <- run.GENOTYPE.Y(Dat)
  #BMA
  BMA <- run.BMA.Mult(Dat,T,F,F)
  #All
  list(CC,CO,Marg,BMA)
})

cc.tab <- do.call(rbind,unlist( lapply(1:4,function(x) tab.pm2.5[[x]][1]),recursive=F))  
co.tab <- do.call(rbind,unlist( lapply(1:4,function(x) tab.pm2.5[[x]][2]),recursive=F))  
marg.tab <- do.call(rbind,unlist( lapply(1:4,function(x) tab.pm2.5[[x]][3]),recursive=F))  
bma.tab <- do.call(rbind,unlist( lapply(1:4,function(x) tab.pm2.5[[x]][4]),recursive=F))  

write.csv(cc.tab,"cc_PM2.5.csv",row.names=F)
write.csv(co.tab,"cc_PM2.5.csv",row.names=F)
write.csv(marg.tab,"cc_PM2.5.csv",row.names=F)
write.csv(bma.tab,"cc_PM2.5.csv",row.names=F)

#HISPANICITY
##############################
#Re-Analyze for HISPANICITY
##############################
mainDir <- "/auto/pmd-01/chemenya/CHS/"
#source(paste0(mainDir,"Real.functions_CapPP.R"))
epi_path <- "/auto/pmd-01/chemenya/CHS/CHS_GxEScan/Data/" #Path from which to read cov/phen files
fam_path <- "/home/pmd-01/chemenya/CHS/txtDosage/" #Path from which to read whole chromosome data, fam, map
dimension_path <- "/home/pmd-01/chemenya/CHS/txtDosage/dimensions/"

#Read in Epi files
phen <- read.table(paste0(epi_path,"chs.pheno"),header=T) #Read in phenotype #3000 samples
fam <- read.table(paste0(fam_path,"chs3000.fam"),header=F) #Read in fam file #3000 Samples
names(fam) <- c("FID","IID","PID","MID","SEX","PHEN")
cov <- read.table(paste0(epi_path,"chs.cov"),header=T) #Read in Covariates #3000 Samples

#Match all data up to fam order
cov.matched <- cov[match(fam$IID,cov$IID),]
phen.matched <- phen[match(fam$IID,phen$IID),]

#Create useful Covariate file with only relevant variables
#With Hispanicity
Cov <- data.frame( cbind(cov.matched$male,cov.matched$afr,cov.matched$natam,
                         cov.matched$asian,cov.matched$ses,cov.matched$hw) )
names(Cov) <- c("male","afr","natam","asian","ses","hw")

#Assign Vectors
E <- cov.matched$hw #Environmental Factor = hispanic whites
#E <- ifelse(cov.matched$pm25 <= median(cov.matched$pm25),0,1) #Environmental Factor = dichotomized pm25
Y <- phen.matched$asthma #Phenotype = asthma  

#Categorize continuous variables into dichotomies
cov1 <- Cov$male
Cov$natam1 <- ifelse(Cov$natam<=0.05,1,0) #European
Cov$natam2 <- ifelse(Cov$natam<=0.5 & Cov$natam>0.05,1,0) #HW
Cov$natam3 <- ifelse(Cov$natam>0.5,1,0) 
cov2 <- Cov$natam2
cov3 <- Cov$natam3
cov4 <- rep(0,nrow(Cov))

tab.hisp <- lapply(1:4,function(i){
  dat <- readRDS(paste0(G_path,"G.chr",hisp.sig[i,1],".",hisp.sig[i,3]))
  g.rs <- as.character( hisp.sig[i,2] )
  r <- which(dat[,1]==g.rs)
  G <- t( dat[r,4:3003] )
  G[G==-9] <- NA
  #dat.1 <- as.data.frame( cbind(Y,G,E,cov1,cov2,cov3,cov4,ses2,ses3,ses4,ses5) )
  dat.1 <- as.data.frame( cbind(Y,G,E,cov1,cov2,cov3,cov4) )
  colnames(dat.1)[2] <- "G" 
  Dat <- dat.1[complete.cases(dat.1),]
  
  #Get Cell Counts
  cell.counts <- as.data.frame(ftable(Dat$Y,Dat$G,Dat$E,Dat$cov1,Dat$cov2,Dat$cov3))
  names(cell.counts) <- c("Y","G","E","Male","NA1","NA2","Freq")
  write.csv(cell.counts,paste0("cell.counts_HISP_",g.rs,"csv"))
  #Analysis
  #CC
  CC <- run.CC(Dat)
  #CO
  CO <- run.CO(Dat)
  #MARG
  Marg <- run.GENOTYPE.Y(Dat)
  #BMA
  BMA <- run.BMA.Mult(Dat,T,F,F)
  #All
  list(CC,CO,Marg,BMA)
})

cc.tab <- do.call(rbind,unlist( lapply(1:4,function(x) tab.pm2.5[[x]][1]),recursive=F))  
co.tab <- do.call(rbind,unlist( lapply(1:4,function(x) tab.pm2.5[[x]][1]),recursive=F))  
marg.tab <- do.call(rbind,unlist( lapply(1:4,function(x) tab.pm2.5[[x]][1]),recursive=F))  
bma.tab <- do.call(rbind,unlist( lapply(1:4,function(x) tab.pm2.5[[x]][1]),recursive=F))  

write.csv(cc.tab,"cc_PM2.5.csv",row.names=F)
write.csv(co.tab,"cc_PM2.5.csv",row.names=F)
write.csv(marg.tab,"cc_PM2.5.csv",row.names=F)
write.csv(bma.tab,"cc_PM2.5.csv",row.names=F)










#Run Analyses
#Assign Prior Probabilities
I_ge=0.5;I_g=0.5;I_gxe=0.5 #Set prior probabilities
BMA_DF2.MVW <- matrix(run.BMA.Mult(Dat,T,F,F),ncol=7) #BMA_DF2.MVW will output pval, PPM1(CC),PPM2(CO),PPSNP,PPSNPxE,GlimEst.SNP,GlimEst.SNPxE
colnames(BMA_DF2.MVW) <- c("BMA_P","BMA_PPM1","BMA_PPM2","BMA_GlimEst.SNP1","BMA_GlimEst.SNP2",  #BMA_DF2 (1)
                        "BMA_GlimEst.SNPxE1","BMA_GlimEst.SNPxE2")
#Case-Control check
CC <- matrix(run.CC(Dat),ncol=4)
#r.hw <- summary(glm(Y~E+G+E*G+cov1+cov2+cov3+cov4+ses2+ses3+ses4+ses5,data=dat, family="binomial"))$coef["E:G",]
r <- summary(glm(Y~E+G+E*G+cov1+cov2+cov3+cov4+ses2+ses3+ses4+ses5,data=Dat, family="binomial"))
r.hw <- glm(Y~G+cov1+cov2+cov3+cov4+ses2+ses3+ses4+ses5,subset=E==1,data=dat, family="binomial")
r.w <- glm(Y~G+cov1+cov2+cov3+cov4+ses2+ses3+ses4+ses5,subset=E==0,data=dat, family="binomial")
effects.hw <- exp( c(coef(r.hw)["G"],confint(r.hw)[2,]) )
effects.w <- exp( c(coef(r.w)["G"],confint(r.w)[2,]) )
t <- table(Dat$G,Dat$E)

#Check other method results
chr.8.BMA_DF2.rds

test <- readRDS("chr.8.CC.rds")
r <- which(test[,2]=="rs4738376")
test[r,]
new.cc <- summary(glm(Y~))


chr.8.CO.rds
chr.8.GENOTYPE.Y.rds






