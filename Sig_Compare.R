###################################################################
# Re-analyze the SNPs from the significant chromosome region
# from BMA_DF2 GxHispanicity interaction findin (most sig: rs4738376)
###################################################################
#PM2.5
#Which ones are significant?
sig22 <- c("rs11705133","rs67355122","rs72619560","rs11703393","rs62227671","rs2267611","rs62227671")
sig5 <- "rs6866110"
sig19 <- c("rs12980282","rs2335626","rs2288890")

#Marginal Compare
marg22 <- readRDS("chr.22.GENOTYPE.Y.rds")
marg5 <- readRDS("chr.5.GENOTYPE.Y.rds")
marg19 <- readRDS("chr.19.GENOTYPE.Y.rds")
marg22.res <- marg22[which(marg22$SNP.map %in% sig22),]
marg5.res <- marg5[which(marg5$SNP.map %in% sig5),]
marg19.res <- marg19[which(marg19$SNP.map %in% sig19),]
write.csv(marg22.res,"PM25.marg22.comp.csv",row.names=F)
write.csv(marg5.res,"PM25.marg5.comp.csv",row.names=F)
write.csv(marg19.res,"PM25.marg19.comp.csv",row.names=F)

#Case Control Compare chr.21.CC.rds
cc22 <- readRDS("chr.22.CC.rds")
cc5 <- readRDS("chr.5.CC.rds")
cc19 <- readRDS("chr.19.CC.rds")
cc22.res <- cc22[which(cc22$SNP.map %in% sig22),]
cc5.res <- cc5[which(cc5$SNP.map %in% sig5),]
cc19.res <- cc19[which(cc19$SNP.map %in% sig19),]
write.csv(cc22.res,"PM25.cc22.comp.csv",row.names=F)
write.csv(cc5.res,"PM25.cc5.comp.csv",row.names=F)
write.csv(cc19.res,"PM25.cc19.comp.csv",row.names=F)

#Case Only Compare
co22 <- readRDS("chr.22.CO.rds")
co5 <- readRDS("chr.5.CO.rds")
co19 <- readRDS("chr.19.CO.rds")
co22.res <- co22[which(co22$SNP.map %in% sig22),]
co5.res <- co5[which(co5$SNP.map %in% sig5),]
co19.res <- co19[which(co19$SNP.map %in% sig19),]
write.csv(co22.res,"PM25.co22.comp.csv",row.names=F)
write.csv(co5.res,"PM25.co5.comp.csv",row.names=F)
write.csv(co19.res,"PM25.co19.comp.csv",row.names=F)


#HISPANICITY
#Which ones are significant?
sig2 <-	"rs4672623"
sig8 <- "chr8:3242613:D"
sig14 <- c("rs73351224","rs79326859","rs73351229",
           "rs8021015","rs73351228","rs73351227")
sig1 <- c("rs75265348","rs2878713")

#Marginal Compare
marg2 <- readRDS("chr.2.GENOTYPE.Y.rds")
marg8 <- readRDS("chr.8.GENOTYPE.Y.rds")
marg14 <- readRDS("chr.14.GENOTYPE.Y.rds")
marg1 <- readRDS("chr.1.GENOTYPE.Y.rds")
marg2.res <- marg2[which(marg2$SNP.map %in% sig2),]
marg8.res <- marg8[which(marg8$SNP.map %in% sig8),]
marg14.res <- marg14[which(marg14$SNP.map %in% sig14),]
marg1.res <- marg1[which(marg1$SNP.map %in% sig1),]
write.csv(marg2.res,"HISP.marg2.comp.csv",row.names=F)
write.csv(marg8.res,"HISP.marg8.comp.csv",row.names=F)
write.csv(marg14.res,"HISP.marg14.comp.csv",row.names=F)
write.csv(marg1.res,"HISP.marg1.comp.csv",row.names=F)

#Case Control Compare chr.21.CC.rds
cc2 <- readRDS("chr.2.CC.rds")
cc8 <- readRDS("chr.8.CC.rds")
cc14 <- readRDS("chr.14.CC.rds")
cc1 <- readRDS("chr.1.CC.rds")
cc2.res <- cc2[which(cc2$SNP.map %in% sig2),]
cc8.res <- cc8[which(cc8$SNP.map %in% sig8),]
cc14.res <- cc14[which(cc14$SNP.map %in% sig14),]
cc1.res <- cc1[which(cc1$SNP.map %in% sig1),]
write.csv(cc2.res,"HISP.cc2.comp.csv",row.names=F)
write.csv(cc8.res,"HISP.cc8.comp.csv",row.names=F)
write.csv(cc14.res,"HISP.cc14.comp.csv",row.names=F)
write.csv(cc1.res,"HISP.cc1.comp.csv",row.names=F)

#Case Only Compare
co2 <- readRDS("chr.2.CO.rds")
co8 <- readRDS("chr.8.CO.rds")
co14 <- readRDS("chr.14.CO.rds")
co1 <- readRDS("chr.1.CO.rds")
co2.res <- co2[which(co2$SNP.map %in% sig2),]
co8.res <- co8[which(co8$SNP.map %in% sig8),]
co14.res <- co14[which(co14$SNP.map %in% sig14),]
co1.res <- co1[which(co1$SNP.map %in% sig1),]
write.csv(co2.res,"HISP.co2.comp.csv",row.names=F)
write.csv(co8.res,"HISP.co8.comp.csv",row.names=F)
write.csv(co14.res,"HISP.co14.comp.csv",row.names=F)
write.csv(co1.res,"HISP.co1.comp.csv",row.names=F)
















##############
# OLD
##############

#Read the relevant file
G_path <- "/home/pmd-01/chemenya/CHS/Split_Imputed_Results/" #Path where rds files are
dat <- readRDS(paste0(G_path,"G.chr8.0146"))
#Check that all the significant SNPs are in this file



sig %in% dat[,1] #Yes, all are in here

#For each of the SNPs, do an analysis to make sure I get the same thing
#Load in phen and cov
mainDir <- "/auto/pmd-01/chemenya/CHS/"
source(paste0(mainDir,"Real.functions_CapPP.R"))
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
Cov <- data.frame( cbind(cov.matched$male,cov.matched$afr,cov.matched$natam,
                         cov.matched$asian,cov.matched$ses) )
names(Cov) <- c("male","afr","natam","asian","ses")
#Assign Vectors
E <- cov.matched$hw #Environmental Factor = hispanic whites
E <- cov$hw
#E <- ifelse(cov.matched$pm25 <= median(cov.matched$pm25),0,1) #Environmental Factor = dichotomized pm25
Y <- phen.matched$asthma #Phenotype = asthma  

#Categorize continuous variables into dichotomies
Cov$afr <- ifelse(Cov$afr<median(Cov$afr),0,1)
Cov$natam <- ifelse(Cov$natam<median(Cov$natam),0,1)
Cov$asian <- ifelse(Cov$asian<median(Cov$asian),0,1)
#Cov$ses1 <- If all are zero
Cov$ses2 <- ifelse(Cov$ses==2,1,0)
Cov$ses3 <- ifelse(Cov$ses==3,1,0)
Cov$ses4 <- ifelse(Cov$ses==4,1,0)
Cov$ses5 <- ifelse(Cov$ses==5,1,0)
cov1 <- Cov$male
cov2 <- Cov$afr
cov3 <- Cov$natam
cov4 <- Cov$asian
ses2 <- Cov$ses2
ses3 <- Cov$ses3
ses4 <- Cov$ses4
ses5 <- Cov$ses5

#Get relevant G
g.rs <- sig[1]
r <- which(dat[,1]==g.rs)
G <- t( dat[r,4:3003] )
G[G==-9] <- NA
dat.1 <- as.data.frame( cbind(Y,G,E,cov1,cov2,cov3,cov4,ses2,ses3,ses4,ses5) )
Dat <- dat.1[complete.cases(dat.1),]
colnames(Dat)[2] <- "G" 

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






