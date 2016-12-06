
#Read the relevant file
G_path <- "/home/pmd-01/chemenya/CHS/Split_Imputed_Results/" #Path where rds files are
dat <- readRDS(paste0(G_path,"G.chr8.0146"))
#Check that all the significant SNPs are in this file
sig <- c("rs4738376",
         "rs11774483",
         "rs12235061",
         "rs2383914",
         "rs4237005",
         "rs4738366")


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
#E <- cov$hw
E <- cov.matched$hw
E2 <- ifelse(cov.matched$pm25 <= median(cov.matched$pm25),0,1) #Environmental Factor = dichotomized pm25
Y <- phen.matched$asthma #Phenotype = asthma  

#Categorize continuous variables into dichotomies
Cov$eu <- 1-(Cov$afr+Cov$natam+Cov$asian)
Cov$afr <- ifelse(Cov$afr<median(Cov$afr),0,1)
Cov$natam <- ifelse(Cov$natam<median(Cov$natam),0,1)
Cov$asian <- ifelse(Cov$asian<median(Cov$asian),0,1)
Cov$eu <- ifelse(Cov$eu<median(Cov$eu),0,1)
#Cov$ses1 <- If all are zero
Cov$ses2 <- ifelse(Cov$ses==2,1,0)
Cov$ses3 <- ifelse(Cov$ses==3,1,0)
Cov$ses4 <- ifelse(Cov$ses==4,1,0)
Cov$ses5 <- ifelse(Cov$ses==5,1,0)
cov1 <- Cov$male
cov2 <- Cov$afr
cov3 <- Cov$natam
cov4 <- Cov$asian
cov5 <- Cov$eu
ses2 <- Cov$ses2
ses3 <- Cov$ses3
ses4 <- Cov$ses4
ses5 <- Cov$ses5

#Get relevant G
g.rs <- sig[1]
r <- which(dat[,1]==g.rs)
G <- t( dat[r,4:3003] )
G[G==-9] <- NA
dat.1 <- as.data.frame( cbind(Y,G,E,E2,cov1,cov2,cov3,cov4,cov5,ses2,ses3,ses4,ses5) )
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






