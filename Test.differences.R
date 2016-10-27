#######################################################################
# WHY ARE THERE DIFFERENT NUMBER OF SNPS BASED ON METHOD?
#######################################################################

#######################################################################
# CONCLUSION:
# There is a difference in numbers of non-missing results produced usually
# between CC and all other methods. While the differences in CC and CO are to
# to be expected because they analyze different sets of data, this difference 
# is to be expected and is random. The real reason for observing systematic
# differences between CC/CO and all other methods is that when there is a 
# dataset created by a particular SNP which doesn't have enough G*E=1,
# or enough cases with G*E=1, we get an NA. But since other models (DF2, 
# BMA_DF2, CO_DF2, GENOTYPE.Y, and GENOTYPE.E) don't have these restrictions
# because they don't analyze the interaction term only, they have NA's while
# the other methods produce non-missing values. The only unresolved issue
# is why the results for the missing SNPs are repeated throughout BMA_DF2.
# Later: Check why results are repeated, and whether there are any significant
# SNPs included in these (Garbage) results for SNPs missing in CC.
#######################################################################
output_path <- "D:/PM25_Imputed_Results/Imputed/Summary/"

#Define all methods for summary
methods <- c("BMA_DF2","BMA2","CC","CO","GENOTYPE.Y","GENOTYPE.E","DF2","CO_DF2")
i=methods[5]
j=methods[3]

#Set Chromosome
chr=22

#Read in the data for 1 chromosome
m1 <- readRDS(paste0(output_path,"chr.",chr,".",i,".rds"))
m2 <- readRDS(paste0(output_path,"chr.",chr,".",j,".rds"))

#Compare dimensions between methods for the same chromosome
dim(m1);dim(m2)

#Remove missings from both of the files
#m1
which.name.m1 <- "GENY_P" #Only for GENOTYPE.Y
which.col.m1 <- which(names(m1)==which.name.m1)
m1$P <- as.numeric(levels(m1[,which.col.m1]))[m1[,which.col.m1]]
m1$CHR <- as.numeric(m1$CHR)
m1$BP <- as.numeric(m1$BP)
m1.2 <- m1[m1$P>-99,]
#m2
which.name.m2 <- "CC_P" #Only for GENOTYPE.Y
which.col.m2 <- which(names(m2)==which.name.m2)
m2$P <- as.numeric(levels(m2[,which.col.m2]))[m2[,which.col.m2]]
m2$CHR <- as.numeric(m2$CHR)
m2$BP <- as.numeric(m2$BP)
m2.2 <- m2[m2$P>-99,]

#Compare dimensions after missings are removed
dim(m1.2);dim(m2.2)

#Check which are missing
m1.m <- m1[m1$P==-99,]
m2.m <- m2[m2$P==-99,]
head(m1.m)
head(m2.m)
dim(m1.m);dim(m2.m)

#Check which missing SNPs do not overlap
m1.only <- m1.m[!m1.m$SNP %in% m2.m$SNP,] #missing only in m1
m2.only <- m2.m[!m2.m$SNP %in% m1.m$SNP,] #missing only in m2
both.m <- m1.m[m1.m$SNP %in% m2.m$SNP,]   #missing from both
dim(m1.only);dim(m2.only);dim(both.m)

#Check instances of m2.only SNPs in m1
m1.not.m2 <- m1[m1$SNP %in% m2.only$SNP, ] 
head(m1.not.m2)
head(m1.2)

#Check what the data looks like for the first SNP missing from CC but in GENY
#Read in cov, fam, and phen data
phen <- read.table("C:/Users/Lilith Moss/Documents/MATH/RESEARCH/BMA_Project/#Real_Data/GEN2-Code/Imputed Data/chs.pheno",header=T)
fam <- read.table("C:/Users/Lilith Moss/Documents/MATH/RESEARCH/BMA_Project/#Real_Data/GEN2-Code/Imputed Data/chs3000.fam",header=F)
cov <- read.table("C:/Users/Lilith Moss/Documents/MATH/RESEARCH/BMA_Project/#Real_Data/GEN2-Code/Imputed Data/chs.cov",header=T)
names(fam) <- c("FID","IID","PID","MID","SEX","PHEN")
#Match all data up to fam order
cov.matched <- cov[match(fam$IID,cov$IID),]
phen.matched <- phen[match(fam$IID,phen$IID),]
#Create useful Covariate file with only relevant variables
Cov <- data.frame( cbind(cov.matched$male,cov.matched$afr,cov.matched$natam,
                         cov.matched$asian,cov.matched$ses) )
names(Cov) <- c("male","afr","natam","asian","ses")
#fh = family history
#hw = hispanic whites
#Assign Vectors
#E <- cov$hw #Environmental Factor = hispanic whites
E <- ifelse(cov.matched$pm25 <= median(cov.matched$pm25),0,1) #Environmental Factor = dichotomized pm25
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

#Test with the first missing SNP from CC rs149862772
#Extract SNP from the data
library(data.table)
map <- read.table("D:/Genomic_Data/chr22.map")
chr.row <- which(map$V2=="rs149862772") #17th row
d <- read.table("D:/Genomic_Data/chr22.dose",nrow=1, skip=(chr.row-1))
D <- d[,4:ncol(d)]

#Put the SNP in dominant format
#Define dominant G function
dom <- function(D){
  G <- vector(mode="numeric", length=3000)
  for(k in 0:2999){
    if(D[(3*k)+1]>=0.9){
      G[k+1]=0
     } else if(D[(3*k)+2]>0.9){
      G[k+1]=1  
      } else if(D[(3*k)+3]>0.9){
      G[k]=1  
      } else {
      G[k+1]=-9
      }
  }
  return(G)
} 
#Call dominant function
G <- dom(D) #Only one non-zero

#Run analysis
G[G==-9] <- NA
dat <- as.data.frame( cbind(Y,G,E,cov1,cov2,cov3,cov4,ses2,ses3,ses4,ses5) )
Dat <- dat[complete.cases(dat),]
source("Real.functions_CapPP.R")
#Assign Prior Probabilities
I_ge=0.5;I_g=0.5;I_gxe=0.5 #Set prior probabilities
BMA_DF2.MVW <- matrix(run.BMA.Mult(Dat,T,F,F),ncol=7) #BMA_DF2.MVW will output pval, PPM1(CC),PPM2(CO),PPSNP,PPSNPxE,GlimEst.SNP,GlimEst.SNPxE
BMA2 <- matrix(run.BMA2(Dat,T,F,F),ncol=4)
CC <- matrix(run.CC(Dat),ncol=4)
CO <- matrix(run.CO(Dat),ncol=4)      
GENOTYPE.Y <- matrix(run.GENOTYPE.Y(Dat),ncol=4)
GENOTYPE.E <- matrix(run.GENOTYPE.E(Dat),ncol=4)
DF2 <- matrix(run.DF2(Dat),ncol=1)
CO_DF2 <- matrix(t(t(run.CO_DF2(GENOTYPE.Y[1],GENOTYPE.Y[2],CO[1],CO[2]))),ncol=1)

#Figure out which method actually produces the right missing value
#CC method
    m <- glm(Y~E+G+E*G+cov1+cov2+cov3+cov4+ses2+ses3+ses4+ses5,data=Dat, family="binomial")
#CO method
    cases=Dat[Y==1,]
    m <- glm(E~G+cov1+cov2+cov3+cov4+ses2+ses3+ses4+ses5,data=cases, family="binomial")
#GENOTYPE.Y
    m <- glm(Y~G+cov1+cov2+cov3+cov4+ses2+ses3+ses4+ses5,data=Dat, family="binomial")
EG <- Dat$E*Dat$G

#So what is the problem with BMA_DF2
run.BMA.Mult <- function(dat, returnResults=NULL, output.P.only=NULL, DSL=NULL) {
    
    Dat3 <- as.data.frame(ftable(dat$Y,dat$G,dat$E,dat$cov1,dat$cov2,dat$cov3,dat$cov4,dat$ses2,dat$ses3,dat$ses4,dat$ses5)) #Without Covariates
    names(Dat3) <- c("Y","G","E","cov1","cov2","cov3","cov4","ses2","ses3","ses4","ses5","Count") #Without Covariates
    X3 <- as.data.frame(model.matrix(as.formula(Count~E+G+E:G+Y+E:Y+G:Y+E:G:Y+cov1+Y:cov1+cov2+Y:cov2+cov3+Y:cov3+cov4+Y:cov4+ses2+Y:ses2+ses3+Y:ses3+ses4+Y:ses4+ses5+Y:ses5),data=Dat3))[,-1]
    len <- ncol(X3)
    co <- which(colnames(X3)=="E1:G1")
    main <- which(colnames(X3)=="G1:Y1")
    int <- which(colnames(X3)=="E1:G1:Y1")
    models3 <- rbind( c(rep(1,len)), c(rep(1,(co-1)),0,rep(1,(len-co))) ) # These indicators correspond to c(E1,G,Y1,E1:G,E1:Y1,G:Y1,E1:G:Y1)
    mat0 <- matrix(models3[,c(co,main,int)],ncol=3)  #E1:G1    #G1:Y1   #E1:G1:Y1
    ge <- ifelse(mat0[,1]==1,I_ge,1-I_ge)
    g <- ifelse(mat0[,2]==1,I_g,1-I_g)
    gxe <- ifelse(mat0[,3]==1,I_gxe,1-I_gxe)
    mod_prior <- matrix(cbind(ge,g,gxe),ncol=3) #No need to condition since p(M)=1 because all 8 are used
    pmw <- mod_prior[,1]*mod_prior[,2]*mod_prior[,3] #Set the model priors; don't have to add up to 1
    alt.models<- c(1:2)
    pmw <- pmw/sum(pmw)  
    
    r3 <- run.glib3(Dat3, models=models3, pmw=pmw, alt.models=alt.models)
    
    Ebeta <- matrix(c(r3$posterior$mean[14],r3$posterior$mean[23]),ncol=1) #Global expected values of main and interaction effects
    Gamma1 <- matrix(r3$posterior.bymodel$var[[1]],ncol=24)
    Gamma2 <- matrix(r3$posterior.bymodel$var[[2]],ncol=23)
    
    Gamma1.1 <- Gamma1[c(15,24),c(15,24)]
    Gamma2.2 <- Gamma2[c(14,23),c(14,23)]
    
    beta.hat1 <- matrix(unlist(r3$posterior.bymodel$mean[1])[c(15,24)],ncol=1)
    beta.hat2 <- matrix(unlist(r3$posterior.bymodel$mean[2])[c(14,23)],ncol=1)
    
    pr1 <- r3$bf$postprob[1]
    pr2 <- r3$bf$postprob[2]
    
    Var1 <- (Gamma1.1+(beta.hat1%*%t(beta.hat1)))*pr1
    Var2 <- (Gamma2.2+(beta.hat2%*%t(beta.hat2)))*pr2
    Gamma <- Var1+Var2 - Ebeta%*%t(Ebeta)
    inv.Gamma <- solve(Gamma)
    W <- t(Ebeta)%*%inv.Gamma%*%Ebeta
    pval <- pchisq(W,2,lower.tail=F)
    
    #Posterior probabilities of models 1(CC) and 2(CO)
    PPM1 <- r3$bf$postprob[1]
    PPM2 <- r3$bf$postprob[2]
    
    #Posterior GLIM effect estimates for the SNP and the interaction for each model
    GlimEst.SNP1 <- r3$glim.est$coef[[1]]["G1"] #CC
    GlimEst.SNP2 <- r3$glim.est$coef[[2]]["G1"] #CO
    
    GlimEst.SNPxE1 <- r3$glim.est$coef[[1]]["E1:G1:Y1"] #CC
    GlimEst.SNPxE2 <- r3$glim.est$coef[[2]]["E1:G1:Y1"] #CO
    
    #Combine pval, PPM1(CC),PPM2(CO),GlimEst.SNP(1&2),GlimEst.SNPxE(1&2)
    output <- cbind(pval,PPM1,PPM2,GlimEst.SNP1,GlimEst.SNP2,GlimEst.SNPxE1,GlimEst.SNPxE2)
    
    #if(returnResults) { pval }
    if(returnResults) { output }
  } else {
    #Null <- -99
    Null <- as.vector( rep(-99,7) )
    if(returnResults) { Null }
  }


run.glib3 <- function(Dat=NULL, models=rbind( c(1,1,1,1,1,1,1), c(1,1,1,0,1,1,1)), pmw=rep(1, nrow(models)), alt.models=c(1)){
  pmw <- pmw/sum(pmw)  
  X3 <- as.data.frame(model.matrix(as.formula(Count~E+G+E:G+Y+E:Y+G:Y+E:G:Y+cov1+Y:cov1+cov2+Y:cov2+cov3+Y:cov3+cov4+Y:cov4+ses2+Y:ses2+ses3+Y:ses3+ses4+Y:ses4+ses5+Y:ses5),data=Dat))[,-1]
  r3 <- glib(X3,y=Dat$Count, error="poisson", link = "log", phi=c(1), psi=c(1000),models=models, pmw=pmw, priormean=rep(0, (ncol(X3)+1)),output.postvar=T)
  return(r3)
}





