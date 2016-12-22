library(BMA)
library(lmtest)
########################################################################

run.CC <- function(dat){
  if(sum(dat$Y)!=0 & sum(dat$G)!=0 & sum(dat$E)!=0 & sum(dat$G)!=nrow(dat) & nrow(dat)>150){
    m <- glm(Y~E+G+E*G+cov1+cov2+cov3+cov4,data=dat, family="binomial")
    a <- summary(glm(Y~E+G+E*G+cov1+cov2+cov3+cov4,data=dat, family="binomial"))$coef[1,]
    a[1:4] <- c(-99,-99,-99,-99) 
    #if(is.na( m$coefficients[12])){r <- a}
    if(is.na( m$coefficients["E:G"])){r <- a}
    #if(!is.na( m$coefficients[12])){
    if(!is.na( m$coefficients["E:G"])){
      r <- summary(glm(Y~E+G+E*G+cov1+cov2+cov3+cov4,data=dat, family="binomial"))$coef["E:G",] 
    }
    return(r)
    
  } else {
    Null <- as.vector( rep(-99,4) )
    return(Null)
    #if(returnResults) { Null }
  }
}

run.CO <- function(dat){
  if(sum(dat$Y)!=0 & sum(dat$G)!=0 & sum(dat$E)!=0 & sum(dat$G)!=nrow(dat) & nrow(dat)>150){
    cases <- dat[dat$Y==1,]
    m <- glm(E~G+cov1+cov2+cov3+cov4,data=cases, family="binomial")
    a <- summary(glm(E~G+cov1+cov2+cov3+cov4,data=cases, family="binomial"))$coef[1,]
    a[1:4] <- c(-99,-99,-99,-99) 
    if(is.na( m$coefficients["G"])){r <- a}
    if(!is.na( m$coefficients["G"])){
      r <- summary(glm(E~G+cov1+cov2+cov3+cov4,data=cases, family="binomial"))$coef["G",]
    }
    return(r)
  } else {
    Null <- as.vector( rep(-99,4) )
    return(Null)
  }
}

run.GENOTYPE.Y <- function(dat){
  if(sum(dat$Y)!=0 & sum(dat$G)!=0 & sum(dat$E)!=0 & sum(dat$G)!=nrow(dat) & nrow(dat)>150){
    r <- summary(glm(Y~G+cov1+cov2+cov3+cov4,data=dat, family="binomial"))$coef["G",]
    return(r)
  } else {
    Null <- as.vector( rep(-99,4) )
    return(Null)
  }
}

run.GENOTYPE.Y.ADJUST <- function(dat){ #Further adjusts for PM25 - No need for HISPANICITY ANALYSIS
  if(sum(dat$Y)!=0 & sum(dat$G)!=0 & sum(dat$E)!=0 & sum(dat$G)!=nrow(dat) & nrow(dat)>150){
    r <- summary(glm(Y~G+E+cov1+cov2+cov3+cov4,data=dat, family="binomial"))$coef["G",]
    return(r)
  } else {
    Null <- as.vector( rep(-99,4) )
    return(Null)
  }
}

run.GENOTYPE.E <- function(dat){
  if(sum(dat$Y)!=0 & sum(dat$G)!=0 & sum(dat$E)!=0 & sum(dat$G)!=nrow(dat) & nrow(dat)>150){
    r <- summary(glm(E~G+cov1+cov2+cov3+cov4,data=dat, family="binomial"))$coef["G",]
    return(r)
  } else {
    Null <- as.vector( rep(-99,4) )
    return(Null)
    #if(returnResults) { Null }
  }
}

run.DF2 <- function(dat){
  if(sum(dat$Y)!=0 & sum(dat$G)!=0 & sum(dat$E)!=0 & sum(dat$G)!=nrow(dat) & nrow(dat)>150){
    mod.full <- glm(Y~E+G+E*G+cov1+cov2+cov3+cov4,data=dat, family="binomial")
    mod.red <- glm(Y~E+cov1+cov2+cov3+cov4,data=dat, family="binomial")
    r <- lrtest(mod.full,mod.red)[2,5] 
    return(r)
  } else {
    Null <- -99
    return(Null)
    #if(returnResults) { Null }
  }
}

run.CO_DF2 <- function(CC1,VCC,CO1,VCO){ 
  if(CC1!=-99 & CO1!=-99){
    W <- (CC1^2)/(VCC^2) + (CO1^2)/(VCO^2)   
    pval <- pchisq(W, 2, lower.tail=F)
    names(pval)=NULL
    return(pval)
  } else {
    Null <- -99
    return(Null)
    #if(returnResults) { Null }
  }
}

#(Joint-BMA)
run.glib3 <- function(Dat=NULL, models=rbind( c(1,1,1,1,1,1,1), c(1,1,1,0,1,1,1)), pmw=rep(1, nrow(models)), alt.models=c(1)){
  pmw <- pmw/sum(pmw)  
  X3 <- as.data.frame(model.matrix(as.formula(Count~E+G+E:G+Y+E:Y+G:Y+E:G:Y+cov1+Y:cov1+cov2+Y:cov2+cov4+Y:cov4),data=Dat))[,-1]
  r3 <- glib(X3,y=Dat$Count, error="poisson", link = "log", phi=c(1), psi=c(1000),models=models, pmw=pmw, priormean=rep(0, (ncol(X3)+1)),output.postvar=T)
  return(r3)
}

run.BMA.Mult <- function(dat, returnResults=NULL, output.P.only=NULL, DSL=NULL) {
  if(sum(dat$Y)!=0 & sum(dat$G)!=0 & sum(dat$E)!=0 & sum(dat$G)!=nrow(dat) & nrow(dat)>150){
    
    Dat3 <- as.data.frame(ftable(dat$Y,dat$G,dat$E,dat$cov1,dat$cov2,dat$cov4)) 
    names(Dat3) <- c("Y","G","E","cov1","cov2","cov4","Count") 
    X3 <- as.data.frame(model.matrix(as.formula(Count~E+G+E:G+Y+E:Y+G:Y+E:G:Y+cov1+Y:cov1+cov2+Y:cov2+cov4+Y:cov4),data=Dat3))[,-1]
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
    
    # Ebeta <- matrix(c(r3$posterior$mean[10],r3$posterior$mean[15]),ncol=1) #Global expected values of main and interaction effects
    # Gamma1 <- matrix(r3$posterior.bymodel$var[[1]],ncol=16)
    # Gamma2 <- matrix(r3$posterior.bymodel$var[[2]],ncol=15)
    #Translate to not using column number inputs
    Ebeta <- matrix(c(r3$posterior$mean[main],r3$posterior$mean[int]),ncol=1) #Global expected values of main and interaction effects
    Gamma1 <- matrix(r3$posterior.bymodel$var[[1]],ncol=(len+1)) #Simply the covariance matrix for model 1
    Gamma2 <- matrix(r3$posterior.bymodel$var[[2]],ncol=(len))   #Simply the covariance matrix for model 2
    
    # Gamma1.1 <- Gamma1[c(15,24),c(15,24)]
    # Gamma2.2 <- Gamma2[c(14,23),c(14,23)]
    #Translate to not using column number inputs
    Gamma1.1 <- Gamma1[c((main+1),(int+1)),c((main+1),(int+1))] #Covariance of Main and interaction effects from model 1
    Gamma2.2 <- Gamma2[c(main,int),c(main,int)] #Covariance of Main and interaction effects from model 1
    
    # beta.hat1 <- matrix(unlist(r3$posterior.bymodel$mean[1])[c(15,24)],ncol=1)
    # beta.hat2 <- matrix(unlist(r3$posterior.bymodel$mean[2])[c(14,23)],ncol=1)
    #Translate to not using column number inputs
    beta.hat1 <- matrix(unlist(r3$posterior.bymodel$mean[1])[c((main+1),(int+1))],ncol=1)
    beta.hat2 <- matrix(unlist(r3$posterior.bymodel$mean[2])[c(main,int)],ncol=1)
    
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
}

##########################
# HISPANICITY ANALYSIS
#########################
run.CC.hisp <- function(dat){
  if(sum(dat$Y)!=0 & sum(dat$G)!=0 & sum(dat$E)!=0 & sum(dat$G)!=nrow(dat) & nrow(dat)>150){
    m <- glm(Y~E+G+E*G+cov1+cov2+cov3,data=dat, family="binomial")
    a <- summary(glm(Y~E+G+E*G+cov1+cov2+cov3,data=dat, family="binomial"))$coef[1,]
    a[1:4] <- c(-99,-99,-99,-99) 
    #if(is.na( m$coefficients[12])){r <- a}
    if(is.na( m$coefficients["E:G"])){r <- a}
    #if(!is.na( m$coefficients[12])){
    if(!is.na( m$coefficients["E:G"])){
      r <- summary(glm(Y~E+G+E*G+cov1+cov2+cov3,data=dat, family="binomial"))$coef["E:G",] 
    }
    return(r)
    
  } else {
    Null <- as.vector( rep(-99,4) )
    return(Null)
    #if(returnResults) { Null }
  }
}


run.CO.hisp <- function(dat){
  if(sum(dat$Y)!=0 & sum(dat$G)!=0 & sum(dat$E)!=0 & sum(dat$G)!=nrow(dat) & nrow(dat)>150){
    cases <- dat[dat$Y==1,]
    m <- glm(E~G+cov1+cov2+cov3,data=cases, family="binomial")
    a <- summary(glm(E~G+cov1+cov2+cov3,data=cases, family="binomial"))$coef[1,]
    a[1:4] <- c(-99,-99,-99,-99) 
    if(is.na( m$coefficients["G"])){r <- a}
    if(!is.na( m$coefficients["G"])){
      r <- summary(glm(E~G+cov1+cov2+cov3,data=cases, family="binomial"))$coef["G",]
    }
    return(r)
  } else {
    Null <- as.vector( rep(-99,4) )
    return(Null)
  }
}


run.GENOTYPE.Y.hisp <- function(dat){
  if(sum(dat$Y)!=0 & sum(dat$G)!=0 & sum(dat$E)!=0 & sum(dat$G)!=nrow(dat) & nrow(dat)>150){
    r <- summary(glm(Y~G+cov1+cov2+cov3,data=dat, family="binomial"))$coef["G",]
    return(r)
  } else {
    Null <- as.vector( rep(-99,4) )
    return(Null)
  }
}

run.GENOTYPE.E.hisp <- function(dat){
  if(sum(dat$Y)!=0 & sum(dat$G)!=0 & sum(dat$E)!=0 & sum(dat$G)!=nrow(dat) & nrow(dat)>150){
    r <- summary(glm(E~G+cov1+cov2+cov3,data=dat, family="binomial"))$coef["G",]
    return(r)
  } else {
    Null <- as.vector( rep(-99,4) )
    return(Null)
  }
}

run.DF2.hisp <- function(dat){
  if(sum(dat$Y)!=0 & sum(dat$G)!=0 & sum(dat$E)!=0 & sum(dat$G)!=nrow(dat) & nrow(dat)>150){
    mod.full <- glm(Y~E+G+E*G+cov1+cov2+cov3,data=dat, family="binomial")
    mod.red <- glm(Y~E+cov1+cov2+cov3,data=dat, family="binomial")
    r <- lrtest(mod.full,mod.red)[2,5] 
    return(r)
  } else {
    Null <- -99
    return(Null)
    #if(returnResults) { Null }
  }
}

#(Joint-BMA)
run.glib3.hisp <- function(Dat=NULL, models=rbind( c(1,1,1,1,1,1,1), c(1,1,1,0,1,1,1)), pmw=rep(1, nrow(models)), alt.models=c(1)){
  pmw <- pmw/sum(pmw)  
  X3 <- as.data.frame(model.matrix(as.formula(Count~E+G+E:G+Y+E:Y+G:Y+E:G:Y+cov1+Y:cov1+cov2+Y:cov2+cov3+Y:cov3),data=Dat))[,-1]
  r3 <- glib(X3,y=Dat$Count, error="poisson", link = "log", phi=c(1), psi=c(1000),models=models, pmw=pmw, priormean=rep(0, (ncol(X3)+1)),output.postvar=T)
  return(r3)
}

run.BMA.Mult.hisp <- function(dat, returnResults=NULL, output.P.only=NULL, DSL=NULL) {
  if(sum(dat$Y)!=0 & sum(dat$G)!=0 & sum(dat$E)!=0 & sum(dat$G)!=nrow(dat) & nrow(dat)>150){
    
    Dat3 <- as.data.frame(ftable(dat$Y,dat$G,dat$E,dat$cov1,dat$cov2,dat$cov3)) 
    names(Dat3) <- c("Y","G","E","cov1","cov2","cov3","Count") 
    X3 <- as.data.frame(model.matrix(as.formula(Count~E+G+E:G+Y+E:Y+G:Y+E:G:Y+cov1+Y:cov1+cov2+Y:cov2+cov3+Y:cov3),data=Dat3))[,-1]
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
    
    r3 <- run.glib3.hisp(Dat3, models=models3, pmw=pmw, alt.models=alt.models)
    
    # Ebeta <- matrix(c(r3$posterior$mean[10],r3$posterior$mean[15]),ncol=1) #Global expected values of main and interaction effects
    # Gamma1 <- matrix(r3$posterior.bymodel$var[[1]],ncol=16)
    # Gamma2 <- matrix(r3$posterior.bymodel$var[[2]],ncol=15)
    #Translate to not using column number inputs
    Ebeta <- matrix(c(r3$posterior$mean[main],r3$posterior$mean[int]),ncol=1) #Global expected values of main and interaction effects
    Gamma1 <- matrix(r3$posterior.bymodel$var[[1]],ncol=(len+1)) #Simply the covariance matrix for model 1
    Gamma2 <- matrix(r3$posterior.bymodel$var[[2]],ncol=(len))   #Simply the covariance matrix for model 2
    
    
    # Gamma1.1 <- Gamma1[c(15,24),c(15,24)]
    # Gamma2.2 <- Gamma2[c(14,23),c(14,23)]
    #Translate to not using column number inputs
    Gamma1.1 <- Gamma1[c((main+1),(int+1)),c((main+1),(int+1))] #Covariance of Main and interaction effects from model 1
    Gamma2.2 <- Gamma2[c(main,int),c(main,int)] #Covariance of Main and interaction effects from model 1
    
    # beta.hat1 <- matrix(unlist(r3$posterior.bymodel$mean[1])[c(15,24)],ncol=1)
    # beta.hat2 <- matrix(unlist(r3$posterior.bymodel$mean[2])[c(14,23)],ncol=1)
    #Translate to not using column number inputs
    beta.hat1 <- matrix(unlist(r3$posterior.bymodel$mean[1])[c((main+1),(int+1))],ncol=1)
    beta.hat2 <- matrix(unlist(r3$posterior.bymodel$mean[2])[c(main,int)],ncol=1)
    
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
}
