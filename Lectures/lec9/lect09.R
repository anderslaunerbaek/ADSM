rm(list=ls())
library(nlme)
#######################################################################
# Example: Ramus bone length
ramus<-read.table("ramus.txt",header=T)

## Plot data
plot(ramus[ ,-1],pch=19,col=as.numeric(ramus[ ,1]))
I <- as.numeric(ramus[ ,1])
for(i in 1:max(I)){
    lines(ramus[I==i,2],ramus[I==i,3],col=i,lty=2)
}

## Modelling
ramus$agered <- ramus$age - 8.75

## The model
ramus.lme1 <- lme(ramus ~ agered, 
                  random=~agered|Boy,data=ramus)
ramus.lme1

## Fixed effects
fixed.effects(ramus.lme1)

## Estimated randomeffects
random.effects(ramus.lme1) # not actual parameters

## Variance-covariance of random effects
VarCorr(ramus.lme1)

(var.ramus<-as.numeric(VarCorr(ramus.lme1)[1:2,1]))
(rho <- as.numeric(VarCorr(ramus.lme1)[2,3]))

## As in the slides
sig.Psi <- diag(var.ramus)
sig.Psi[1,2] <- sig.Psi[2,1] <-  rho * sqrt(prod(var.ramus))
sig.Psi

fitted(ramus.lme1)

##################################################
## Models for each boy
plot(ramus[ ,c(4,3)],pch=19,col=as.numeric(ramus[ ,1]))
lines(ramus[1:4,4],fitted(ramus.lme1)[1:4],col=1,lty=1)
lines(ramus[5:8,4],fitted(ramus.lme1)[5:8],col=2,lty=1)
lines(ramus[9:12,4],fitted(ramus.lme1)[9:12],col=3,lty=1)
lines(ramus[13:16,4],fitted(ramus.lme1)[13:16],col=4,lty=1)
lines(ramus[17:20,4],fitted(ramus.lme1)[17:20],col=5,lty=1)
##################################################

## Plot of intercept v
pairs(ramus.lme1,pch=19)

residuals(ramus.lme1,type="p")
residuals(ramus.lme1,type="n")


#####################################################
## Random effects
#####################################################

ranef(ramus.lme1)
sig.Psi
(var.ramus<-as.numeric(VarCorr(ramus.lme1)[1:2,1]))

sig.Psi <- diag(var.ramus)
sig.Psi[1,2] <- sig.Psi[2,1] <-  rho * sqrt(prod(var.ramus))
sig.Psi

Z <- matrix(0,ncol=10,nrow=20)
Z[1:4,1:2]<-Z[5:8,3:4]<-Z[9:12,5:6]<-
    Z[13:16,7:8]<-Z[17:20,9:10]<-
    cbind(1,ramus[1:4,4])
Psi <- kronecker(diag(5),sig.Psi)
Sigma <- diag(20)*ramus.lme1$sigma^2
X <- cbind(1,ramus$agered)

## Estimation of random effects (by matrix calculations)
solve(t(Z)%*%solve(Sigma)%*%Z+solve(Psi))%*%t(Z)%*%solve(Sigma)%*%
  (ramus[ ,3] - X%*%ramus.lme1$coefficients$fixed)
ranef(ramus.lme1)

## Information of random effect (see p. 183)
Iu <- t(Z)%*%solve(Sigma)%*%Z + solve(Psi)
sqrt(diag(solve(Iu)))


##################################################
##
##################################################

##################################################
# Rat example 
rats.tmp<-read.table("rats.csv",sep=";",header=T)
rats<-c()
for(i in 1:dim(rats.tmp)[1]){
  rats<-rbind(rats,cbind(rats.tmp[i,1],rats.tmp[i,2],
                         1:10,t(rats.tmp[i,3:12])))
}

rats<-data.frame(treatm=rats[ ,1],cage=rats[ ,2],month=rats[ ,3],
                 lnc=log(rats[ ,4]))
library(lattice)

head(rats,3)

# Plot the data
 plot(rats$lnc~rats$month, type='n')
 by(rats,rats$cage,
    function(X){
      lines(X$month,X$lnc, col=X$treatm[1])
    }
 )->out
 legend("topright", bty='n', legend=1:3, lty='solid', col=1:3)
################################################################

rats$treatm<-factor(rats$treatm)
rats$month<-factor(rats$month)

 
##################################################
## Simple lme
 library(nlme)
 fit.mm<-lme(lnc~month+treatm+month:treatm, random = ~1|cage, data=rats)
 
fit.mm
anova(fit.mm)
##################################################

fit.cs<-gls(lnc~month+treatm+month:treatm,
            correlation=corCompSymm(form=~1|cage),
            data=rats)

logLik(fit.cs)
logLik(fit.mm)
## exactly the same

#################################################
# Likelihhod (for comparing model with anova)
 fit.cs<-gls(lnc~month+treatm+month:treatm,
             correlation=corCompSymm(form=~1|cage),
             data=rats, method="ML")
 logLik(fit.cs)
 fit.mm<-lme(lnc~month+treatm+month:treatm, 
             random = ~1|cage, 
             data=rats, method="ML")
 logLik(fit.mm)
##################################################

##################################################
## Gaussian correlation structure
fit.gau <- lme(lnc~month+treatm+month:treatm,
               random=~1|cage,
               correlation=corGaus(form=~as.numeric(month)|cage,nugget=TRUE),
               data=rats)
 ##################################################


##################################################
fit.gau
##################################################

##################################################
 nu.sq<-0.1404056^2
 sigma.sq<-0.2171559^2*0.2186743
 tau.sq<-0.2171559^2-sigma.sq
 rho.sq<-2.3863954
 c(nu.sq=nu.sq, sigma.sq=sigma.sq, tau.sq=tau.sq, rho.sq=rho.sq)

##################################################

####################################################
## Comparing models 
##################################################
# Note: Wrong analysis.....
fit.mm <- lme(lnc~month+treatm+month:treatm, 
              random = ~1|cage, data=rats)
fit.mmR <- lme(lnc~month+treatm, 
              random = ~1|cage, data=rats)
anova(fit.mm,fit.mmR) # Note that you need to use ML.. and we get a warning.
########################################################


 fit.id <- lm(lnc~month+treatm+month:treatm, data=rats) # WRONG independent model 
 fit.mm <- lme(lnc~month+treatm+month:treatm, 
              random = ~1|cage, 
               data=rats, method="ML")
 fit.gau<- lme(lnc~month+treatm+month:treatm,
               random=~1|cage,
               correlation=corGaus(form=~as.numeric(month)|cage,nugget=TRUE),
               data=rats, method="ML")
 anova(fit.gau,fit.mm,fit.id)
##################################################

 ## Non nested models
 fit.exp<- lme(lnc~month+treatm+month:treatm, random=~1|cage,
               correlation=corExp(form=~as.numeric(month)|cage,nugget=TRUE),
               data=rats, method="ML")
 AIC(fit.exp)
 AIC(fit.gau)
 anova(fit.gau,fit.exp,test=FALSE) 
 
 ## So we prefer exp struct.

 fit.exp
 fit.exp2<- lme(lnc~month+treatm+month:treatm, random=~1|cage,
               correlation=corExp(form=~as.numeric(month)|cage,nugget=FALSE),
               data=rats, method="ML")
 
 anova(fit.exp,fit.exp2)
 
##################################################
 # Analysis using REML (for reporting results)
 fit.gau<- lme(lnc~month+treatm+month:treatm, random=~1|cage,
               correlation=corGaus(form=~as.numeric(month)|cage,nugget=TRUE),
               data=rats, method="REML")
 fit.exp<- lme(lnc~month+treatm+month:treatm, random=~1|cage,
               correlation=corExp(form=~as.numeric(month)|cage,nugget=FALSE),
               data=rats, method="REML")

 
fit.exp 
 
#graph(width=7,height=7,type="eps",file="gau")
 plot(Variogram(fit.gau), main='Gaussian')
#dev.off()
#graph(width=7,height=7,type="eps",file="exp")
 plot(Variogram(fit.exp), main='Exponential')
#dev.off()

##################################################




# Calculation of variogram...
N <- numeric(9); V <- numeric(9)
for(i in 1:max(as.numeric(rats$treatm))){
    for(j in 1:max(as.numeric(rats$cage))){
        I <- as.numeric(rats$treatm)==i & as.numeric(rats$cage)==j
        for(k in 1:9){
            V[k] <- V[k] + sum(diff(residuals(fit.gau)[I],lag=k)^2)
            N[k] <- N[k] + length(diff(residuals(fit.gau)[I],lag=k))
        }
    }
}
V
N
V/(2*N)

(Variogram(fit.gau,resType="response")[ ,1]-V/(2*N)) ## These are the same

###########################################################
## End
# Som plots for slides
x<-seq(-3,3,by=0.01)
gaussCor<-function(x){
  exp(-x*x)
}
y<-0.1+gaussCor(x)

#source("~/tsmodels/graph.R")
#graph(width=14,height=7,type="eps",file="xspgau")
plot(x,y,type="l",ylim=c(0,1.2),axes=F,
     xlab=expression(paste("Distance ",t[i[2]],"-",t[i[1]])),ylab="Covariance")
lines(range(x)*2,c(0.1,0.1),lty=2)
lines(range(x)*2,max(y)*c(1,1),lty=2)
box()
axis(1,at=c(-0.83,0,0.83),labels=c(expression(paste(-0.83,rho)),0,expression(paste(0.83,rho))))
axis(2,at=c(0,0.1,0.1+gaussCor(-0.83),max(y)),
     labels=c(0,expression(nu^2),expression(nu^2+tau^2/2),expression(nu^2+tau^2)))
lines(-0.83*c(1,1),c(-1,0.1+gaussCor(-0.83)),lty=2)
lines(0.83*c(1,1),c(-1,0.1+gaussCor(-0.83)),lty=2)
lines(c(-4,0.83),0.1+gaussCor(-0.83)*c(1,1),lty=2)
#dev.off()

##########################################################
# End....
##########################################################
