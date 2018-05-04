##################################################
# Example 6.1
# data 
z <- c(0,1,2,3)
d <- c(803,100,14,3)
##################

# moment est.
mu <- sum(z * d) / sum(d)
s2 <- sum((z - mu) ^ 2 * d) / sum(d)
mu;s2
###################################

###############################################
# glm possion
z.tmp <- c(rep(0,d[1]), rep(1,d[2]), rep(2,d[3]), rep(3,d[4]))
fit1 <- glm(z.tmp ~ 1, family=poisson())
summary(fit1)
anova(fit1)

exp(coef(fit1))
mu
############################################


nll <- function(lambda,z,d){
    -sum((d*dpois(z[1:4],exp(lambda),log=TRUE))) }


nll(0,z,d)
op1 <- nlminb(0,nll,z=z,d=d)
op1

## Taking care of censoring
nll <- function(lambda,z,d){
    -sum((d[1:3]*dpois(z[1:3],exp(lambda),log=TRUE))+d[4]*log(1-ppois(z[3],exp(lambda))))
    }

nll(0,z,d)
op2 <- nlminb(0,nll,z=z,d=d)

c(dpois(0:2,exp(op1$par)),1-ppois(2,exp(op1$par)))*sum(d)
exp.pois <- c(dpois(0:2,exp(op2$par)),
              1-ppois(2,exp(op2$par)))*sum(d)
exp.pois

##################################
## Goodness of fit test.. (terms collected)
Q <- (803-791.85)^2/791.85+(100-118.78)^2/118.78 +
               (17-8.91-0.46)^2/(8.91+0.46)

1-pchisq(Q,df=2)
######################################

##################################################
# glm quasipossion (oversidpersion)
fit2 <- glm(z.tmp ~ 1, family=quasipoisson())
summary(fit2)
##################################################

##################################################
# Compound Poisson Gamma
# negative log-likelihood

nll <- function(pars){
  alpha <- exp(pars[1])
  p <- 1/(1+exp(pars[2]))
  -sum(c(dnbinom(0:2,size=alpha,prob=p,log=TRUE),
         log(1-pnbinom(2,alpha,p)))* c(803,100,14,3)
       )
}


OPt <- nlminb(c(log(0.1),0), nll)
OPt
c(exp(OPt$par[1]),1/(1+exp(OPt$par[2]))) ## Pars for NBinom model

# expected from model:
exp.nbin <- c(dnbinom(0:2,exp(OPt$par[1]),1/(1+exp(OPt$par[2]))),
  1-pnbinom(2,exp(OPt$par[1]),1/(1+exp(OPt$par[2])))) *
    sum(d)
# Note: estimate in book use l(Y) = l(0) + l(1)+ l(2) + l(3)

plot(z,log10(d),pch=19,ylim=c(-1,3),xlim=c(-0.2,3.2),type="p")
points(z-0.1,log10(exp.pois),type="p")
points(z+0.1,log10(exp.nbin),type="p",pch=3)

# unceartainty by Wald statistics
library(numDeriv)
H <- hessian(nll,OPt$par)

V <- solve(H)
cov2cor(V)
exp(OPt$par[1]+2*c(-1,1)*sqrt(diag(V)[1]))
1/(1+exp(OPt$par[2]-2*c(-1,1)*sqrt(diag(V)[2])))
##################################################

# unceartainty by profile likelihood
Pnll <- function(p1,p2){
    pars <- c(p1,p2)
    nll(pars)    
}


OPt$par[2]+3*c(-1,1)*sqrt(diag(V)[2])

beta <- seq(-2.8,-0.3,by=0.05)
pb <- beta*0
for(i in 1:length(beta)){
    pb[i] <- nlminb(OPt$par[1],Pnll,p2=beta[i])$objective
}

plot(1/(1+exp(beta)),exp(-(pb-OPt$objective)),type="l")
lines(range(1/(1+exp(beta))),
      exp(-0.5*qchisq(0.95,df=1))*c(1,1),
      col=2,lwd=2)
1/(1+exp(OPt$par[2]-2*c(-1,1)*sqrt(diag(V)[2])))

##################################################
#
##################################################




##################################################
# Example 6.3
##################################################
lids <- data.frame(
  def = 0:9,
  samples = c(131, 38, 28, 11, 4, 5, 5, 2, 3, 2),
  n = rep(770,10)
)
lids

p <- sum(lids$samples*(0:9)) /
  sum(lids$samples)/770
p

lids2 <- data.frame(y=rep(lids$def,lids$samples),n=770)
fit.glm<- glm(cbind(y,n-y)~1,data=lids2,family=binomial())
summary(fit.glm)
1/(1+exp(-coef(fit.glm)))

################################################################
plot(lids$def,lids$samples,log="y",type="h",ylim=c(0.01,120))
points(lids$def,lids$samples,pch=19)

lines(lids$def-0.2,dbinom(lids$def,prob=p,size=lids$n)*sum(lids$samples),
      type="h")
points(lids$def-0.2,dbinom(lids$def,prob=p,size=lids$n)*sum(lids$samples))

## Likelihood
ldCBb <- function(y,alpha,beta,n){
  lchoose(n,y)+lgamma(alpha + y) + lgamma(beta + n - y) -
    lgamma(alpha + beta + n) - lgamma(alpha) - 
    lgamma(beta)+ lgamma(alpha+beta)}

## Calculate from slide formulation of mu and phi
CBinBeta <- function(pars,y,n){
  mu <- 1/(1+exp(-(pars[1])))
  phi <- exp(pars[2])
  alpha <- mu*phi
  beta <- (1-mu)*phi
  -sum(ldCBb(y,alpha,beta,n))
}



OPt <- nlminb(c(1,0.5), CBinBeta, y=lids2$y, n=lids2$n)
mu <- 1/(1+exp(-(OPt$par[1])))
phi <- exp(OPt$par[2])
alpha <- mu*phi
beta <- (1-mu)*phi
alpha
beta


# plot
lines(lids$def+0.2,exp(ldCBb(lids$def,alpha,beta,lids$n))*sum(lids$samples),
      type="h")
points(lids$def+0.2,exp(ldCBb(lids$def,alpha,beta,lids$n))*sum(lids$samples),
       pch=2)

# Just to compare with book ######################
mpi<- alpha/(alpha+beta)
gamma<- 1/(alpha+beta)
mpi;gamma
##################################################
##################################################

########################################################
# Beetle example
conc <- c(24.8, 24.6, 23, 21, 20.6, 18.2, 16.8, 15.8, 14.7, 10.8)
n <- c(30, 30, 31, 30, 26, 27, 31, 30, 31, 24)
y <- c(23, 30, 29, 22, 23, 7, 12, 17, 10, 0)

plot(log(conc),y/n,pch=19)

resp <- cbind(y, n - y)
fit <- glm(resp ~ I(log(conc)), family = binomial())
 
fit
  
1 - pchisq(fit$deviance, fit$df.residual)
  
pred<-predict(fit,se=TRUE)
points(log(conc),1/(1+exp(-pred$fit)),col=2,pch=19,type="b")
for(i in 1:length(pred$fit)){
  lines(log(conc[i])*c(1,1),1/(1+exp(-pred$fit[i]+c(-2,2)*pred$se.fit[i])),
        col=2,pch=19)
}


CBinBeta <- function(pars,y,n,x){
    N <- length(y)
    mu <- 1/(1+exp(-(pars[1])))
    phi <- exp(pars[2])
    alpha <- mu*phi
    beta <- (1-mu)*phi
    -sum(ldCBb(y,alpha,beta,n))
}

    
OPt <- nlminb(c(0.5,0.5), CBinBeta, y=y, n=n)

fit <- glm(resp ~ 1, family = binomial())
fit

# parameters
1/(1+exp(-coef(fit)))
1/(1+exp(-OPt$par[1]))
sum(y)/sum(n)


##################################################
# regression model
CBinBetaReg <- function(pars,y,n,x){
  N <- length(y)
  mu <- 1/(1+exp(-(pars[1]+pars[2] * x)))
  phi <- exp(pars[3])
  alpha <- mu*phi
  beta <- (1-mu)*phi
  -sum(ldCBb(y,alpha,beta,n))
  }

OPt1 <- nlminb(c(0.5,0.5,0.5), CBinBetaReg, y=y, n=n,x=conc)

mu <- 1/(1+exp(-(OPt1$par[1]+OPt1$par[2] * conc)))
phi <- exp(OPt1$par[3])

alpha <- mu * phi
beta <- (1-mu) * phi

plot(log(conc),y/n,pch=19)
lines(log(conc),mu)

resp <- cbind(y, n - y)
fit <- glm(resp ~ I(log(conc)), family = binomial())
 
fit
  
1 - pchisq(fit$deviance, fit$df.residual)
  
pred<-predict(fit,se=TRUE)
points(log(conc),1/(1+exp(-pred$fit)),col=2,pch=19,type="b")
for(i in 1:length(pred$fit)){
 lines(log(conc[i])*c(1,1), 1/(1+exp(-pred$fit[i]+c(-2,2)*pred$se.fit[i])),
        col=2,pch=19)
}

##################################################
# Conditional pars see Th. 6.6
alpha.z <- alpha + y
beta.z <- beta + n - y
mu.z <- alpha.z/(alpha.z+beta.z)
sigma.z <- alpha.z+beta.z

# Plot results
lines(log(conc), mu.z, pch=19, col=3, type="b")

lines(log(conc),
      qbeta(0.025,alpha,beta),
      col=1)
lines(log(conc),
      qbeta(0.975,alpha,beta),
      col=1)

lines(log(conc),
      qbeta(0.025,alpha.z,beta.z),
      col=3)
lines(log(conc),
      qbeta(0.975,alpha.z,beta.z),
      col=3)

## Compare with TMB
library(glmmTMB)
dat <-  data.frame(y=y,n=n,conc=conc,obs=1:10)
fitTMB<-glmmTMB(cbind(y,n-y)~I(log(conc))+(1|obs),family=binomial,data=dat)
## plot random effects
qqnorm(fitTMB$fit$parfull[3:12])
qqline(fitTMB$fit$parfull[3:12]) 

## Compare estimated mean values 
eta.tmb<-fitTMB$fit$parfull[1]+fitTMB$fit$parfull[2]*log(conc)+fitTMB$fit$parfull[3:12]
mu.tmb <- 1/(1+exp(-eta.tmb))
round(cbind(mu.z,mu.tmb),digits=2)

#####################################################################
#
#####################################################################
