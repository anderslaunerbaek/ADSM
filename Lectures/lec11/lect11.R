setwd("~/DTU/Courses/ADSM/Lectures/lec11/")
# rm(list=ls())


#################################################################
# Beetle example
conc <- c(24.8, 24.6, 23, 21, 20.6, 18.2, 16.8, 15.8, 14.7, 10.8)
n <- c(30, 30, 31, 30, 26, 27, 31, 30, 31, 24)
y <- c(23, 30, 29, 22, 23, 7, 12, 17, 10, 0)
dat <- data.frame(y=y,n=n,logc=log(conc))

plot(y/n~logc,data=dat,pch=19)

fit0 <- glm(cbind(y, n - y) ~ logc, family = binomial(),data=dat)
fit0  
1 - pchisq(fit0$deviance, fit0$df.residual)
  
pred<-predict(fit0,se=TRUE)
points(log(conc),1/(1+exp(-pred$fit)),col=2,pch=19,type="b")
for(i in 1:length(pred$fit)){
  lines(log(conc[i])*c(1,1),1/(1+exp(-pred$fit[i]+c(-2,2)*pred$se.fit[i])),
        col=2,pch=19)
}

## Implementation of a miced effect model
X <- matrix(0,ncol=2,nrow=dim(dat)[1])
X[ ,1]<-1
X[ ,2]<-dat$logc

library(numDeriv)
nll <- function(u,beta,sigma.u,X){
    eta <- X%*%beta + u
    mu <- 1/(1+exp(-eta))
    -sum(dbinom(dat$y,prob=mu,size=dat$n,log=TRUE)+
           dnorm(u,sd=sigma.u,log=TRUE))
}


nll.LA <- function(theta,X){
    beta <- theta[1:dim(X)[2]]
    sigma.u <- exp(theta[dim(X)[2]+1])
    est <- nlminb(rep(0,length(dat$y)),objective = nll, beta=beta, sigma.u=sigma.u,X=X)
    u <- est$par
    l.u <- est$objective
    H <- hessian(func = nll, x = u, beta = beta, sigma.u = sigma.u, X=X)
    l.u + 0.5 * log(det(H/(2*pi)))
}


system.time((fit1 <- nlminb(c(0,0,1),nll.LA,X=X)))
-fit1$objective
logLik(fit0)
1-pchisq(2*(-fit1$objective-logLik(fit0)),df=1)
## Hence a significant improvement


##################################################
## improve by using independence of u's in nlminb
nll.LA3 <- function(theta,X){
    beta <- theta[1:dim(X)[2]]
    sigma.u <- exp(theta[dim(X)[2]+1])
    fun.tmp <- function(ui,u,beta,sigma.u,X,i){
        u <- u*0
        u[i]<-ui
        nll(u,beta,sigma.u,X)
    }
    u <- numeric(length(dat$y))
    for(i in 1:length(u)){
        u[i] <- nlminb(0,objective = fun.tmp, u=u,beta=beta, sigma.u=sigma.u,X=X,i=i)$par
    }
    l.u <- nll(u,beta,sigma.u,X)
    H <- numeric(length(u))
    for(i in 1:length(u)){
        H[i] <- hessian(func = fun.tmp, x = u[i],u=u, beta = beta, 
                        sigma.u = sigma.u, X=X,i=i)}
    l.u + 0.5 * log(prod(H/(2*pi)))
}

system.time(fit2 <- nlminb(c(0,0,0),nll.LA3,X=X))

fit1$objective
fit2$objective

##################################################
## Improve by implementing the derivative

## Implement the theoretical first and second derivative
## wrt u
Dnll <- function(u,beta,sigma.u,X){
    y <- dat$y
    n <- dat$n
    eta <- X%*%beta + u
    mu <- 1/(1+exp(-eta))
    Dmu <- mu^2*exp(-eta)
    (n-y)/(1-mu)*Dmu-y/mu*Dmu+u/sigma.u^2
}

Hnll <- function(u,beta,sigma.u,X){
    y <- dat$y
    n <- dat$n
    eta <- X%*%beta + u
    mu <- 1/(1+exp(-eta))
    Dmu <- mu^2*exp(-eta)
    D2mu <- mu^2*exp(-eta)*(2*mu*exp(-eta)-1)
    diag(as.numeric(D2mu * ((n-y)/(1-mu)-y/mu) +
                    Dmu^2 * ((n-y)/(1-mu)^2 + y/mu^2) + 1/sigma.u^2))
}


nll.LA4 <- function(theta,X){
    beta <- theta[1:dim(X)[2]]
    sigma.u <- exp(theta[dim(X)[2]+1])
    est <- nlminb(rep(0,length(dat$y)),objective = nll,
                  gradient = Dnll, hessian = Hnll,
                  beta=beta, sigma.u=sigma.u,X=X)
    u <- est$par
    l.u <- est$objective
    H <- diag(Hnll(u,beta,sigma.u,X))
    l.u + 0.5 * log(prod(H/(2*pi)))
}
    
system.time(fit3 <- nlminb(c(0,0,0),nll.LA4,X=X))
## Very large improvement in timing
fit3$objective

##################################################
## Check accuracy by importance sampling
nll2 <- function(u,beta,sigma.u,X,uhat){
    eta <- X%*%beta + u
    mu <- 1/(1+exp(-eta))
    -sum(dbinom(dat$y,prob=mu,size=dat$n,log=TRUE)+dnorm(u,mean=uhat,sd=sigma.u,log=TRUE))
}


nll.LA4Rexsamp <- function(theta,X,k,seed){
    set.seed(seed)
    beta <- theta[1:dim(X)[2]]
    sigma.u <- exp(theta[dim(X)[2]+1])
    est <- nlminb(rep(0,length(dat$y)),objective = nll,
                  gradient = Dnll, hessian = Hnll,
                  beta=beta, sigma.u=sigma.u,X=X)
    u <- est$par
    l.u <- est$objective
    H <- diag(Hnll(u,beta,sigma.u,X))
    L <- numeric(k)
    s <- sqrt(1/H)
    for(i in 1:k){
        u.sim <- rnorm(length(u),mean=u,sd=s)
        L[i] <- exp(-nll(u.sim,beta,sigma.u,X))/prod(dnorm(u.sim,mean=u,sd=s))
    }
    -log(mean(L))
}


k <- 100000
L <- nll.LA4Rexsamp(fit3$par,X,k=k,seed=16345)
c(L,fit3$objective)

n<-10000
L.sim <- numeric(n)
for(i in 1:n){
    L.sim[i] <- exp(-nll.LA4Rexsamp(fit3$par,X,k=1,seed=i))
}

cs <- cumsum(L.sim)/(1:n)
css <- cumsum(L.sim^2)/1:n
csse <- css/(1:n)-(cs/1:n)^2


plot(log10(1:n),-log(cs),type="l")#,ylim=range(-log(L)))
lines(log10(3:n),-log((cs+2*sqrt(csse))[-(1:2)]),type="l",col=2)#,ylim=range(-log(L)))
lines(log10(3:n),-log((cs-2*sqrt(csse))[-(1:2)]),type="l",col=2)#,ylim=range(-log(L)))
lines(log10(c(1,n)),fit3$objective*c(1,1),col=3)

# using resampling for estimation
system.time((fitRe1 <- nlminb(c(0,0,0),nll.LA4Rexsamp,X=X,
                              k=100,seed=1)))
fitRe1
system.time(fitRe2<-nlminb(c(0,0,0),nll.LA4Rexsamp,X=X,
                           k=1000,seed=1))
fitRe2

## use numerical integration
nllInt <- function(u,eta,sigma.u,n,y){
    eta <- eta + u
    mu <- 1/(1+exp(-eta))
    -(dbinom(y,prob=mu,size=n,log=TRUE)+
        dnorm(u,sd=sigma.u,log=TRUE))
}

ll <- function(u,eta,sigma.u,n,y){
    exp(-nllInt(u,eta,sigma.u,n,y))
}

nll.MInt <- function(theta,X){
    beta <- theta[1:dim(X)[2]]
    sigma.u <- exp(theta[dim(X)[2]+1])
    L <- numeric(dim(X)[1])
    eta <- as.numeric(X%*%beta)
    for(i in 1:length(L)){
        L[i] <- integrate(ll,-Inf,Inf,eta=eta[i],sigma.u=sigma.u,
                          n=dat$n[i],y=dat$y[i])$value
    }
    -sum(log(L))
}


L.True <- nll.MInt(fit3$par,X)
c(fit3$objective,L)-L.True

lines(log(c(1,k)),L.True*c(1,1),col=4)

##################################################
## glmmTMB
library(glmmTMB)

dat <- cbind(obs = 1:10,dat)

(fitTMB <- glmmTMB(cbind(y,n-y) ~ logc +(1|obs),
                   family=list(family="binomial",link="logit"),
                   data=dat))
fit3$objective

## or the default
system.time(fitTMB <- glmmTMB(cbind(y,n-y) ~ logc +(1|obs),
                              family="binomial",data=dat))


##################################################
## Implement in TMB
library(TMB)                          

compile("beetle.cpp")
dyn.load(dynlib("beetle"))
parameters <- list(B = rep(0, 10),
                   beta = c(0,0),
                   sigma_b = 1
                   )


obj <- MakeADFun(data = dat,
                 parameters = parameters,
                 random = c("B"),
                 DLL = "beetle"
                 )



system.time(opt <- nlminb(obj$par, obj$fn, obj$gr))
opt$par
opt$objective




##################################################
## Orange tree examples 1
##################################################

library(datasets)
data(Orange)
Orange

tree <- factor(rep(1:5,each=7))
data <- list(tree = tree, y = Orange$circumference, t=Orange$age)

compile("orange.cpp")
dyn.load(dynlib("orange"))


parameters <- list(u = rep(1, nlevels(Orange$Tree)),
                   beta = c(mean(data$y),mean(data$t),diff(range(data$t))),
                   sigma_u= 1,
                   sigma = 1
                   )

obj <- MakeADFun(data = data,
                 parameters = parameters,
                 random = "u",
                 DLL = "orange"
                 )


## Test eval function and gradient
obj$fn()
obj$gr()

## Fit model
opt <- nlminb(obj$par, obj$fn, obj$gr,lower=0.01)
opt$par
opt$objective

## report on result
sdreport(obj)


rap <- sdreport(obj,getJointPrecision = TRUE)
summary(rap,"random")
rap$par.random
rap$diag.cov.random
names(rap)

## Plot the results
mean.fun <- function(beta,u,t){
    (beta[1] + u)/(1 + exp(-(t-beta[2])/beta[3]))
}
names(rap)
rap$par.fixed
mean.fun(rap$par.fixed[1:3],rap$par.random[1],data$t[1:7])

plot(data$t[1:7],mean.fun(rap$par.fixed[1:3],rap$par.random[1],data$t[1:7]),type="l",ylim=c(0,220))
lines(data$t[1:7],mean.fun(rap$par.fixed[1:3],rap$par.random[2],data$t[1:7]),lty=2)
lines(data$t[1:7],mean.fun(rap$par.fixed[1:3],rap$par.random[3],data$t[1:7]),lty=3)
lines(data$t[1:7],mean.fun(rap$par.fixed[1:3],rap$par.random[4],data$t[1:7]),lty=4)
lines(data$t[1:7],mean.fun(rap$par.fixed[1:3],rap$par.random[5],data$t[1:7]),lty=5)
points(data$t,data$y,pch=19,col=data$tree)


##################################################
## Orange 2
##################################################
tree <- factor(rep(1:5,each=7))
times <- factor(Orange$age)
data <- list(tree = tree, y = Orange$circumference, 
             t=Orange$age,times=times)

compile("orange2.cpp")
dyn.load(dynlib("orange2"))


parameters <- list(u1 = rep(1, nlevels(Orange$Tree)),
                   u2 = rep(1, nlevels(data$times)),
                   beta = c(mean(data$y),mean(data$t),
                            diff(range(data$t))),
                   sigma_u1= 1,
                   sigma_u2= 1,
                   sigma = 1
                   )

obj <- MakeADFun(data = data,
                 parameters = parameters,
                 random = c("u1","u2"),
                 DLL = "orange2"
                 )


## Fit model
opt <- nlminb(obj$par, obj$fn, obj$gr,lower=0.01)
opt$par
opt$objective

## Random effects
rap <- sdreport(obj,getJointPrecision = TRUE)
rap$par.random
rap$diag.cov.random

## Covarive of random effects
round(rap$jointPrecision[6:12,1:5],digits=3)
round(rap$jointPrecision[1:5,1:5],digits=3)

## Plot model
mean.fun <- function(beta,u1,u2,t){
    (beta[1]+u1+u2)/(1+exp(-(t-beta[2])/beta[3]))
}
names(rap)
rap$par.fixed
mean.fun(rap$par.fixed[1:3],rap$par.random[1],
         rap$par.random[6:12],data$t[1:7])

plot(data$t[1:7],mean.fun(rap$par.fixed[1:3],rap$par.random[1],rap$par.random[6:12],data$t[1:7]),type="l",ylim=c(0,220))
lines(data$t[1:7],mean.fun(rap$par.fixed[1:3],rap$par.random[2],rap$par.random[6:12],data$t[1:7]),lty=2)
lines(data$t[1:7],mean.fun(rap$par.fixed[1:3],rap$par.random[3],rap$par.random[6:12],data$t[1:7]),lty=3)
lines(data$t[1:7],mean.fun(rap$par.fixed[1:3],rap$par.random[4],rap$par.random[6:12],data$t[1:7]),lty=4)
lines(data$t[1:7],mean.fun(rap$par.fixed[1:3],rap$par.random[5],rap$par.random[6:12],data$t[1:7]),lty=5)
points(data$t,data$y,pch=19,col=data$tree)


##################################################
## Orange 3
##################################################
tree <- factor(rep(1:5,each=7))
times <- factor(Orange$age)
season <- rep(c(-1,-1,1,1,-1,1,-1),5)/2
data <- list(tree = tree, y = Orange$circumference, t=Orange$age,times=times,season=season)

compile("orange3.cpp")
dyn.load(dynlib("orange3"))


parameters <- list(u = rep(1, nlevels(Orange$Tree)),
                   beta = c(mean(data$y),mean(data$t),
                            diff(range(data$t)),0),
                   sigma_u= 1,
                   sigma = 1
                   )

obj <- MakeADFun(data = data,
                 parameters = parameters,
                 random = c("u"),
                 DLL = "orange3"
                 )


## Fit model
opt <- nlminb(obj$par, obj$fn, obj$gr,lower=c(-Inf,-Inf,-Inf,-Inf,0.01,0.01))
opt$par
opt$objective


rap <- sdreport(obj,getJointPrecision = TRUE)
rap$par.random
rap$diag.cov.random
names(rap)

## Plot model
mean.fun <- function(beta,u,t){
    (u)/(1+exp(-(t-beta[2])/beta[3]))
}
names(rap)
rap$par.fixed
mean.fun(rap$par.fixed[1:3],rap$par.random[1],data$t[1:7])

plot(data$t[1:7],mean.fun(rap$par.fixed[1:3],rap$par.random[1],data$t[1:7]),type="l",ylim=c(0,220))
lines(data$t[1:7],mean.fun(rap$par.fixed[1:3],rap$par.random[2],data$t[1:7]),lty=2)
lines(data$t[1:7],mean.fun(rap$par.fixed[1:3],rap$par.random[3],data$t[1:7]),lty=3)
lines(data$t[1:7],mean.fun(rap$par.fixed[1:3],rap$par.random[4],data$t[1:7]),lty=4)
lines(data$t[1:7],mean.fun(rap$par.fixed[1:3],rap$par.random[5],data$t[1:7]),lty=5)
points(data$t,data$y,pch=19,col=data$tree)


##################################################
## Orange 4
##################################################

tree <- factor(rep(1:5,each=7))
times <- factor(Orange$age)
data <- list(tree = tree, y = Orange$circumference, 
             t=Orange$age,times=times,n=35)

compile("orange4.cpp",libtmb = FALSE,libinit=FALSE,openmp=FALSE)
dyn.load(dynlib("orange4"))
parameters <- list(u1 = rep(1, nlevels(Orange$Tree)),
                   u2 = rep(1, nlevels(data$times)),
                   beta = c(mean(data$y),mean(data$t),
                            diff(range(data$t))),
                   sigma_u1= 1,
                   sigma_u2= 1,
                   sigma = 1,
                   phi = 1
                   )


parameters <- list(u1 = rep(1, nlevels(Orange$Tree)),
                   u2 = rep(1, nlevels(data$times)),
                   beta = c(192.1,857.5,348.1),
                   sigma_u1= 32.7,
                   sigma_u2= 12,
                   sigma = 6.12,
                   phi = 0
                   )

obj <- MakeADFun(data = data,
                 parameters = parameters,
                 random = c("u1","u2"),
                 DLL = "orange4"
                 )



opt <- nlminb(obj$par, obj$fn, obj$gr,lower=0.01)
opt$par
opt$objective



rap <- sdreport(obj,getJointPrecision = TRUE)
rap$par.random
rap$diag.cov.random
names(rap)


rap$jointPrecision[6:12,1:5]




## Plot model
mean.fun <- function(beta,u1,u2,t){
    (beta[1]+u1+u2)/(1+exp(-(t-beta[2])/beta[3]))
}
names(rap)
rap$par.fixed
mean.fun(rap$par.fixed[1:3],rap$par.random[1],
         rap$par.random[6:12],data$t[1:7])

plot(data$t[1:7],mean.fun(rap$par.fixed[1:3],rap$par.random[1],rap$par.random[6:12],data$t[1:7]),type="l",ylim=c(0,220))
lines(data$t[1:7],mean.fun(rap$par.fixed[1:3],rap$par.random[2],rap$par.random[6:12],data$t[1:7]),lty=2)
lines(data$t[1:7],mean.fun(rap$par.fixed[1:3],rap$par.random[3],rap$par.random[6:12],data$t[1:7]),lty=3)
lines(data$t[1:7],mean.fun(rap$par.fixed[1:3],rap$par.random[4],rap$par.random[6:12],data$t[1:7]),lty=4)
lines(data$t[1:7],mean.fun(rap$par.fixed[1:3],rap$par.random[5],rap$par.random[6:12],data$t[1:7]),lty=5)
points(data$t,data$y,pch=19,col=data$tree)
##################################################
##
##################################################




##################################################
## Orange 5
##################################################

tree <- factor(rep(1:5,each=7))
times <- factor(Orange$age)
season <- rep(c(-1,-1,1,1,-1,1,-1),5)/2
data <- list(tree = tree, y = Orange$circumference, 
             t=Orange$age,times=times,season=season,ntimes=7)

compile("orange5.cpp",libtmb = FALSE,libinit=FALSE,openmp=FALSE)
dyn.load(dynlib("orange5"))
parameters <- list(u = rep(1, nlevels(Orange$Tree)),
                   beta = c(mean(data$y),mean(data$t),
                            diff(range(data$t)),0),
                   sigma_u= 1,
                   sigma = 1,
                   phi = 0)


parameters <- list(u = rep(1, nlevels(Orange$Tree)),
                   beta = c(216,859,437,0.3),
                   sigma_u= 36.7,
                   sigma = 5.8,
                   phi = 1
                   )



obj <- MakeADFun(data = data,
                 parameters = parameters,
                 random = c("u"),
                 DLL = "orange5"
                 )


## Fit model
opt <- nlminb(obj$par, obj$fn, obj$gr)
opt$par
opt$objective


rap <- sdreport(obj,getJointPrecision = TRUE)
rap$par.random
rap$diag.cov.random
names(rap)


rap$jointPrecision[1:5,1:5]




## Plot model
mean.fun <- function(beta,u,t,s){
    (u)/(1+exp(-(t-beta[2])/beta[3]+beta[4]*s))
}
names(rap)
rap$par.fixed
mean.fun(rap$par.fixed[1:4],rap$par.random[1],
         data$t[1:7],data$season[1:7])

plot(data$t[1:7],mean.fun(rap$par.fixed[1:4],
                          rap$par.fixed[1],data$t[1:7],
                          data$season[1:7]),type="l",ylim=c(0,220))
for(i in 1:5){
  lines(data$t[1:7],mean.fun(rap$par.fixed[1:4],rap$par.random[i],
                             data$t[1:7],
                             data$season[1:7]),lty=i+1)
}
points(data$t,data$y,pch=19,col=data$tree)


##################################################
##
##################################################
