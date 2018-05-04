## Solution to exercise....
##################################################
## Example 5.7 (textbook)
##################################################
setwd("~/Kurser/kursus02424/2018/lectures/lect11/")
dat <- read.table("seed.dat", sep = ";",  header = TRUE)

summary(fit1<-glm(cbind(y,n-y)~factor(variety)+factor(root),family=binomial,data=dat))
1-pchisq(39.7,df=18)

X <- model.matrix(fit1)

nll <- function(u,beta,sigma.u,X){
    eta <- X%*%beta + u
    mu <- 1/(1+exp(-eta))
    -sum(dbinom(dat$y,prob=mu,size=dat$n,log=TRUE)+dnorm(u,sd=sigma.u,log=TRUE))
}

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
    
system.time(fit5 <- nlminb(c(0,0,0,0),nll.LA4,X=X))
fit5

library(numDeriv)
H <- hessian(nll.LA4,x=fit5$par,X=X)
V <- solve(H)
se <- sqrt(diag(V))

summary(fit1)
## unceartainty of parameters
z <- fit5$par/se
round(cbind(fit5$par,se,z , 2*(1-pnorm(abs(z)))),digits=3)

##################################################
## Importance sampling
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
L <- nll.LA4Rexsamp(fit5$par,X,k=k,seed=16345)
c(L,fit5$objective)

n<-10000
L.sim <- numeric(n)
for(i in 1:n){
    L.sim[i] <- exp(-nll.LA4Rexsamp(fit5$par,X,k=1,seed=i))
}

cs <- cumsum(L.sim)/(1:n)
css <- cumsum(L.sim^2)/1:n
csse <- css/(1:n)-(cs/1:n)^2

lower <- -log((cs+2*sqrt(csse))[-(1:2)])
upper <- -log((cs-2*sqrt(csse))[-(1:2)])
plot(log10(1:n),-log(cs),type="l",ylim=range(c(upper[-c(1:1000)],lower[-c(1:1000)])))
lines(log10(3:n),lower,type="l",col=2)#,ylim=range(-log(L)))
lines(log10(3:n),upper,type="l",col=2)#,ylim=range(-log(L)))
lines(log10(c(1,n)),fit5$objective*c(1,1),col=3)


## use numerical integration
nllInt <- function(u,eta,sigma.u,n,y){
    eta <- eta + u
    mu <- 1/(1+exp(-eta))
    -(dbinom(y,prob=mu,size=n,log=TRUE)+dnorm(u,sd=sigma.u,log=TRUE))
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
        L[i] <- integrate(ll,-Inf,Inf,eta=eta[i],sigma.u=sigma.u,n=dat$n[i],y=dat$y[i])$value
    }
    -sum(log(L))
}


L.True <- nll.MInt(fit5$par,X)
c(fit5$objective,L)-L.True

lines(log(c(1,k)),L.True*c(1,1),col=4)
##################################################
## glmmTMB
dat <- cbind(plate=1:21,dat)
library(glmmTMB)
glmmTMB(cbind(y,n-y)~factor(variety) + factor(root)+(1|plate),family="binomial",data=dat)
fit5$objective
## same model 


##################################################
colnames(dat) <- c("plate","extract","seed","r","n")

library(TMB)                          

compile("seed.cpp")
dyn.load(dynlib("seed"))
parameters <- list(B = rep(0, 21),
                   beta = c(0,0,0),
                   sigma_b = 1
                   )


obj <- MakeADFun(data = dat,
                 parameters = parameters,
                 random = c("B"),
                 DLL = "seed"
                 )



opt <- nlminb(obj$par, obj$fn, obj$gr)
opt$par
opt$objective

##################################################
## Done!
##################################################
