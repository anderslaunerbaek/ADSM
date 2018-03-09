##################################################
# Example Binomial
rm(list=ls())
n <- 10; y <- 3
theta <- seq(0,1,by=0.01)
plot(theta,dbinom(y,size=n,prob=theta),type="l")

# Likelihood function
Lik.fun <- function(theta,y,n){
    dbinom(y,size=n,prob=theta)}

# Optimize likelihood function (in general not recommented)
optimize(Lik.fun,c(0,1),y=y,n=n,maximum=TRUE)

# log-Likelihood function
logLik.fun <- function(theta,y,n){
  dbinom(y,size=n,prob=theta,log=TRUE)}

optimize(logLik.fun,c(0,1),y=y,n=n,maximum=TRUE)




##################################################
# Example: Normal
##################################################
y <- c(4.6,6.3,5)

Lik.fun <- function(theta,y){
    prod(dnorm(y,mean=theta,sd=1))
}

optimize(Lik.fun,c(-10,10),y=y,
         maximum=TRUE)



theta <- seq(3.5,7,by=0.01)
plot(theta,sapply(theta,Lik.fun,y=y),type="l")
lines(c(1,1)*mean(y),c(0,0.1),lty=2)
rug(y)

# Numerical derivatives (you will need this package)
library(numDeriv)

# log likelihood
ll.fun <- function(theta,y){
    sum(dnorm(y,mean=theta,sd=1,log=TRUE))
}

# Observed Information
-hessian(ll.fun,5.3,y=y)


##################################################
# Estimate both mean and var
ll.fun2 <- function(theta,y){
     -sum(dnorm(y,mean=theta[1],sd=sqrt(theta[2]),log=TRUE))
}



nlminb(c(0,1),ll.fun2,y=y)

# Invariance
ll.fun2 <- function(theta,y){
     -sum(dnorm(y,mean=theta[1],sd=exp(0.5*theta[2]),log=TRUE))
}


OPt <- nlminb(c(0,1),ll.fun2,y=y)
OPt$par
exp(OPt$par[2])
hessian(ll.fun2,OPt$par,y=y)
# Note independence between sigma and mu


##################################################
## Invariance
(fit1 <- glm(cbind(3,7)~1,family=binomial))
coef(fit1)
1/(1+exp(-coef(fit1)))


##################################################
# Example Poisson regression
##################################################

### Data
x <- 1:10
y <- c(3, 0, 4, 5, 6, 4, 9, 7, 4, 10)
plot(x,y,pch=19)

#### Likelihood poisson
log.lik1 <- function(th) {
    sum(dpois(y, exp(th[1]+th[2] * x ), log = TRUE))
}

#### Find MLE and hessian
fit1 <- optim(par = c(0, 0), fn = log.lik1, hessian = TRUE,
              control=list(fnscale=-1))
fit1

# Wald statistics and and parameter s.e.
se <- sqrt(diag(solve(-fit1$hessian))) 
Z <- fit1$par/se
pv <- 2 * (1 - pnorm(abs(Z)))
round(cbind(Par=fit1$par, se, Z, pv), digits=5)

# Directly in R
fit.glm <- glm(y ~ x, family = poisson)
summary(fit.glm)

# plot it..
plot(x,y,pch=19,cex=1)
curve(exp(fit1$par[1]+fit1$par[2] * x),add=TRUE,col=2,lwd=2)

##################################################
##
log.lik2 <- function(th) {
    sum(dpois(y, exp(th[1]), log = TRUE))
}

#### Find MLE and hessian
fit2 <- optim(par = c(0), fn = log.lik2, hessian = TRUE,
              control=list(fnscale=-1))
fit2

## Deviance
D <- 2*(fit1$value - fit2$value)

## p-value
1-pchisq(D,df=1)

## Directly in R
fit2.glm <- glm(y ~ 1, family = poisson)
anova(fit2.glm,fit.glm,test="Chisq")


################################################################
## Profile likelihood
################################################################

#### Poisson regression example
Plog.lik <- function(beta1) {
    fun.tmp <- function(beta0,beta1){
        -sum(dpois(y, exp(beta0+beta1 * x ), log = TRUE))
    }
    optimize(fun.tmp,c(-10,10),beta1=beta1)$objective   
}

## Plot profile log likelihood
confint(fit.glm,level=0.99)
beta1 <- seq(0.012,0.27,length=100)
plot(beta1,-sapply(beta1,Plog.lik),type="l")

## subtract mle
plot(beta1,-sapply(beta1,Plog.lik)-fit1$value,type="l")
lines(range(beta1),-0.5*qchisq(0.95,df=1)*c(1,1))
confint(fit.glm)
##################################################
# End.......
##################################################
