setwd("~/Kurser/kursus02424/2018/lectures/lect10")
##################################################
#  Example 5.6, not the same result though (misprint in the book)
flow <- read.table("flow.txt",header=TRUE)
pch <- c(3 ,4 ,15 ,19,17,1)
plot(flow[flow$flow==0.1,"dev"],flow[flow$flow==0.5,"dev"],
     ylim=c(-7,12),xlim=c(-7,12),
     pch=pch[flow$meter[flow$flow==0.1]-40])

## Moment estimates eqs. (5.85)-(5.92):
## Average (eq. 5.85)
mu <- c(mean(flow$dev[flow$flow==0.1]),
        mean(flow$dev[flow$flow==0.5]))

## SSE (eq (5.87)), mu_i (eq (5.86))
SSE <- matrix(0,ncol=2,nrow=2)
mui <- matrix(0,ncol=6,nrow=2)
meter <- 41:46
flows <- c(0.1,0.5)
for(i in 1:6){
    mui[,i] <- c(mean(flow$dev[flow$flow==0.1 &
                               flow$meter==meter[i]]),
                 mean(flow$dev[flow$flow==0.5 &
                               flow$meter==meter[i]]))
    SSE <- SSE + 2*var(cbind(flow$dev[flow$meter==meter[i] &
                                      flow$flow==0.1],
                             flow$dev[flow$meter==meter[i] &
                                      flow$flow==0.5]))
}
mui
SSE

## SST (eq. (5.89))
SST<-17 * var(cbind(flow$dev[flow$flow==0.1],
                    flow$dev[flow$flow==0.5]))
SST

## Between group variation (eq (5.88))
SSB <- 5 * 3 * var(t(mui))
SSB

## Moment estimates (Theo. 5.11)
Sigma.m <- SSE/12
Sigma0.m <- (SSB/5 - SSE/12)/3

## Are these cavariance matrices?
cov2cor(Sigma.m)
cov2cor(Sigma0.m)

## Within meter covariance matrix
SigBlock <- kronecker(matrix(1,ncol=3,nrow=3),Sigma0.m) +
    kronecker(diag(3),Sigma.m)
Sig <- kronecker(diag(6),SigBlock)

SigBlock

## Within meter corelation matrix
round(cov2cor(SigBlock),digits=TRUE)
#
##################################################

####################################
# Likelihood (from Theorem 5.12)
ll <- function(mu,Sigma,Sigma0,SSE,N,y,ni){
    ll <- 0
    k <- dim(y)[2]    
    for( i in 1:dim(y)[2]){
        ll <- ll -  0.5 * log(det(Sigma/ni[i]+Sigma0)) -
            0.5 * t(mui[ ,i] - mu) %*% solve(Sigma/ni[i]+Sigma0) %*%
                (mui[ ,i] -mu)
    }
    ll <- ll - (N - k) / 2 * log(det(Sigma)) - 0.5 * sum(diag(SSE %*% solve(Sigma)))
    #
    return(-ll)
}

# E.g set up parameters
X <- cbind(1,flow$flow==0.5)
sig1 <- 1; sig2 <- 2; rho <- 0.25
sig01 <- 3; sig02 <- 4; rho0 <- 0.5

Sigma <- Sigma0 <- matrix(ncol=2,nrow=2)
Sigma[1,1] <- sig1; Sigma[2,2] <- sig2
Sigma[1,2] <- Sigma[2,1] <- rho * sqrt(sig1*sig2)
Sigma0[1,1] <- sig01; Sigma0[2,2] <- sig02
Sigma0[1,2] <- Sigma0[2,1] <- rho0 * sqrt(sig01*sig02)
ni <- rep(3,6)
y <- mui
y
N <- sum(ni)

ll(mu,Sigma,Sigma0,SSE,N,y,ni)
##################################################

# Parametrize likelihood (to be function of parameter vector
nll <- function(pars){
    rho0 <- pars["rho0"]
    Sigma <- Sigma0 <- matrix(ncol=2,nrow=2)
    Sigma[1,1] <- pars["sig1"]
    Sigma[2,2] <- pars["sig2"]
    Sigma[1,2] <- Sigma[2,1] <- pars["rho"] * sqrt(pars["sig1"] * pars["sig2"])
    Sigma0[1,1] <- pars["sig01"]
    Sigma0[2,2] <- pars["sig02"]
    Sigma0[1,2] <- Sigma0[2,1] <- rho0 * sqrt(pars["sig01"] * pars["sig02"])
    ll(mu,Sigma,Sigma0,SSE,N,y,ni)
}

pars <- c(sig1=1,sig2=1,rho=0,sig01=1,sig02=1,rho0=1)
nll(pars)
OP1<-nlminb(pars,nll,
            lower=c(0.01,0.01,-0.99,0.01,0.01,-Inf),
            upper = c(Inf,Inf,0.99, Inf, Inf, Inf))
OP1

## Compare vith moment estimates
OP1$par
Sigma.m; cov2cor(Sigma.m)[1,2]
Sigma0.m; cov2cor(Sigma0.m)[1,2]

Sigma.ml <- diag(2); Sigma.ml[1,2] <- Sigma.ml[2,1] <- OP1$par["rho"]
Sigma.ml <- sqrt(diag(OP1$par[1:2])) %*% Sigma.ml %*% sqrt(diag(OP1$par[1:2]))
Sigma.ml; cov2cor(Sigma.ml)

Sigma0.ml <- diag(2); Sigma0.ml[1,2] <- Sigma0.ml[2,1] <- OP1$par["rho0"]
Sigma0.ml <- sqrt(diag(OP1$par[4:5])) %*% Sigma0.ml %*% sqrt(diag(OP1$par[4:5]))
Sigma0.ml; cov2cor(Sigma0.ml)


##################################################
# Likelihood parametrized from correlation structure (1D obs.)
##################################################
nll<- function(pars){
    ## Setting up cavariance matrix
    rho0 <- pars["rho0"]
    Sigma <- Sigma0 <- matrix(ncol=2,nrow=2)
    Sigma[1,1] <- pars["sig1"]
    Sigma[2,2] <- pars["sig2"]
    Sigma[1,2] <- Sigma[2,1] <- pars["rho"] *
        sqrt(pars["sig1"] * pars["sig2"])
    Sigma0[1,1] <- pars["sig01"]
    Sigma0[2,2] <- pars["sig02"]
    Sigma0[1,2] <- Sigma0[2,1] <- rho0 *
        sqrt(pars["sig01"] * pars["sig02"])

    Sig <- matrix(0,ncol=36,nrow=36)
    SigBlock <- matrix(0,ncol=6,nrow=6)
    
    SigBlock <- kronecker(matrix(1,ncol=3,nrow=3),Sigma0) +
        kronecker(diag(3),Sigma)
    Sig <- kronecker(diag(6),SigBlock)
    ##################################################

    ## Measurements and expected value
    y <- flow$dev
    mu <- flow$dev
    mu[flow$flow==0.1] <- mean(flow$dev[flow$flow==0.1])
    mu[flow$flow==0.5] <- mean(flow$dev[flow$flow==0.5])
    ##################################################

    ## Likelihood
    ll <- - 0.5 * log(det(Sig)) - 0.5 * t(y-mu) %*% solve(Sig) %*%
        (y-mu)
    ###########################################################

    return(-ll)
}


OP2 <- nlminb(pars,nll,lower=c(0.01,0.01,-0.99,0.01,0.01,-Inf),
       upper = c(Inf,Inf,0.99, Inf, Inf, Inf))

OP1
OP2
OP1$par-OP2$par # These are the same estimates...

##################################################
# Restricted MLE
nllRe<- function(pars){
    ## Setting up cavariance matrix
    rho0 <- pars["rho0"]
    Sigma <- Sigma0 <- matrix(ncol=2,nrow=2)
    Sigma[1,1] <- pars["sig1"]
    Sigma[2,2] <- pars["sig2"]
    Sigma[1,2] <- Sigma[2,1] <- pars["rho"] *
        sqrt(pars["sig1"] * pars["sig2"])
    Sigma0[1,1] <- pars["sig01"]
    Sigma0[2,2] <- pars["sig02"]
    Sigma0[1,2] <- Sigma0[2,1] <- rho0 *
        sqrt(pars["sig01"] * pars["sig02"])

    SigBlock <- kronecker(matrix(1,ncol=3,nrow=3),Sigma0) +
        kronecker(diag(3),Sigma)
    Sig <- kronecker(diag(6),SigBlock)
    ##################################################

    ## Measurements and expected value
    y <- flow$dev
    mu <- flow$dev
    mu[flow$flow==0.1] <- mean(flow$dev[flow$flow==0.1])
    mu[flow$flow==0.5] <- mean(flow$dev[flow$flow==0.5])
    ##################################################

    ## Design matrix
    X <- matrix(0,ncol=2,nrow=36)
    X[flow$flow==0.1 ,1] <-1
    X[flow$flow==0.5 ,2] <-1
    ##################################################


    ## Restricted likelihoof function
    ll <- - 0.5 * log(det(Sig)) - 0.5 * t(y-mu) %*% solve(Sig) %*%
        (y-mu) - 0.5 * log(det(t(X)%*%solve(Sig)%*%X))
    ##################################################

    return(-ll)
}

pars <- c(sig1=1,sig2=1,rho=0,sig01=1,sig02=1,rho0=0)


OP3 <- nlminb(pars,nllRe,lower=c(0.01,0.01,-0.99,0.01,0.01,-Inf),
       upper = c(Inf,Inf,0.99, Inf, Inf, Inf))

OP3$par
OP2$par

OP3$par
Sigma.m; cov2cor(Sigma.m)[1,2]
Sigma0.m; cov2cor(Sigma0.m)[1,2]
## Restricted maximum likelihood and moment estimates are exactly
## the same
##################################################

##################################################################
# Example: Orange tree
library(datasets)
data(Orange)
head(Orange)
summary(Orange)
plot(circumference~age,data=Orange,col=Tree,pch=19)


## Simple model, no random effects (but nonlinear regression)
Mod0 <- nls(circumference~SSlogis(age, Asym, xmid, scal),data=Orange)
summary(Mod0)

lines(1:2000,SSlogis(1:2000, coef(Mod0)["Asym"], coef(Mod0)["xmid"],
                     coef(Mod0)["scal"]))

plot(Orange$age,resid(Mod0),pch=19,col=Orange$Tree)
# Something is wrong....


######################################################
# First mixed effect model
Orange$Tree<-sort(as.numeric(Orange$Tree)) # To avoid "funny"
                    ## sorting effects

## Mean of observation (include random effect).
mu <- function(beta, u, tree, age) {
  (beta[1] + u[tree])/(1 + exp(-(age - beta[2])/beta[3]))
  }

## log-likelihood from random effect
l.u <- function(u, s.u) {
  sum(dnorm(u, mean = 0, sd = s.u, log = TRUE))
}

## Conditional log-likelihood of obs. given random effect
l.cir <- function(cir, u, b, s) {
  mv <- mu(b, u, Orange$Tree, Orange$age)
  sum(dnorm(cir, mean = mv, sd = s, log = TRUE))
}

## Negative joint log likelihood of theta 
nl <- function(th, u, cir) {
  -l.cir(cir, u, th[1:3], exp(th[4])) - l.u(u, exp(th[5]))
  }

 library(numDeriv)
 l.LA <- function(th) {
   # Note implementation differ slightly from book
   u.init <- rep(0, 5)
       obj <- function(u) nl(th, u, Orange$cir)
       est <- nlminb(u.init, obj)
       lval <- est$obj
       u <- est$par
       H <- hessian(obj, u)
       nll <- lval + 0.5 * log(det(H)) - length(u)/2 * log(2 * pi)
       return(list(nll=nll,u=u))
   }

obj.LA <- function(th){l.LA(th)$nll} 
 
fit <- nlminb(c(200, 700, 200, 0, 0), obj.LA)
H <- hessian(obj.LA, fit$par)

fit
round(cbind(est = fit$par, sd = sqrt(diag(solve(H)))),digits=2)

exp(fit$par[4:5]) # To compare with Table 5.6

u <- l.LA(fit$par)$u

trees<-1:5
plot(Orange$age,Orange$circumference,col=rep(1:5,each=7),pch=19,ylim=c(0,250),type="p")
for(i in 1:5){
  lines(1:2000,mu(fit$par,u,i,1:2000),col=i)
}

res <- function(tree){
  Orange$circumference[Orange$Tree==tree] -
    mu(fit$par,u,tree,Orange$age[Orange$Tree==tree])}

plot(Orange$age[Orange$Tree==1],res(1),col=1,pch=19,ylim=c(-15,15),type="b")
for(i in 2:5){
  points(Orange$age[Orange$Tree==i],res(i),col=i,pch=19,type="b")
} # Still something to be done...
## See book example for modellng of this.. 

#####################################################
# Profile likelihood
Pl.LA <- function(th,sig.u){
  obj.LA(c(th,sig.u))
}

n <-7
sig.u <- seq(fit$par[5]-3*sqrt(diag(solve(H)))[5],
             fit$par[5]+3*sqrt(diag(solve(H)))[5],
             length=n)
system.time({
pl <- numeric(n)
for(i in 1:n){
  pl[i] <- nlminb(c(200, 700, 200, 0), Pl.LA,sig.u=sig.u[i])$objective
  print(i)
}
})


plot(exp(sig.u),exp(-(pl-fit$objective)),type="l")
lines(range(exp(sig.u)),exp(-0.5*qchisq(0.95,df=1))*c(1,1))
# Spline approximation
sp.f <- splinefun(sig.u,-(pl-fit$objective))
x <- seq(min(sig.u),max(sig.u),length=200)
lines(exp(x),exp(sp.f(x)),lty=2,lwd=2)

# Compare with Wald confidence interval
round(
    exp(c(fit$par[5]-2*sqrt(diag(solve(H)))[5],
          fit$par[5]+2*sqrt(diag(solve(H)))[5])),
  digits=2)
##########################################################


#################################################################
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
## Over dispersion can be dealt with by hierarchical models (see lecture 12)
##########################################################
#
##########################################################
