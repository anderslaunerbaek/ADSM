
#rm(list=ls())
library(numDeriv)

##################################################
# Data
##################################################

plate<-1:21
extract<-c(rep(0,11),rep(1,10))
seed<-c(rep(0,5),rep(1,6),rep(0,5),rep(1,5))
exse<-extract*seed
r<-c(10,23,23,26,17,5,53,55,32,46,10,8,10,8,23,0,3,22,15,
     32,3)
n<-c(39,62,81,51,39,6,74,72,51,79,13,16,30,28,45,4,12,41,
     30,51,7)

#source("~/tsmodels/graph.R")
#graph(width=14,height=9,type="eps",file="Germ")
plot(r/n,ylim=c(0,1),xlab="Plate",pch=19)
I<-1:21
lines(c(1,1)*(max(I[(seed+extract)==0])+0.5),c(-1,2),
      lty=2)
text(x=1,y=0.9,labels="Extract=0\nSeed=0",pos=4)
lines(c(1,1)*(max(I[seed==1 & extract==0])+0.5),c(-1,2),
      lty=2)
text(x=7,y=0.9,labels="Extract=0\nSeed=1",pos=4)
lines(c(1,1)*(max(I[seed==0 & extract==1])+0.5),c(-1,2),
      lty=2)
text(x=12,y=0.9,labels="Extract=1\nSeed=0",pos=4)
text(x=17,y=0.9,labels="Extract=1\nSeed=1",pos=4)
#dev.off()

##################################################
# numeric optimzation
##################################################
ll<-function(B,theta,x1,x2,x3,r,n){
  lin.pred<-theta["mu"]+x1*theta["beta1"]+
    x2*theta["beta2"]+x3*theta["beta3"]+B
  p<-1/(1+exp(-lin.pred))
  nll<--log(p)*r-log(1-p)*(n-r)+
    0.5*(log(theta["sig2"]^2)+B^2/theta["sig2"]^2)
  return(nll)
}

lTheta.num<-function(theta,B,extract,seed,exse,r,n,ll){
  # return Laplace approximation of negative log likelihood
  jnll<-0
  for(i in 1:length(B)){
   nll <- optim(0, ll, theta=theta, x1=extract[i],
                x2=seed[i], x3=exse[i], r=r[i],
                n=n[i], method="BFGS", hessian=TRUE) ## numerical app. of hessian
   jnll<-jnll+nll$value+0.5*log(nll$hessian)
  }
  return(jnll)
}

theta<-c(mu=0,beta1=0,beta2=0,beta3=0,sig2=1)
B<-numeric(length(n))

## Timing using numerical optimization and numerical approximation
## of hessian
system.time(
  OPnum<-nlminb(theta,lTheta.num,B=B,extract=extract,
                seed=seed,exse=exse,r=r,n=n,ll=ll))


##################################################
## log likelihood and derivatives (automatic differentiation in R))
##################################################
## log likelihood

## Differentialtion with R
## Condiional likelihood
(obsll<-"-log(%s)*r-log(1-%s)*(n-r)+0.5*(log(sig2^2)+
         B^2/sig2^2)")
## formula for p
(p<-"exp(mu+x1*beta1+x2*beta2+x3*beta3+B)/
         (1+exp(mu+x1*beta1+x2*beta2+x3*beta3+B))")

(obsll <- sprintf(obsll,p,p)) # replacing string (insert p)
(obsll <- parse(text=obsll))  # Conveting to expression

## derivative of ll w.r.t. random effects (as R-function)
(dobsll.B <- deriv(obsll,c("B"),
                  function(B, mu, beta1, beta2, beta3,
                           sig2, x1, x2, x3, r, n){}))

## See the output
dobsll.B(c(0,1),0,0,0,0,1,0,0,0,10,30)

# Hessian of ll w.r.t random effects
dobsll.BEx <- D(obsll,c("B"))
(Hobsll.B <- deriv(dobsll.BEx, c("B"),
                  function(B, mu, beta1, beta2, beta3,
                           sig2, x1, x2, x3, r,
                           n){}))


## Likelihood
lB<-function(B,mu,beta1,beta2,beta3,sig2,x1,x2,x3,r,n,
             dobsll.B){
  lin.pred<-mu+x1*beta1+x2*beta2+x3*beta3+B
  p<-1/(1+exp(-lin.pred))
  nll<--log(p)*r-log(1-p)*(n-r)+0.5*(log(sig2^2)+
                                     B^2/sig2^2)
  return(nll)
}

## Derivative
DlB<-function(B ,mu, beta1, beta2, beta3, sig2, x1,
              x2, x3, r, n, dobsll.B){
  # return derivative of negative log-likelihod of B |
  # fixed effects
    attr(dobsll.B(B, mu, beta1, beta2, beta3, sig2,
                  x1, x2,x3,r,n), "gradient")
}

## Hessian
HlB<-function(B ,mu, beta1, beta2, beta3, sig2, x1,
              x2, x3, r, n, dobsll.B){
    attr(Hobsll.B(B, mu, beta1, beta2, beta3, sig2,
                  x1, x2, x3, r, n), "gradient")
}

## Marginal lielihood (laplace)
lTheta <- function(theta, B, lB, DlB, extract, seed,
                   exse, r, n, dobsll.B){
  # return Laplace approximation of negative log
  # likelihood
    jnll<-0
    for(i in 1:length(B)){
        ll<-nlminb(0, lB, gradient = DlB,hessian=HlB,
                   mu = theta["mu"],
                   beta1 = theta["beta1"],
                   beta2 = theta["beta2"],
                   beta3 = theta["beta3"],
                   sig2 = theta["sig2"],
                   x1 = extract[i], x2 = seed[i],
                   x3 = exse[i], r = r[i],
                   n = n[i], dobsll.B = dobsll.B)
        jnll<-jnll+ll$objective
        B[i]<-ll$par
    }
    jnll <- jnll +
        0.5 * sum(log(HlB(B=B, mu = theta["mu"],
                          beta1 = theta["beta1"],
                          beta2 = theta["beta2"],
                          beta3 = theta["beta3"],
                          sig2 = theta["sig2"],
                          x1 = extract,x2 = seed,
                          x3 = exse,r = r,
                          n = n)))
    return(jnll)
}


##################################################
## Estimate parameters
##################################################

theta<-c(mu=0,beta1=0,beta2=0,beta3=0,sig2=1)
B<-numeric(length(n))

system.time(OP1<-nlminb(theta,lTheta,B=B,lB=lB,DlB=DlB,
                        extract=extract,seed=seed,exse=exse,
                        r=r,n=n,dobsll.B=dobsll.B,
                        lower=c(-Inf,-Inf,-Inf,-Inf,0)))


OPnum$objective-OP1$objective
OPnum$par-OP1$par


H<-hessian(lTheta,OP1$par, B = B, lB = lB,
           DlB = DlB, extract = extract,
           seed=seed, exse = exse, r = r,
           n = n, dobsll.B = dobsll.B)

SE <- sqrt(diag(solve(H)))

Par <- cbind(OP1$par, SE, OP1$par - 2 * SE,
             OP1$par + 2 * SE)
## Summary statistics based on Wald statistics
round(Par,digits=3)
## all pars different from 0


##################################################
## present results
##################################################

theta<-OP1$par
B<-numeric(length(n))
for(i in 1:length(B)){
  ll <- nlminb(0, lB, gradient = DlB,
               mu = theta["mu"],
               beta1 = theta["beta1"],
               beta2 = theta["beta2"],
               beta3 = theta["beta3"],
               sig2 = theta["sig2"],
               x1 = extract[i], x2 = seed[i],
               x3 = exse[i], r = r[i], n = n[i],
               dobsll.B = dobsll.B)
  B[i]<-ll$par
}

qqnorm(B)
qqline(B)

#graph(width=14,height=12,type="eps",file="GermModel")
plot(plate,r/n,ylim=c(0,1))
phat<-exp(theta["mu"]+theta["beta1"]*extract+
          theta["beta2"]*seed+theta["beta3"]*exse)/
  (1+exp(theta["mu"]+theta["beta1"]*extract+
         theta["beta2"]*seed+theta["beta3"]*exse))
phat2<-exp(theta["mu"]+theta["beta1"]*extract+
           theta["beta2"]*seed+theta["beta3"]*exse+B)/
  (1+exp(theta["mu"]+theta["beta1"]*extract+
         theta["beta2"]*seed+theta["beta3"]*exse+B))
points(plate,phat,col="red",cex=2,pch=19)
points(plate,phat2,col="blue",pch=19)
points(plate,r/n,pch=19)
#dev.off()
n[16]

##################################################
# Profile likelihood
##################################################

## proofile likelihood for parameter theta
lThetaP<-function(theta,thetaP,B,lB,DlB,extract,seed,
                  exse,r,n,dobsll.B){
  # return Laplace approximation of negative log likelihood
  jnll<-0
  theta<-c(theta,thetaP)
  for(i in 1:length(B)){
    ll<-nlminb(0,lB,gradient=DlB,mu=theta["mu"],
               beta1=theta["beta1"],beta2=theta["beta2"],
               beta3=theta["beta3"],sig2=theta["sig2"],
               x1=extract[i],x2=seed[i],x3=exse[i],r=r[i],
               n=n[i],dobsll.B=dobsll.B)
    jnll<-jnll+ll$objective
    B[i]<-ll$par
  }
  jnll<-jnll+0.5*sum(log(attr(Hobsll.B(mu=theta["mu"],
               beta1=theta["beta1"],beta2=theta["beta2"],
               beta3=theta["beta3"],sig2=theta["sig2"],
               B=B,x1=extract,x2=seed,x3=exse,r=r,n=n),
                              "grad")))
  return(jnll)
}

## mu
thetaP<-c()
theta<-OP1$par
mu.vec<-sort(c(Par["mu",1],
               seq(Par["mu",1]-4*Par["mu",2],
                   Par["mu",1]+4*Par["mu",2],length=7)))
theta<-theta[-1]
PL.mu<-numeric(length(mu.vec))
for(i in 1:length(mu.vec)){
  thetaP["mu"]<-mu.vec[i]
  PL.mu[i]<-nlminb(theta,lThetaP,thetaP=thetaP,B=B,
                   lB=lB,DlB=DlB,extract=extract,seed=seed,
                   exse=exse,r=r,n=n,
                   dobsll.B=dobsll.B)$objective
  print(i)
}

## spline approximation
S.mu<-spline(mu.vec,exp(-PL.mu+OP1$objective),n=100)## app profile likelihood

## beta 1
thetaP<-c()
theta<-OP1$par
beta1.vec<-sort(c(Par["beta1",1],
               seq(Par["beta1",1]-4*Par["beta1",2],
                   Par["beta1",1]+
                   4*Par["beta1",2],length=8)))
theta<-theta[-2]
PL.beta1<-numeric(length(beta1.vec))
for(i in 1:length(beta1.vec)){
  thetaP["beta1"]<-beta1.vec[i]
  PL.beta1[i]<-nlminb(theta,lThetaP,thetaP=thetaP,
                   B=B,lB=lB,DlB=DlB,extract=extract,
                   seed=seed,exse=exse,r=r,n=n,
                      dobsll.B=dobsll.B)$objective
  print(i)
}

## spline approximation
S.beta1<-spline(beta1.vec,exp(-PL.beta1+OP1$objective),
                n=100)


## beta 2
thetaP<-c()
theta<-OP1$par
beta2.vec<-sort(c(Par["beta2",1],
               seq(Par["beta2",1]-4*Par["beta2",2],
                   Par["beta2",1]+4*Par["beta2",2],
                   length=8)))
theta<-theta[-3]
PL.beta2<-numeric(length(beta2.vec))
for(i in 1:length(beta2.vec)){
  thetaP["beta2"]<-beta2.vec[i]
  PL.beta2[i]<-nlminb(theta,lThetaP,thetaP=thetaP,B=B,
                      lB=lB,DlB=DlB,extract=extract,
                      seed=seed,exse=exse,r=r,n=n,
                   dobsll.B=dobsll.B)$objective
  print(i)
}

## spline approximation
S.beta2<-spline(beta2.vec,exp(-PL.beta2+OP1$objective),
                 n=100)

## beta 3
thetaP<-c()
theta<-OP1$par
beta3.vec<-sort(c(Par["beta3",1],
               seq(Par["beta3",1]-4*Par["beta3",2],
                   Par["beta3",1]+4*Par["beta3",2],
                   length=8)))
theta<-theta[-4]
PL.beta3<-numeric(length(beta3.vec))
for(i in 1:length(beta3.vec)){
  thetaP["beta3"]<-beta3.vec[i]
  PL.beta3[i]<-nlminb(theta,lThetaP,thetaP=thetaP,B=B,
                      lB=lB,DlB=DlB,extract=extract,
                      seed=seed,exse=exse,r=r,n=n,
                      dobsll.B=dobsll.B)$objective
  print(i)
}

## spline approximation
S.beta3<-spline(beta3.vec,exp(-PL.beta3+OP1$objective),
                n=100)

## sig 2
thetaP<-c()
theta<-OP1$par
sig2.vec<-sort(seq(1E-2,Par["sig2",1]+4*Par["sig2",2],length=8))
theta<-theta[-5]
PL.sig2<-numeric(length(sig2.vec))
for(i in 1:length(sig2.vec)){
  thetaP["sig2"]<-sig2.vec[i]
  PL.sig2[i]<-nlminb(theta,lThetaP,thetaP=thetaP,B=B,
                     lB=lB,DlB=DlB,extract=extract,
                     seed=seed,exse=exse,r=r,n=n,
                     dobsll.B=dobsll.B)$objective
  print(i)
}

## spline approximation
S.sig2<-spline(sig2.vec,exp(-PL.sig2+OP1$objective),n=100)


## Fixed effects
#graph(width=14,height=11,type="eps",file="ProfileFixed")
par(mfrow=c(2,2),mar=c(3,1,0,0),oma=c(0,2,1,1))
plot(mu.vec,exp(-PL.mu+OP1$objective),axes=FALSE)
lines(S.mu$x,S.mu$y)
lines(range(S.mu$x),c(1,1)*0.1465)
lines(range(S.mu$x),c(1,1)*0.03625)
lines((Par["mu",1]-2*Par["mu",2])*c(1,1),c(-1,2),
      col="blue")
lines((Par["mu",1]+2*Par["mu",2])*c(1,1),c(-1,2),
      col="blue")
axis(1);axis(2);box()
mtext(expression(mu),side=1,line=2)
mtext("Relative profile likelihood",side=2,line=2)

plot(beta1.vec,exp(-PL.beta1+OP1$objective),axes=FALSE)
lines(S.beta1$x,S.beta1$y)
lines(range(S.beta1$x),c(1,1)*0.1465)
lines(range(S.beta1$x),c(1,1)*0.03625)
lines((Par["beta1",1]-2*Par["beta1",2])*c(1,1),c(-1,2),
      col="blue")
lines((Par["beta1",1]+2*Par["beta1",2])*c(1,1),c(-1,2),
      col="blue")
mtext(expression(beta[1]),side=1,line=2)
axis(1);box()

plot(beta2.vec,exp(-PL.beta2+OP1$objective),axes=FALSE)
lines(S.beta2$x,S.beta2$y)
lines(range(S.beta2$x),c(1,1)*0.1465)
lines(range(S.beta2$x),c(1,1)*0.03625)
lines((Par["beta2",1]-2*Par["beta2",2])*c(1,1),c(-1,2),
      col="blue")
lines((Par["beta2",1]+2*Par["beta2",2])*c(1,1),c(-1,2),
      col="blue")
mtext(expression(beta[2]),side=1,line=2)
axis(1);box();axis(2)
mtext("Relative profile likelihood",side=2,line=2)

plot(beta3.vec,exp(-PL.beta3+OP1$objective),axes=FALSE)
lines(S.beta3$x,S.beta3$y)
lines(range(S.beta3$x),c(1,1)*0.1465)
lines(range(S.beta3$x),c(1,1)*0.03625)
lines((Par["beta3",1]-2*Par["beta3",2])*c(1,1),c(-1,2),
      col="blue")
lines((Par["beta3",1]+2*Par["beta3",2])*c(1,1),c(-1,2),
      col="blue")
mtext(expression(beta[3]),side=1,line=2)
axis(1);box()
#dev.off()


## random effects
## graph(width=14,height=11,type="eps",file="ProfileRandom")
par(mfrow=c(1,1),mar=c(3,1,0,0),oma=c(0,2,1,1))
plot(sig2.vec,exp(-PL.sig2+OP1$objective),axes=FALSE,
     xlab="",ylab="")
lines(S.sig2$x,S.sig2$y)
lines(range(S.sig2$x),c(1,1)*0.1465)
lines(range(S.sig2$x),c(1,1)*0.03625)
lines((Par["sig2",1]-2*Par["sig2",2])*c(1,1),c(-1,2),
      col="blue")
lines((Par["sig2",1]+2*Par["sig2",2])*c(1,1),c(-1,2),
      col="blue")
mtext(expression(sigma),side=1,line=2)
axis(1);box();axis(2)
mtext("Relative profile likelihood",side=2,line=2)
#dev.off()
##################################################


##################################################
## Compare with and without random effect 
library(glmmTMB)

fit1 <- glm(cbind(r,n-r)~factor(extract)*factor(seed),family="binomial")

fit2 <- glmmTMB(cbind(r,n-r)~factor(extract)*factor(seed)+(1|plate),family="binomial")

## Same pars 
OP1$par[1:4]-fit2$fit$par[1:4]
OP1$par[5]-exp(fit2$fit$par[5])## Different scale for random eff par


(ll1 <- logLik(fit1))
(ll2 <- logLik(fit2))
OP1$objective ## why is this different 


1-pchisq(2*(ll2-ll1),df=1) ## So random effect not significant


##################################################
## Hierachical model
##################################################

## adapted version of code from lecture 12
ldCBb <- function(y,alpha,beta,n){
  lchoose(n,y)+lgamma(alpha + y) + lgamma(beta + n - y) -
    lgamma(alpha + beta + n) - lgamma(alpha) - 
    lgamma(beta)+ lgamma(alpha+beta)}


CBinBetaReg <- function(pars,y,n,X){
    N <- length(y)
    mu <- 1/(1+exp(-( X%*%pars[-length(pars)])))
    phi <- exp(pars[length(pars)])
    alpha <- mu*phi
    beta <- (1-mu)*phi
    -sum(ldCBb(y,alpha,beta,n))
}

X <- model.matrix(fit2)

## Optimize log-likelihood
OPH <- nlminb(c(0,0,0,0,0), CBinBetaReg, y=r, n=n,X=X)


logLik(fit1)
logLik(fit2)
-OPH$objective


## Unconditional mean and precision
mu <- 1/(1+exp(-(X%*%OPH$par[-5])))
(phi <- exp(OPH$par[5])) ## large value

alpha <- mu * phi
beta <- (1-mu) * phi

## Conditional mean and precision (Theorem 6.6)
alpha.z <- alpha + r
beta.z <- beta + n - r
mu.z <- alpha.z/(alpha.z+beta.z)
sigma.z <- alpha.z+beta.z


## compare the random effects
mu.tmb <- 1/(1+exp(-X%*%fit2$fit$parfull[1:4]-fit2$fit$parfull[5:25]))


plot(mu.z,mu.tmb)
lines(c(0,1),c(0,1))
## So in the case the model are the same (at least for the mean values).

##################################################
## End
##################################################
