dat <- read.table("seed.dat", sep = ";",  header = TRUE)
dat

summary(fit1<-glm(cbind(y,n-y)~factor(variety)+factor(root),family=binomial,data=dat))
1-pchisq(39.7,df=18)
## overdispersion

##################################################
## Laplace approximation implemented by direct use
## numerical hessian

## design matrix
X <- matrix(0,ncol=3,nrow=dim(dat)[1])
X[ ,1]<-1
X[dat$variety==2,2]<-1
X[dat$root==2,3]<-1

library(numDeriv)
nll <- function(u,beta,sigma.u,X){
    eta <- X%*%beta + u
    mu <- 1/(1+exp(-eta))
    -sum(dbinom(dat$y,prob=mu,size=dat$n,log=TRUE)+dnorm(u,sd=sigma.u,log=TRUE))
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


system.time(fit2 <- nlminb(c(0,0,0,1),nll.LA,X=X))
-fit2$objective
logLik(fit1)
1-pchisq(2*(-fit2$objective-logLik(fit1)),df=1)
## Hence a significant improvement

##################################################
## improve by using independence of u's
nll.LA2 <- function(theta,X){
    beta <- theta[1:dim(X)[2]]
    sigma.u <- exp(theta[dim(X)[2]+1])
    est <- nlminb(rep(0,length(dat$y)),objective = nll, beta=beta, sigma.u=sigma.u,X=X)
    u <- est$par
    l.u <- est$objective
    fun.tmp <- function(ui,u,beta,sigma.u,X,i){
        u[i]<-ui
        nll(u,beta,sigma.u,X)        
    }
    H <- numeric(length(u))
    for(i in 1:length(u)){
        H[i] <- hessian(func = fun.tmp, x = u[i],u=u, beta = beta, sigma.u = sigma.u, X=X,i=i)}
    l.u + 0.5 * log(prod(H/(2*pi)))
}

system.time(fit3 <- nlminb(c(0,0,0,0),nll.LA2,X=X))
## so slightly faster 
-fit2$objective
-fit3$objective
## and same likelihood 
logLik(fit1)


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
        H[i] <- hessian(func = fun.tmp, x = u[i],u=u, beta = beta, sigma.u = sigma.u, X=X,i=i)}
    l.u + 0.5 * log(prod(H/(2*pi)))
}


system.time(fit4 <- nlminb(c(0,0,0,0),nll.LA3,X=X))
## So a large improvement in timing

-fit2$objective
-fit3$objective
-fit4$objective
logLik(fit1)


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
## Very large improvement in timing
fit5

H <- hessian(nll.LA4,x=fit5$par,X=X)
V <- solve(H)
se <- sqrt(diag(V))

summary(fit1)
## unceartainty of parameters
z <- fit5$par/se
round(cbind(fit5$par,se,z , 2*(1-pnorm(abs(z)))),digits=3)

## compare with model with overdispersion
summary(fit1.overD <- glm(cbind(y,n-y)~factor(variety)+factor(root),
                        family=quasibinomial,data=dat))
##################################################
## End
##################################################


