##################################################
# Ex 3.4
##################################################
x <- rep(1:3,2)
method <- c(rep(0,3),rep(1,3))
y <- c(0.22,0.38,0.72,0.31,0.66,0.99)

plot(y ~ x, col=method+1, pch = method+1)

v <- x^((method<=1)+(method==0))
Sigma <- diag(v)

# Parameter estimates
X0 <- cbind(1, method, x, method * x)
ISigma <- solve(Sigma)
beta0 <- solve(t(X0) %*% ISigma %*% X0) %*%
    t(X0) %*% ISigma %*% y

# Estimated (hat-)values
yh <- X0%*%beta0
sigma.h2 <- t(y-yh)%*%ISigma%*%(y-yh)/
    (6-4)
sigma.h2

# With R (lm)
M0 <- lm(y ~ factor(method)*x,weights=1/v)
beta0
coef(M0)

# Simplify model
M1 <- lm(y ~ factor(method):x,weights=1/v)
M2 <- lm(y ~ x,weights=1/v) 

# And test
anova(M1,M0)
anova(M2,M0)
# The is only weak evidence of difference  between the two methods. I.e. we will
# accept that the concentrations are the same.


##################################################
# Exercise 3.7
##################################################
data(anscombe)
M1 <- lm(y1 ~ x1, data=anscombe)
M2 <- lm(y2 ~ x2, data=anscombe)
M3 <- lm(y3 ~ x3, data=anscombe)
M4 <- lm(y4 ~ x4, data = anscombe)

summary(M1)
summary(M2)
summary(M3)
summary(M4)
# All stats are the same

# but have a look at residual plots
par(mfrow=c(2,2))
plot(M1) # ok.

par(mfrow=c(2,2))
plot(M2) # Not ok, missing a second order term
plot(anscombe$x2,residuals(M2))
M22 <- update(M2,.~.+I(x2^2))
par(mfrow=c(2,2))
plot(M22)

par(mfrow=c(2,2))
plot(M3) # Not ok. one outlier
plot(anscombe$x3,anscombe$y3)

par(mfrow=c(2,2))
plot(M4) # Not ok, no exertation of covariate 
plot(y4~x4,data = anscombe)


plot(y1~x1,data = anscombe)
plot(y2~x2,data = anscombe)
plot(y3~x3,data = anscombe)
plot(y4~x4,data = anscombe)
##################################################


## 3.6 (verification by simulation)
n <- 10
x <- seq(0.1,10,length=n)
sigma<-2
y <- rnorm(n,mean=x,sd=sigma/x)
plot(x,y)

## Q1:
(sum(x^2))^(-1)*sum(x*y)
summary(lm(y~-1+x))

## var of estimator
sigma^2/sum(x^2)^2

k <- 1000000
beta <- numeric(k)
for(i in 1:k){
    y <- rnorm(n,mean=x,sd=sigma/x)
    beta[i]<- (sum(x^2))^(-1)*sum(x*y)
}
mean(beta)
var(beta)
sigma^2*n/sum(x^2)^2
hist(beta)



## Q2:
Sigma <- diag(1/x^2)
sum(x^4)^(-1)*sum(x^3*y)

summary(lm(y~-1+x,weights=1/diag(Sigma)))

#k <- 100000
beta2 <- numeric(k)
for(i in 1:k){
    y <- rnorm(n,mean=x,sd=sigma/x)
    beta2[i]<- sum(x^4)^(-1)*sum(x^3*y)
}
mean(beta2)
var(beta2)-
sigma^2/sum(x^4)

hist(beta2)

var(beta2)/var(beta)
(sum(x^2)^2/n)/(sum(x^4))


## Q4:
library(mvtnorm)
rho<-0.5
Sigma <- diag(10)
Sigma1 <- (row(Sigma) == col(Sigma) + 1) 
Sigma[Sigma1] <- rho
Sigma1 <- (row(Sigma) +1== col(Sigma) ) 
Sigma[Sigma1] <- rho
Sigma

y <- as.numeric(rmvnorm(1,mean=x,sigma=sigma^2*Sigma))

lm(y~-1+x)

k <- 10000
beta <- numeric(k)
for(i in 1:k){
    y <- as.numeric(rmvnorm(1,mean=x,sigma=sigma^2*Sigma))
    beta[i]<- (sum(x^2))^(-1)*sum(x*y)
}

mean(beta)
var(beta)
sigma^2/sum(x^2)^2*x%*%Sigma%*%x
##################################################
##
##################################################
