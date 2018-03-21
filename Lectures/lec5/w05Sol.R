##################################################
# Exercise 1
##################################################
y1 <- 4; y2 <- 14
## see table 4.2

par(mfrow=c(2,2))
# Normal
mu <- seq(0,30,by=0.01)
matplot(mu, cbind((mu-y1)^2,(mu-y2)^2),type="l",
        ylim=c(0,10))
rug(y1,col=1,lwd=2)
rug(y2,col=2,lwd=2)
# var = 1

# Poisson
d <- function(y,mu){2*(y*log(y/mu)-(y-mu))}
matplot(mu,cbind(d(y1,mu),d(y2,mu)),type="l",
        ylim=c(0,10))
rug(y1,col=1,lwd=2)
rug(y2,col=2,lwd=2)
# Var = 4, 14


# Gamma
d <- function(y,mu){2*(y/mu-log(y/mu)-1)}
matplot(mu,cbind(d(y1,mu),d(y2,mu)),type="l",
        ylim=c(0,2))
rug(y1,col=1,lwd=2)
rug(y2,col=2,lwd=2)
# Var = 4^2,14^2

# I Gaus
d <- function(y,mu){(y-mu)^2/(y*mu^2)}
matplot(mu,cbind(d(y1,mu),d(y2,mu)),type="l",
        ylim=c(0,2))
rug(y1,col=1,lwd=2)
rug(y2,col=2,lwd=2)
# Var = 4^3,14^3
##################################################

##################################################
## Exercise 2
##################################################
setwd("~/Kurser/kursus02424/2018/lectures/lect05/")
dat <- read.table("AIDS.csv",sep=",",header=TRUE)

plot(dat$cases) ## clear increase in time

plot(dat$year,dat$cases) ## clear increase in time (measured in years)

plot(dat$quarter,dat$cases) ## not a clear seasonal effect

dat <- cbind(dat,time=1:20)

## An initial model could be (this is a satuated model)
fit0 <- glm(cases~ factor(time),family=poisson,data=dat)
summary(fit0)

## can we reduce this model?

fit1 <- glm(cases~ time,family=poisson,data=dat)
summary(fit1)

AIC(fit0)
AIC(fit1)
## This reduction is too much

## Try higher order polynomials
fit2 <- glm(cases~ time + I((time-10)^2),family=poisson,data=dat)
summary(fit2)

AIC(fit0)
AIC(fit2)
## This model is prefered
## Plot model
par(mfrow=c(1,2))
plot(dat$time,coef(fit2)[1]+coef(fit2)[2]*dat$time+coef(fit2)[3]*(dat$time-10)^2,type="l")
## this might look like log(time)
plot(log(dat$time),type="l")

## Try log time
fit3 <- glm(cases~ log(time),family=poisson,data=dat)
summary(fit3)

AIC(fit0)
AIC(fit2)
AIC(fit3)
## Hence a much better fit, so final model fit3

## Hence final model is:

## log(lambda_i)=beta_0+beta_1*log(time_i) or
## lambda = e^(beta_0)*time_i^beta_1

## plot result:

pred <- predict(fit3,type="response")

par(mfrow=c(1,1))
plot(dat$time,dat$cases) ## clear increase in time
lines(dat$time,pred)
##################################################
