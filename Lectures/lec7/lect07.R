setwd("~/DTU/Courses/ADSM/Lectures/lec7")
##################################################
## 1: Seeds (binomial)
##################################################
rm(list=ls())
dat <- read.table("seed.dat", sep = ";",  header = TRUE)
head(dat)
str(dat)

dat$variety <- as.factor(dat$variety)
dat$root <- as.factor(dat$root)
dat$resp <- cbind(dat$y, (dat$n - dat$y))
fit1 <- glm(resp ~ variety * root,
            family = binomial,
            data = dat)
fit1

par(mfrow=c(2,2))
plot(fit1)
## No concerns here

# just to illustrate calculations....
w <- dat$n
X <- model.matrix(fit1)
eta <- predict(fit1,type="link")
mu <- predict(fit1,type="response")
Vmu <- mu*(1-mu)
W <- diag(1/(1+exp(-eta))^4*exp(-2*eta)*1/Vmu)*w
H <- sqrt(W)%*%X%*%solve(t(X)%*%W%*%X)%*%t(X)%*%sqrt(W)
h <- diag(H)
stdDevRes <- residuals(fit1,type="deviance")*sqrt(w)/(sqrt(1-diag(H)))

par(mfrow=c(2,2)) 
plot(fit1,which=1)
points(predict(fit1,type="link"),residuals(fit1,type="deviance"),
       pch=19)

plot(fit1,which=2)
points(qnorm((1:21-0.5)/21),sort(stdDevRes),
       pch=19)

plot(fit1,which=3)
points(predict(fit1,type="link"),sqrt(abs(stdDevRes)),
       pch=19)

plot(fit1,which=5)
rps <- (dat$y/dat$n-mu)/sqrt(mu*(1-mu)*(1-h)/dat$n)
points(h,rps,pch=19)

# Goodness of fit test
summary(fit1)
(pval <- 1 - pchisq(33.28, 17))
## So the assumption is rejected

par(mfrow=c(1,2))
resDev <- residuals(fit1,type="deviance")
plot(jitter(as.numeric(dat$variety), amount=0.1), resDev,
     xlab="Variety", ylab="Deviance residuals", cex=0.6,
     axes=FALSE)
box()
axis(1,label=c("O.a. 75", "O.a. 73" ),at=c(1,2))
axis(2)
plot(jitter(as.numeric(dat$root),
            amount=0.1), resDev, xlab= "Root" ,
ylab="Deviance residuals", cex=0.6, axes=FALSE)
box()
axis(1,label=c("Bean","Cucumber"),at=c(1,2))
axis(2)
## Doesn't seems to be systematic effects here..

# Including overdispersion
fit2 <- glm(resp ~ variety * root,
          family = quasibinomial, data = dat)
summary(fit2)

# JUST TO COMPARE THIS MODEL IS CONSIDERED WRONG HERE
confint(fit1)
confint(fit2) # Wider CI's

## Model reduction
drop1(fit2, test = "F")

fit3 <- glm(resp ~ variety + root,
            family = quasibinomial,
            data = dat)
drop1(fit3, test = "F")

fit4<-glm(resp ~ root,
          family = quasibinomial, data = dat)
drop1(fit4, test="F")
summary(fit4)

par <- coef(fit4)
par ## Germination more likeli on root2

std<-sqrt(diag(vcov(fit4)))
std

par+std%o%c(lower=-1,upper=1)*qt(0.975,19)

confint.default(fit4)
# same as above but with quantile qnorm(0.975)
confint(fit4) ## profile likelihood based

# Odds ratio
exp(coef(fit4)[2])
exp(confint(fit4)[2, ])



##################################################
## 2: Accident rates (Posison)
##################################################
dat <- data.frame(sex = c("F", "M"),
                  years = c(17.3,21.4),
                  y = c(175, 320))

fit1<-glm(y~offset(log(years))+sex,family=poisson,data=dat)
## Note included offset.
anova(fit1,test= "Chisq")

summary(fit1)
## More likely that elderly males are involved in accidents

## Odds ratio:
exp(coef(fit1)[2])
## or
(320/21.4)/(175/17.3)

####################################################################
## 3: Challenger (Binomial)
####################################################################
dat<-read.table("challenger_data.txt",header=TRUE)
dat
par(mfrow=c(1,2))
plot(dat$temp, dat$failed, xlab= ' Temperature ' ,
     ylab= ' No damaged (out of 6) ' )
plot(dat$pres, dat$failed, xlab= ' Pressure ' ,
     ylab= ' No damaged (out of 6) ' )

dat$resp<-cbind(dat$failed,dat$n-dat$failed)

fit0<-glm(formula = resp ~ temp+pres,
          family = binomial(link = logit),
          data = dat)
summary(fit0)

## Different parametrization
fit02<-glm(formula = I(failed/n) ~ temp+pres,
          family = binomial(link = logit),
          data = dat,weights=dat$n)
summary(fit02)
summary(fit0)

1-pchisq(9.4,df=20)
par(mfrow=c(2,2))
plot(fit0)
## many zeros (small samples) imply "bad" dianostic plots

drop1(fit0, test= "Chisq" )


fit1<-glm(formula = resp ~ temp,
          family = binomial(link = logit),
          data = dat)
drop1(fit1, test= "Chisq" )

summary(fit1)
1-pchisq(9.5,df=21) ## accept, but the chi^2 app. is not good in this case

par(mfrow=c(1,2))
tmp <- 31:85
pred<-predict(fit1, type= "response" , newdata=data.frame(temp=tmp),se=TRUE)
plot(tmp, pred$fit, type= "l" , lwd=3, col= "red" , xlab= "Temperature" , ylab= "P(damage)" )
points(dat$temp,dat$fail/dat$n)
lines(tmp,pred$fit+2*pred$se.fit)
lines(tmp,pred$fit-2*pred$se.fit)
coef(fit1)

logit <- function(x){  exp(x)/(1+exp(x))    }
pred<-predict(fit1, type= "link" , newdata=data.frame(temp=tmp),se=TRUE)
plot(tmp, logit(pred$fit), type= "l" , lwd=3, col= "red" , xlab= "Temperature" , ylab= "P(damage)" )
points(dat$temp,dat$fail/dat$n)
lines(tmp,logit(pred$fit+2*pred$se.fit))
lines(tmp,logit(pred$fit-2*pred$se.fit))
coef(fit1)
##################################################
##
##################################################



####################################################
## 4: Example 4.12 (Multinomial)
####################################################
rm(list=ls())
vdiss <- c( 234,  41,  42, 35)
diss  <- c( 559, 100,  76, 48)
neu   <- c(1157, 145,  89, 39)
sat   <- c(5826, 602, 254, 95)
vsat  <- c(2553, 237,  72, 27)

(tab <- cbind(vdiss, diss, neu, sat, vsat))
tot <- rowSums(tab)

tabp <- tab/tot
(accp <- t(apply(tabp,1,cumsum)))
t <- c(0,2,5,7)

par(mfrow=c(1,2))
matplot(t, accp[ ,-5], pch = 19,
        ylim = c(0,1), type = "b")

logisp <- log(accp[ ,-5]/(1 - accp[ ,-5]))

matplot(t, logisp[ ,-5], pch = 19,
        type = "b")

dat <- rbind(data.frame(y=vdiss,t=t,n=tot,level="vdiss"),
             data.frame(y=diss+vdiss,t=t,n=tot,level="diss"),
             data.frame(y=neu+diss+vdiss,t=t,n=tot,level="neu"),
             data.frame(y=sat+neu+diss+vdiss,t=t,n=tot,level="sat"))


## Model
dat$resp <- cbind(dat$y,dat$n-dat$y)
glmTot <- glm(resp~factor(t)*level,data=dat, family=binomial) # Satuated model
pred01 <- predict(glmTot,type="link",se=TRUE)

## Unceartainty on obs (using the satuated model)
levels <- levels(dat$level)
for(i in 1:length(t)){
  for(j in 1:length(levels)){
    I <- dat$level==levels[j] & dat$t==t[i]
    lines(t[i]*c(1,1),pred01$fit[I]+c(-2,2)*pred01$se.fit[I],col=j)
  }
}

## Two models
glm1 <- glm(resp~t*level, data=dat, family=binomial)
glm2 <- glm(resp~t+level, data=dat, family=binomial)

pred11 <- predict(glm1,type="link",se=TRUE)

matplot(t, logisp[ ,-5], pch = 19,
        type = "p")

for(i in 1:4){
  lines(t,pred11$fit[dat$level==levels[i]],col=i)
}


anova(glm1,glmTot,test="Chisq")
anova(glm2,glm1,test="Chisq") ## we need diff. slopes


drop1(glm1,test="Chisq") # I.e. different slopes
## These the deviances are not independent hence these are approximative.

summary(glm1)

par(mfrow=c(2,2))
plot(glm1)

anova(glm1,glm2,test="Chisq") # I.e. different slopes


1-pchisq(2.01,df=8)

## Plot the result
par(mfrow=c(1,1))
pred1 <- predict(glm1,type="response",se=TRUE)
matplot(t, accp[ ,-5], pch = "x",
        ylim = c(0,1), type = "p")

for(i in 1:4){
  lines(t,pred1$fit[dat$level==levels[i]],col=i)
}

## UNceartainties
for(i in 1:length(t)){
  for(j in 1:length(levels)){
    I <- dat$level==levels[j] & dat$t==t[i]
    lines(t[i]*c(1,1),pred1$fit[I]+c(-2,2)*pred1$se.fit[I],col=j)
  }
}

## Extrapolating
par(mfrow=c(1,2))
time <- 0:20
newdat<-data.frame(t=rep(time,4),level=rep(c("vdiss","diss","neu","sat"),each=length(time)))

pred <- predict(glm1,newdata=newdat,type="response",se=TRUE)
plot(time,pred$fit[newdat$level=="vdiss"],type="l",ylim=c(0,1))
lines(time,pred$fit[newdat$level=="diss"],lty=2,ylim=c(0,1))
lines(time,pred$fit[newdat$level=="neu"],lty=3,ylim=c(0,1))
lines(time,pred$fit[newdat$level=="sat"],lty=4,ylim=c(0,1))


time <- 0:20
newdat<-data.frame(t=rep(time,4),level=rep(c("vdiss","diss","neu","sat"),each=length(time)))

pred <- predict(glm1,newdata=newdat,type="link",se=TRUE)
plot(time,pred$fit[newdat$level=="vdiss"],type="l",ylim=range(pred$fit))
lines(time,pred$fit[newdat$level=="diss"],lty=2)
lines(time,pred$fit[newdat$level=="neu"],lty=3)
lines(time,pred$fit[newdat$level=="sat"],lty=4)



##################################################
## 5: Example 4.13 (Gamma)
##################################################
var <- c(1.29,1.223,0.8248)
nsheets <- c(3,5,7)
wgt <- c(2,4,6)/2

fit0 <- glm(var~-1+nsheets,family=Gamma(link=inverse),weights=wgt)
summary(fit0)

par(mfrow=c(2,2))
plot(fit0) ## useless with 3 obs

nshet <- 1:10
pred.link <- predict(fit0, type= "link" ,
                     newdata=data.frame(nsheets=nshet),se=TRUE)
pred.resp <-predict(fit0, type= "response" ,
                    newdata=data.frame(nsheets=nshet),se=TRUE)

pred.link
pred.resp



##################################################
## end
##################################################
