#################################################
## Example (Simulation)

## Simulate
set.seed(123456)
alt<-round(runif(100,0,4000))
tmt<-as.factor(rep(c(0,1), each=50))
sex<-as.factor(rep(c('M','F','M','F'), each=25))
X<-model.matrix(~-1+sex+tmt:alt)
y<-X%*%c(2.3,5,1/1000,3/1000)+rnorm(100,sd=.7)

## Plot
par(mfrow=c(1,1))
# pdf("lect05_simex.pdf",width=20/2.54,height=15/2.54)
plot(y~alt,pch=as.numeric(sex),
     col=rep(c(0,1), each=50)+1)
legend(x="topleft",
       legend=c("F:Placebo", "M:Placebo",
           "F:Treatment","M:Treatment"),
       col=c(1,1,2,2),pch=c(1,2))
# dev.off()

## Estimation
# Directly in R
fit0 <- lm(y~sex*tmt*alt)
par(mfrow=c(2,2))
plot(fit0)
summary(fit0)

# "Hand calculations" in R
x1 <- as.numeric(sex == "M")
x2 <- as.numeric(tmt == 1)
X <- cbind(1,x1,x2,alt,x1*x2,x1*alt,x2*alt,
           x1*x2*alt)
beta <- solve(t(X)%*%X)%*%t(X)%*%y
beta

## Simpler model
fit1 <- lm(y ~ sex + tmt + alt)

## Test comparing simple af large model
anova(fit1, fit0)
# No interactions is too simple, but a numnber of models
# between the two..
##################################################


##################################################
# Hypothesis chains

##################################################
# Type 1
anova(fit0)
fit1 <- lm(y ~ (sex + tmt + alt)^2)
anova(fit1,fit0)## Exactly the same as last line in type I

fit1 <- lm(y ~ sex * tmt * alt - sex:alt)
anova(fit1,fit0)## Different result than the type I

fit01 <- lm(y ~ (sex + tmt + alt)^2 - tmt:alt )
fit11 <- lm(y ~ (sex + tmt + alt)^2 - tmt:alt - sex:alt)
anova(fit11,fit01)

anova(fit0)## Compare with line 4. 
       ## (diff in p-value, due to different baseline)



fit11 <- lm(y ~ sex + tmt + alt + tmt:alt)
anova(fit11,fit0)
anova(fit11)
summary(fit11)
add1(fit11,~sex*tmt*alt,test="F") ## note oly interactions where lower 
                                  ## order included

## Type II and III
library(car)
Anova(fit0,type="II")
Anova(fit0,type="III") ## different result
##################################################

##################################################
# Type II
drop1(fit0,test="F")
fit1 <- update(fit0,.~.-sex:tmt:alt)

drop1(fit1,test="F")
fit2 <- update(fit1,.~.-sex:tmt)

drop1(fit2,test="F")
fit3 <- update(fit2,.~.-sex:alt)

drop1(fit3,test="F")

## Report result:
confint(fit3)

par(mfrow=c(2,2))
plot(fit3) ## Still reasonable

# Plot to interpret
par(mfrow=c(1,1))
plot(y~alt,pch=as.numeric(sex),
     col=rep(c(0,1), each=50)+1)
legend(x="topleft",
       legend=c("F:Placebo", "M:Placebo",
                "F:Treatment","M:Treatment"),
       col=c(1,1,2,2),pch=c(1,2))
lines(sort(alt),predict(fit3,newdata=data.frame(sex="F",tmt=factor(0),alt=sort(alt))))
lines(sort(alt),predict(fit3,newdata=data.frame(sex="F",tmt=factor(1),alt=sort(alt))),col=2)
lines(sort(alt),predict(fit3,newdata=data.frame(sex="M",tmt=factor(0),alt=sort(alt))),
      col=1)
lines(sort(alt),predict(fit3,newdata=data.frame(sex="M",tmt=factor(1),alt=sort(alt))),
      col=2)


## Forward/backward/both (use AIC, not p-values)
(fit.step <- step(fit0,~.+1))
(fit.back <- step(fit0,~+1,direction="back"))
(fitNull <- lm(y ~ 1))
(fit.forw <- step(fitNull,~ sex * tmt * alt,direction="forward"))
(fit.step <- step(fitNull,~.+sex * tmt * alt))
## In this case they all give the same result
##################################################


##################################################
# Another balanced design (orthogonal design)
set.seed(12413)
x1 <- factor(rep(rep(1:3,each=2),2))
x2 <- factor(rep(1:2,each=6))

(X <- model.matrix(~x1*x2))

## Simulated data
y <- X%*%c(1,1,0,0,0,0)+rnorm(12)


M <- lm(y~x1*x2)
X <- model.matrix(M)
t(X)%*%X ## Not off diagonal elements
anova(M)
Anova(M,type=2) 
Anova(M,type=3) 


## Orthorgonal design (parmetrization)
M <- lm(y~x1*x2,contrasts = list(x1=contr.helmert,x2=contr.helmert))
X <- model.matrix(M)
t(X)%*%X
anova(M)
Anova(M,type=2)
Anova(M,type=3)
##################################################

##################################################
## Collinearity and polynomial regression
setwd("~/DTU/Courses/ADSM/Lectures/lec4")
ozone.pollution<-read.table("../lec3/ozone.data.csv",header=T)
## pdf("ozone.pdf",width=20/2.54,height=13/2.54)
pairs(ozone.pollution, panel = panel.smooth, main = "Ozone data")
## dev.off()

## model from last week
model1 <-lm(log(ozone)~temp+wind+rad+I(temp^2)+I(temp^3)+I(wind^2),
            data=ozone.pollution[-17, ])
summary(model1)
X <- model.matrix(model1)
round(cov2cor(t(X)%*%X),digits=2) ## correlations close to zero
round(sort(abs(cov2cor(t(X)%*%X)))[-(43:49)],digits=2) # remove diagonal elements
plot(data.frame(X[ ,-1])) ## Strong collinearities

## subtract mean value
model2 <-lm(log(ozone)~I(temp-mean(temp))+I(wind-mean(wind))+
                I(rad-mean(rad))+I((temp-mean(temp))^2)+
                I((temp-mean(temp))^3)+I((wind-mean(wind))^2),
            data=ozone.pollution[-17, ])
X <- model.matrix(model2)
colnames(X) <- c("intercept","temp","wind","rad","temp2","temp3","wind2")
round(cov2cor(t(X)%*%X),digits=2)
round(sort(abs(cov2cor(t(X)%*%X)))[-(43:49)],digits=2)
plot(data.frame(X[ ,-1])) # not as strong
summary(model2)

## Legendre polynomials
leg1<-function(x){
    x <- 0.5*(x-mean(x))/sd(x)
}
    
leg2<-function(x){
    x <- 0.5*(x-mean(x))/sd(x)
    0.5*(3*x^2-1)}
leg3<-function(x){
    x <- 0.5*(x-mean(x))/sd(x)
    0.5*(5*x^3-3*x)}


model3 <-lm(log(ozone)~I(leg1(temp))+I(leg1(wind))+I(leg1(rad))+I(leg2(temp))+
                        I(leg3(temp))+I(leg2(wind)),data=ozone.pollution[-17, ])

X <- model.matrix(model3)
colnames(X) <- c("intercept","temp","wind","rad","temp2","temp3","wind2")
round(cov2cor(t(X)%*%X),digits=2)
round(sort(abs(cov2cor(t(X)%*%X)))[-(43:49)],digits=2)
plot(data.frame(X[ ,-1])) ## Even less collinaerity
summary(model3)

## Step does not really work here!


##################################################
# Another example (residual analysis)
library(MASS)
?trees
## pdf("tree.pdf",width=20/2.54,height=13/2.54)
pairs(trees, panel = panel.smooth, main = "trees data")
## dev.off()

M0 <- lm(Volume ~ Height*Girth +
         Height*I(Girth^2),data=trees)
summary(M0)
drop1(M0,test="F")

M1 <- update(M0,.~.-Height:I(Girth^2))
drop1(M1,test="F")

M2 <- update(M1,.~.-Height:Girth)
drop1(M2,test="F")

## A physical model?
Ma <- lm(Volume ~ -1 +
         Height:I(Girth^2),data=trees)

summary(Ma)
summary(M2)

## Ma nested with M0
anova(Ma,M0)

add1(Ma, ~ Height * I(Girth^2) + Height * Girth,
     test="F")

confint(Ma)
pi /(3 * 4 * 144) # Cone
pi /(4 * 144)  #Cylinder
## We do not accept any of these

################################################
# Residual analysis
par(mfrow=c(2,2))
plot(Ma)

# Confidence regions for qqplots...
library(car)
par(mfrow=c(1,1))
qqPlot(Ma,simulate=FALSE)

## Outliers?
range(sort(rstudent(Ma)))

##################################################
## an outlier
model1 <-lm(log(ozone)~temp+wind+rad+I(temp^2)+I(temp^3)+I(wind^2),
            data=ozone.pollution)
range(rstudent(model1))
range(rstandard(model1))
## which one is the outlier
which(abs(rstudent(model1))==max(abs(rstudent(model1))))
##################################################

#################################################
# Transformation
M0 <- lm(Volume ~ Height:I(Girth^2) + 
         log(Height) + log(Girth),
         data=trees)
par(mfrow=c(2,2))
plot(M0)
par(mfrow=c(1,1))
boxcox(M0) # Choose a transformation.

M0 <- lm(log(Volume) ~   Height:I(Girth^2) +
         log(Height)  + log(Girth),data=trees)

drop1(M0,test="F")
M1 <- update(M0,.~.-Height:I(Girth^2))


drop1(M1,test="F")

summary(M1)

# Residuals.. 
par(mfrow=c(2,2))
plot(M1)


##################################################
## Confidence reions:
## Inspect the model
## Units:
## Volume: [ft^3]
## Height: [ft]
## Grith : [Inches] (This is actually the diameter here)

## I.e. we may assume that (cylinder hypothesis)
## Volume = pi * Height * (Grith / 2 / 12)^2
##        = (pi /(4 * 144)) * Height * Grith^2
log(pi /(4 * 144))
confint(M1) # could be

# How about a cone?
# Volume = pi /3 * Height * (Grith / 2 / 12)^2
#        = (pi /(3 *4 * 144)) * Height * Grith^2
log(pi /(3 * 4 * 144))
confint(M1) # could be


##################################################
# Confidence region
##install.packages("ellipse")
par(mfrow=c(2,2))
library(ellipse)

plot(ellipse(M1,which=c(1,2)),type="l")
points(log((pi /(4 * 144))),1,pch=19) # Not a cylinder
points(log((pi /(3 * 4 * 144))),1,pch=19) ## Cone
lines(confint(M1)[1, ],mean(confint(M1)[2, ])*c(1,1))
lines(mean(confint(M1)[1, ])*c(1,1),confint(M1)[2, ]) 

plot(ellipse(M1,which=c(1,3)),type="l")
points(log((pi /(4 * 144))),2,pch=19) ## Cylinder
points(log((pi /(3 * 4 * 144))),2,pch=19) ## Cone
lines(confint(M1)[1, ],mean(confint(M1)[3, ])*c(1,1))
lines(mean(confint(M1)[1, ])*c(1,1),confint(M1)[3, ])


plot(ellipse(M1,which=c(2,3)),type="l")
points(1,2,pch=19)
lines(confint(M1)[2, ],mean(confint(M1)[3, ])*c(1,1))
lines(mean(confint(M1)[2, ])*c(1,1),confint(M1)[3, ])

## In conclusion a cone seems like a reasonable model for the
## volume

##################################################
# Tests 
X <- model.matrix(M1)
cov2cor(solve(t(X) %*% X))

betah <- coef(M1)
beta0 <- c(log((pi /(3 * 4 * 144))),1,2)
sigma <- summary(M1)$sigma

t(betah-beta0) %*% t(X) %*% X %*% (betah-beta0)
sigma ^ 2 * 3 * qf(0.95, 3, 28)

1-pf((t(betah-beta0) %*% t(X) %*% X %*% (betah - beta0)/(sigma ^ 2 * 3))[1,1],3,28)
# So not a pefect cone
##################################################


##################################################
# Prediction
##################################################

plot(log(trees))

par(mfrow=c(1,1))
Heights <- mean(trees$Height) # For given Height!
pred <- predict(M1,se=TRUE,
                newdata=data.frame(Girth=trees$Girth,
                    Height=Heights),interval="prediction")


conf <- predict(M1,se=TRUE,
                newdata=data.frame(Girth=trees$Girth,
                    Height=Heights),interval="confidence")


matplot(trees$Girth,pred$fit,type="l",col=c(1,2,2),
        lty=c(1,2,2))


matlines(trees$Girth,conf$fit,type="l",col=c(1,3,3),
        lty=c(1,2,2))


# original domain
matplot(trees$Girth,exp(pred$fit),type="l",col=c(1,2,2),
        lty=c(1,2,2))

matlines(trees$Girth,exp(conf$fit),type="l",col=c(1,3,3),
        lty=c(1,2,2))



######################################################
## Appendix Type I-II anova formulated as projection


#################################################
# Example (Simulation)
set.seed(123456)
alt<-round(runif(100,0,4000))
tmt<-as.factor(rep(c(0,1), each=50))
sex<-as.factor(rep(c('M','F','M','F'), each=25))
X <- model.matrix(~-1+sex+tmt:alt)
y <- X%*%c(2.3,5,1/1000,3/1000)+
    rnorm(100,sd=.7)

X0<- model.matrix(~tmt*sex*alt)
head(X0)

##############################################
## Type I table
M0 <- lm(y~tmt*sex*alt)

m0 <- dim(X0)[2]
n <- dim(X0)[1]

## Projection matrices
H0 <- X0%*%solve(t(X0)%*%X0)%*%t(X0)
X1 <- X0[ ,-dim(X0)[2]]
H1 <- X1%*%solve(t(X1)%*%X1)%*%t(X1)
X2 <- X1[ ,-dim(X1)[2]]
H2 <- X2%*%solve(t(X2)%*%X2)%*%t(X2)
X3 <- X2[ ,-dim(X2)[2]]
H3 <- X3%*%solve(t(X3)%*%X3)%*%t(X3)
X4 <- X3[ ,-dim(X3)[2]]
H4 <- X4%*%solve(t(X4)%*%X4)%*%t(X4)
X5 <- X4[ ,-dim(X4)[2]]
H5 <- X5%*%solve(t(X5)%*%X5)%*%t(X5)
X6 <- X5[ ,-dim(X5)[2]]
H6 <- X6%*%solve(t(X6)%*%X6)%*%t(X6)
X7 <- X6[ ,-dim(X6)[2]]
H7 <- X7%*%solve(t(X7)%*%X7)%*%t(X7)

## Deviance (RSS partionioning)
D0 <- t(y)%*%(diag(100)-H0)%*%y
D1 <- t(y)%*%(H0-H1)%*%y
D2 <- t(y)%*%(H1-H2)%*%y
D3 <- t(y)%*%(H2-H3)%*%y
D4 <- t(y)%*%(H3-H4)%*%y
D5 <- t(y)%*%(H4-H5)%*%y
D6 <- t(y)%*%(H5-H6)%*%y
D7 <- t(y)%*%(H6-H7)%*%y

Df <- c(rep(1,7),n-m0)
D <- c(D7,D6,D5,D4,D3,D2,D1,D0)
F <- D[-8]/(D[8]/(n-m0))
pv <- 1-pf(F,1,n-m0)
anova(M0)
round(cbind(Df,D,c(D/Df),c(F,NA),c(pv,NA)),digits=3)


##############################################
## Type II table

m0 <- dim(X0)[2]
n <- dim(X0)[1]

## 3rd order interaction
H0 <- X0%*%solve(t(X0)%*%X0)%*%t(X0)
X1 <- X0[ ,-8]

## 2nd order interactions
H1 <- X1%*%solve(t(X1)%*%X1)%*%t(X1)
X2 <- X0[ ,-c(7:8)]
H2 <- X2%*%solve(t(X2)%*%X2)%*%t(X2)
X3 <- X0[ ,-c(6,8)]
H3 <- X3%*%solve(t(X3)%*%X3)%*%t(X3)
X4 <- X0[ ,-c(5,8)]
H4 <- X4%*%solve(t(X4)%*%X4)%*%t(X4)

## Main effects
X05 <- X0[ ,-c(6:8)]## Alt
H05 <- X05%*%solve(t(X05)%*%X05)%*%t(X05)
X5 <- X0[ ,-c(4,6:8)]
H5 <- X5%*%solve(t(X5)%*%X5)%*%t(X5)

X06 <- X0[ ,-c(5,7:8)]## Sex
H06 <- X06%*%solve(t(X06)%*%X06)%*%t(X06)
X6 <- X0[ ,-c(3,5,7:8)]
H6 <- X6%*%solve(t(X6)%*%X6)%*%t(X6)

X07 <- X0[ ,-c(5:6,8)]## tmt
H07 <- X07%*%solve(t(X07)%*%X07)%*%t(X07)
X7 <- X0[ ,-c(2,5:6,8)]
H7 <- X7%*%solve(t(X7)%*%X7)%*%t(X7)


## Deviance (RSS partionioning)
D0 <- t(y)%*%(diag(100)-H0)%*%y
D1 <- t(y)%*%(H0-H1)%*%y
D2 <- t(y)%*%(H1-H2)%*%y
D3 <- t(y)%*%(H1-H3)%*%y
D4 <- t(y)%*%(H1-H4)%*%y
D5 <- t(y)%*%(H05-H5)%*%y
D6 <- t(y)%*%(H06-H6)%*%y
D07 <- t(y)%*%(diag(100)-H07)%*%y
D7 <- t(y)%*%(H07-H7)%*%y

## setting up the table
D <- c(D7,D6,D5,D4,D3,D2,D1,D0)
F <- c(D7,D6,D5,D4,D3,D2,D1)/(D0/(n-m0))
pv <- 1-pf(F,df1=1,df2=n-m0)

Anova(M0,type="II")
round(cbind(D,Df,c(F,NA),c(pv,NA)),digits=3)


##############################################
## Type III table

M0 <- lm(y~tmt*sex*alt)
X0 <- model.matrix(M0)
m0 <- dim(X0)[2]
n <- dim(X0)[1]


## Projection matrices
H0 <- X0%*%solve(t(X0)%*%X0)%*%t(X0)
X1 <- X0[ ,-8]
H1 <- X1%*%solve(t(X1)%*%X1)%*%t(X1)
X2 <- X0[ ,-7]
H2 <- X2%*%solve(t(X2)%*%X2)%*%t(X2)
X3 <- X0[ ,-6]
H3 <- X3%*%solve(t(X3)%*%X3)%*%t(X3)
X4 <- X0[ ,-5]
H4 <- X4%*%solve(t(X4)%*%X4)%*%t(X4)
X5 <- X0[ ,-4]
H5 <- X5%*%solve(t(X5)%*%X5)%*%t(X5)
X6 <- X0[ ,-3]
H6 <- X6%*%solve(t(X6)%*%X6)%*%t(X6)
X7 <- X0[ ,-2]
H7 <- X7%*%solve(t(X7)%*%X7)%*%t(X7)
X8 <- X0[ ,-1]
H8 <- X8%*%solve(t(X8)%*%X8)%*%t(X8)


## Deviance (RSS partionioning)
D0 <- t(y)%*%(diag(100)-H0)%*%y
D1 <- t(y)%*%(H0-H1)%*%y
D2 <- t(y)%*%(H0-H2)%*%y
D3 <- t(y)%*%(H0-H3)%*%y
D4 <- t(y)%*%(H0-H4)%*%y
D5 <- t(y)%*%(H0-H5)%*%y
D6 <- t(y)%*%(H0-H6)%*%y
D7 <- t(y)%*%(H0-H7)%*%y
D8 <- t(y)%*%(H0-H8)%*%y
D <- c(D8,D7,D6,D5,D4,D3,D2,D1,D0)

Df <- c(rep(1,8),n-m0)
F<-c(D8,D7,D6,D5,D4,D3,D2,D1)/(D0/(n-m0))
pv <- 1-pf(F,df1=1,df2=n-m0)

Anova(M0,type="III")
round(cbind(D,Df,c(F,NA),c(pv,NA)),digits=3)
##################################################
##
##################################################
