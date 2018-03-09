rm(list=ls())
setwd("~/Kurser/kursus02424/2018/lectures/lect03")
##################################################
#
##################################################
X <- cbind(1,rep(c(1,0),3),rep(c(0,1),3),c(1,1,rep(0,4)),
           c(0,0,1,1,0,0),c(rep(0,4),1,1))
X
solve(t(X)%*%X) # solve(X)=X^(-1)

X <- X[ ,-2]
solve(t(X)%*%X)
X

X <- X[ ,-3]
solve(t(X)%*%X)
X
## We can find the invese
##################################################



##################################################
## Class exercise..

## Original design
x1 <- c(1,0,1)
x2 <- c(0,1,1)
X <- cbind(x1,x2)

## New design
x1 <- c(1,1,2)
x2 <- c(0,1,1)

X2 <- cbind(x1,x2)
X2
X2%*%solve(t(X2)%*%X2)%*%t(X2)

X%*%solve(t(X)%*%X)%*%t(X) # So exactly the same model!
##################################################


##################################################
## Example
##################################################
rm(list=ls())
setwd("~/Kurser/kursus02424/2018/lectures/lect03")

ozone.pollution<-read.table("./ozone.data.csv",header=T)
names(ozone.pollution)
attach(ozone.pollution)

## A first view on data
pairs(ozone.pollution,panel=panel.smooth)

## gam to get an idea of relationsships
library(mgcv)
par(mfrow=c(2,2))
model<-gam(ozone~s(rad)+s(temp)+s(wind),data=ozone.pollution)
plot(model)

## Devellopment of a linear model
model0<-lm(ozone~temp*wind*rad+I(rad^2)+I(temp^2)++I(temp^3)+I(wind^2)+
               I(wind^3),data=ozone.pollution)
summary(model1)
par(mfrow=c(2,2))
plot(model0)



## Try transforming data
model0<-lm(log(ozone)~temp*wind*rad+I(rad^2)+I(temp^2)++I(temp^3)+I(wind^2)+
               I(wind^3),data=ozone.pollution)
summary(model0)
par(mfrow=c(2,2))
plot(model0)

## An odd looking data-point?
sort(ozone)
ozone.pollution[16:18, ]
n<-length(ozone)

pairs(ozone.pollution,panel=panel.smooth,pch=c(rep(1,16),19,rep(1,n-17)),
      col=c(rep(1,16),2,rep(1,n-17)),cex=c(rep(1,16),2,rep(1,n-17)))



## In log domain
pairs(cbind(ozone.pollution[,-4],log(ozone.pollution$ozone)),
      panel=panel.smooth,pch=c(rep(1,16),19,rep(1,n-17)),
      col=c(rep(1,16),2,rep(1,n-17)),cex=c(rep(1,16),2,rep(1,n-17)))


## Delete dat-point no. 17
model1<-lm(log(ozone)~temp*wind*rad+I(rad^2)+I(temp^2)+I(temp^3)+
               I(wind^2)+I(wind^3),data=ozone.pollution[-17, ])

summary(model1)
par(mfrow=c(2,2))
plot(model1)

## We should use a transformation between identity and log
model1<-lm(ozone~temp*wind*rad+I(rad^2)+I(temp^2)++I(temp^3)+I(wind^2)+
               I(wind^3),data=ozone.pollution[-17, ])
library(MASS)
par(mfrow=c(1,1))
boxcox(model1,lambda=seq(0,0.5,by=0.01))

model1<-lm(ozone^(1/4)~temp*wind*rad+I(rad^2)+I(temp^2)++I(temp^3)+I(wind^2)+
               I(wind^3),data=ozone.pollution[-17, ])
par(mfrow=c(2,2))
plot(model1)

summary(model1)

model2<-update(model1,~. -temp:wind:rad)
summary(model2)

model3<-update(model2,~. -wind:rad)
summary(model3)

model4<-update(model3,~. -I(wind^3))
summary(model4)

model5<-update(model4,~. -temp:rad)
summary(model5)

model6<-update(model5,~. - I(rad^2))
summary(model6)


anova(model6,model1)

## Could we reduce it further
model7<-update(model6,~. -temp:wind)
anova(model6,model7)
summary(model7)

##################################################
## Final model
model6
par(mfrow=c(2,2))
plot(model6)

## Can we compare likeliood across transformation? 
logLik(model6)
logLik(model0)
## No!!!!!!!! (at least not directly)!

## Example: comparing transfromations
model0<-lm(ozone~temp*wind*rad+I(rad^2)+I(temp^2)+I(temp^3)+
               I(wind^2)+I(wind^3),data=ozone.pollution[-17, ])
model01<-lm(ozone^0.25~temp*wind*rad+I(rad^2)+I(temp^2)+I(temp^3)+
               I(wind^2)+I(wind^3),data=ozone.pollution[-17, ])

logLik(model0)
logLik(model01)
n <- length(residuals(model01))
logLik(model01) - (sum(0.75*log(ozone.pollution$ozone[-17]))-n*log(0.25))

par(mfrow=c(1,1))
bc <- boxcox(model0,lambda=seq(0,1,by=0.01))
diff(range(bc$y))

##################################################
## Slides (R-code for plots in slides)
##################################################

##################################################
# Multivariate Normal
##################################################
rm(list=ls())
library(mvtnorm)
library(grid)
library(lattice)

Sigma0 <- diag(2)
Sigma1 <- diag(2)
Sigma1[1,2] <- Sigma1[2,1] <- 0.5
Sigma2 <- diag(2)
Sigma2[1,2] <- Sigma2[2,1] <- -0.5
Sigma3 <- diag(2)
Sigma3[1,2] <- Sigma3[2,1] <- 0.9
Sigma0
Sigma1
Sigma2
Sigma3


x <- seq(-4,4,by=0.05)
d <- expand.grid(x=x,y=x)
z0 <- dmvnorm(d,mean=c(0,0),sigma=Sigma0)
z0 <- z0/max(z0)
z1 <- dmvnorm(d,mean=c(0,0),sigma=Sigma1)
z1 <- z1/max(z1)
z2 <- dmvnorm(d,mean=c(0,0),sigma=Sigma2)
z2 <- z2/max(z2)
z3 <- dmvnorm(d,mean=c(0,0),sigma=Sigma3)
z3 <- z3/max(z3)

dat0 <- cbind(z=z0,d,rho=0)
dat1 <- cbind(z=z1,d,rho=0.5)
dat2 <- cbind(z=z2,d,rho=-0.5)
dat3 <- cbind(z=z3,d,rho=0.9)
dat <- rbind(dat0,dat1,dat2,dat3)


pdf("MultivarNormal.pdf",width=12/2.56,height=10/2.56)
levelplot(z ~ x + y|factor(rho), data=dat,layout=c(2,2),colorkey=FALSE,
          col.regions= rev(heat.colors(100))[c(1,10:100)],
          index.cond=list(c(1,4,2,3)),xlab="",ylab="",
          strip= strip.custom(factor.levels=c(expression(paste("C: ",rho==-0.5)),
                                              expression(paste("A: ",rho==0.00)),
                                              expression(paste("B: ",rho==0.5)),
                                              expression(paste("D: ",rho==0.9)))))
dev.off()



##################################################
## Norm of vectors

pdf("MultivarNormalNorms.pdf",width=11/2.56,height=9/2.56)
levelplot(z ~ x + y|factor(rho), data=dat,layout=c(2,2),colorkey=FALSE,
          col.regions= rev(heat.colors(100))[c(1,10:100)],
          index.cond=list(c(1,4,2,3)),xlab="",ylab="",
          strip= strip.custom(factor.levels=c(expression(paste("C: ",rho==-0.5)),
                                              expression(paste("A: ",rho==0.00)),
                                              expression(paste("B: ",rho==0.5)),
                                              expression(paste("D: ",rho==0.9)))),
          panel = function(...){
              panel.levelplot(...)
              panel.xyplot(c(0,1),
                           c(0,1),
                           type="b",
                           col=1,lwd=2)
          }

          )
dev.off()



y <- c(1,1)##/sqrt(2)
c(y%*%solve(Sigma0)%*%y,
  y%*%solve(Sigma1)%*%y,
  y%*%solve(Sigma2)%*%y,
  y%*%solve(Sigma3)%*%y)
##################################################


##################################################
## Projection
##################################################


##################################################
# Example 
##################################################
rm(list=ls())
library(scatterplot3d)
set.seed(123)
X <- matrix(c(1,0,0,1,1,1), ncol = 2, byrow=TRUE)
C <- t(X) %*% X
P <- X %*% solve(C) %*% t(X)
y <- numeric(3)
y[1] <- 3
y[2] <- 6
y[3] <- y[1] + y[2]
e <- 3*rnorm(3)
y1 <- y + e

y1%*%(diag(3)-P)%*%y1

y1%*%(P%*%y1)

##C1 <- t(X1) %*% X1
##P1 <- X1 %*% solve(C1) %*% t(X1)

pdf("ProjectionEks.pdf",width=13/2.56,height=13/2.56)
cex<-0.8
##par(mfrow=c(1,1),mar=c(5,5,1,5))
s3d <- scatterplot3d(x=y1[1],y=y1[2],z=y1[3],
                     pch=19,box=FALSE,
                     color="black",angle=30,
                     xlab=expression(y[1]),
                     ylab=expression(y[2]),
                     zlab=expression(y[3]),
                     xlim=c(0,10),
                     ylim=c(0,10),zlim=c(0,20))

yhat <- as.numeric(P%*%y1) # Projections
(y1-yhat)%*%yhat # Orthogonal?

sum((y1-yhat)*yhat)


for(i in seq(0,8,by=2)){
  y.tmp <- cbind(c(i,0,0),c(i,8,0))
  y.tmp[3, ] <- y.tmp[1, ]+y.tmp[2, ]
  yhat.tmp <- P%*%y.tmp
  s3d$points3d(x=yhat.tmp[1,],y=yhat.tmp[2,],
               z=yhat.tmp[3,],
               pch=19,type="l",col=gray(0.75))
}


for(i in seq(0,8,by=2)){
  y.tmp <- cbind(c(0,i,0),c(8,i,0))
  y.tmp[3, ] <- y.tmp[1, ]+y.tmp[2, ]
  yhat.tmp <- y.tmp
  s3d$points3d(x=yhat.tmp[1,],y=yhat.tmp[2,],
               z=yhat.tmp[3,],
               pch=19,type="l",col=gray(0.75))
}



s3d$points3d(x=c(y1[1],yhat[1]),
             y=c(y1[2],yhat[2]),
             z=c(y1[3],yhat[3]),
             pch=19,type="b",col=c(1,"blue"),lwd=1)

s3d$points3d(x=c(y1[1],0),
             y=c(y1[2],0),
             z=c(y1[3],0),
             pch=19,type="b",col=c(1,"blue"),lwd=2)

s3d$points3d(x=c(0,yhat[1]),
             y=c(0,yhat[2]),
             z=c(0,yhat[3]),
             pch=19,type="b",col=c(1,"blue"),lwd=2)


s3d.coords <- s3d$xyz.convert(y1[1]-1, y1[2], y1[3]+1)
text(s3d.coords$x, s3d.coords$y,             # x and y coordinates
         labels=expression(y),               # text to plot
         cex=cex, pos=4)

s3d.coords <- s3d$xyz.convert(yhat[1], yhat[2], yhat[3]+1)
text(s3d.coords$x, s3d.coords$y,             # x and y coordinates
         labels=expression(paste(hat(y))),               # text to plot
         cex=cex, pos=4)


s3d.coords <- s3d$xyz.convert(yhat[1]/2-0.5, yhat[2]/2, yhat[3]/2+4)
text(s3d.coords$x, s3d.coords$y,             # x and y coordinates
         labels=expression(paste("||",hat(y),"||")),               # text to plot
         cex=cex, pos=4)


s3d.coords <- s3d$xyz.convert(y1[1]/2-2, y1[2]/2, y1[3]/2+2)
text(s3d.coords$x, s3d.coords$y,             # x and y coordinates
         labels=expression(paste("||",y,"||")),               # text to plot
         cex=cex, pos=4)


##s3d.coords <- s3d$xyz.convert((yhat[1]+y1[1])/2-1, (yhat[1]+y1[2])/2-4, (yhat[1]+y1[3])/2+4)
s3d.coords <- s3d$xyz.convert(3, 3, 18)
text(s3d.coords$x, s3d.coords$y,             # x and y coordinates
         labels=expression(paste("||",y-hat(y),"||")),               # text to plot
         cex=cex, pos=4)
dev.off()





pdf("ProjectionEks2.pdf",width=13/2.56,height=13/2.56)
cex<-0.8
##par(mfrow=c(1,1),mar=c(5,5,1,5))
s3d <- scatterplot3d(x=y1[1],y=y1[2],z=y1[3],
                     pch=19,box=FALSE,
                     color="black",angle=30,
                     xlab=expression(y[1]),
                     ylab=expression(y[2]),
                     zlab=expression(y[3]),
                     xlim=c(0,10),
                     ylim=c(0,10),zlim=c(0,20))

yhat <- as.numeric(P%*%y1) # Projections
(y1-yhat)%*%yhat # Orthogonal?

sum((y1-yhat)*yhat)


for(i in seq(0,8,by=2)){
  y.tmp <- cbind(c(i,0,0),c(i,8,0))
  y.tmp[3, ] <- y.tmp[1, ]+y.tmp[2, ]
  yhat.tmp <- P%*%y.tmp
  s3d$points3d(x=yhat.tmp[1,],y=yhat.tmp[2,],
               z=yhat.tmp[3,],
               pch=19,type="l",col=gray(0.75))
}


for(i in seq(0,8,by=2)){
  y.tmp <- cbind(c(0,i,0),c(8,i,0))
  y.tmp[3, ] <- y.tmp[1, ]+y.tmp[2, ]
  yhat.tmp <- y.tmp
  s3d$points3d(x=yhat.tmp[1,],y=yhat.tmp[2,],
               z=yhat.tmp[3,],
               pch=19,type="l",col=gray(0.75))
}



s3d$points3d(x=c(y1[1],yhat[1]),
             y=c(y1[2],yhat[2]),
             z=c(y1[3],yhat[3]),
             pch=19,type="b",col=c(1,"blue"),lwd=1)

s3d$points3d(x=c(y1[1],0),
             y=c(y1[2],0),
             z=c(y1[3],0),
             pch=19,type="b",col=c(1,"blue"),lwd=2)

s3d$points3d(x=c(0,yhat[1]),
             y=c(0,yhat[2]),
             z=c(0,yhat[3]),
             pch=19,type="b",col=c(1,"blue"),lwd=2)


s3d.coords <- s3d$xyz.convert(y1[1]-1, y1[2], y1[3]+1)
text(s3d.coords$x, s3d.coords$y,             # x and y coordinates
         labels=expression(y),               # text to plot
         cex=cex, pos=4)

s3d.coords <- s3d$xyz.convert(yhat[1], yhat[2], yhat[3]+1)
text(s3d.coords$x, s3d.coords$y,             # x and y coordinates
         labels=expression(paste(hat(y))),               # text to plot
         cex=cex, pos=4)


s3d.coords <- s3d$xyz.convert(yhat[1]/2-0.5, yhat[2]/2, yhat[3]/2+4)
text(s3d.coords$x, s3d.coords$y,             # x and y coordinates
         labels=expression(paste("||",hat(y),"||")),               # text to plot
         cex=cex, pos=4)


s3d.coords <- s3d$xyz.convert(y1[1]/2-2, y1[2]/2, y1[3]/2+2)
text(s3d.coords$x, s3d.coords$y,             # x and y coordinates
         labels=expression(paste("||",y,"||")),               # text to plot
         cex=cex, pos=4)


##s3d.coords <- s3d$xyz.convert((yhat[1]+y1[1])/2-1, (yhat[1]+y1[2])/2-4, (yhat[1]+y1[3])/2+4)
s3d.coords <- s3d$xyz.convert(3, 3, 18)
text(s3d.coords$x, s3d.coords$y,             # x and y coordinates
         labels=expression(paste("||",y-hat(y),"||")),               # text to plot
     cex=cex, pos=4)

x1 <- c(1,1,2)
x2 <- c(0,1,1)

X2 <- cbind(x1,x2)
X2
X2%*%solve(t(X2)%*%X2)%*%t(X2)

X%*%solve(t(X)%*%X)%*%t(X) # So exactly the same model!
##################################################

##################################################
# Partioning
beta.hat <- solve(t(X)%*%X)%*%t(X)%*%y1;beta.hat
sigma.hat <- as.numeric(sqrt(t(y1)%*%(diag(3)-P)%*%y1));
                sigma.hat
V.beta.hat <- sigma.hat^2 * solve(t(X)%*%X);
   sqrt(diag(V.beta.hat))

summary(M<-lm(y1~ -1 + X)) # Note the -1


predict(M,se=TRUE)
P%*%y1
sqrt(diag(P)*sigma.hat^2) # Variance of hat(y)



##################################################
# Model reduction?
X2 <- c(1,1,2)
P2 <- X2%*%solve(t(X2)%*%X2)%*%t(X2)
P2

(P-P2)%*%P2 # P1 and P2 orthogonal

I<- diag(3)
sum(diag(diag(3)-P))
sum(diag(P-P2))
sum(diag(P2))

yhat2 <- P2%*%y1

s3d$points3d(x=c(y1[1],yhat2[1]),
             y=c(y1[2],yhat2[2]),
             z=c(y1[3],yhat2[3]),
             pch=19,type="b",col=c(1,"red"))

s3d$points3d(x=c(yhat[1],yhat2[1]),
             y=c(yhat[2],yhat2[2]),
             z=c(yhat[3],yhat2[3]),
             pch=19,type="l",col=c(1))


y.tmp1 <- c(0,0,0)
y.tmp2 <- c(8,8,16)
y.hat.tmp1 <- P2%*%y.tmp1
y.hat.tmp2 <- P2%*%y.tmp2

s3d$points3d(x=c(y.hat.tmp1[1],y.hat.tmp2[1]),
             y=c(y.hat.tmp1[2],y.hat.tmp2[2]),
             z=c(y.hat.tmp1[3],y.hat.tmp2[3]),
             pch=19,type="l",col=c(2),lwd=2)



dev.off()



##################################################
##
##################################################    
