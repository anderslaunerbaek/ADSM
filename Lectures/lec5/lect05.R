# Code used in Lecture No. 6 - Motivating example with mice.
# (Note that the R code on the slides are code for direct input..)
#

mice <- data.frame(
 stillb=c(15, 17, 22, 38, 144),
 total=c(297, 242, 312, 299, 285),
 conc=c(0, 62.5, 125, 250, 500)
)

mice

mice$resp <- cbind(mice$stillb, mice$total-mice$stillb)

mice.glm <- glm(formula=resp ~ conc,
                family=binomial(link=logit),data=mice)

## Deviance table
anova(mice.glm,test="Chisq")
pchisq(253.33,df=1,lower.tail=FALSE)

## Marginal test and coefficients
summary(mice.glm)

# Plot in the linear domain. Logit transformed observatoins and corresponding
# linear dose response assay

## Plot from slides
p <- mice$stillb/mice$total
logit <- log(p/(1-p))
plot(mice$conc, logit, las=1, ylim=c(-3.5,0.1),
     xlab='Conc', ylab="Logit(still born fraction)")
abline(coef(mice.glm))

# Plot in the original scale - Observed fraction stillborn and
# corresponding fitted values under logistic regression for
# dose response assay

## Extrapolate results (for illustration)
xx<-seq(0,1500)
yy<-coef(mice.glm)[1]+coef(mice.glm)[2]*xx
yy2<- exp(yy)/(1+exp(yy))
mice.pred <- predict(mice.glm, newdata=data.frame(conc=xx),
                     type='link', se.fit=TRUE)

## Link pred
plot(xx,yy, las=1, ylim=range(yy), xlab='Concentration',
     ylab="Still born fraction",xlim=range(xx),type="l")
matlines(xx,cbind(mice.pred$fit+2*mice.pred$se.fit,
                  mice.pred$fit-2*mice.pred$se.fit),lty=2,col=2)
points(mice$conc,logit,pch=19)

## Response pred
plot(mice$conc,p, las=1, ylim=range(yy2), xlab='Concentration',
     ylab="Still born fraction",xlim=range(xx),pch=19)
lines(xx,exp(yy)/(1+exp(yy)))


# The fittet values plus standard error
mice.fit <- predict(mice.glm, type='response',
                    se.fit=TRUE, interval="conficence")
mice.fit$fit
mice.fit$se.fit
mice.fit$residual.scale

## Confidence limits 
conc <- seq(0,1500)
mice.fit <- predict(mice.glm, type='response',
                    newdata = data.frame(conc=conc),
                    se.fit=TRUE, interval="conficence")
matlines(conc, cbind(mice.fit$fit + 2 * mice.fit$se.fit,
                     mice.fit$fit - 2 * mice.fit$se.fit),
         col=2,lty=2)
##################################################


# The response residuals
residuals(mice.glm, type='response')

# The deviance residuals
residuals(mice.glm, type='deviance')

# The Pearson residuals
residuals(mice.glm, type='pearson')
########################################################
## End
########################################################n
