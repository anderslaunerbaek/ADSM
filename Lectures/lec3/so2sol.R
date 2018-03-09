#
# Solution in R to Exercise 2 (Question 1 - 3)
#
# Comments / answers are found below
#

# Question 1

pollute<-read.table("./so2.csv",header=T)
attach(pollute)
names(pollute)

pairs(pollute,panel=panel.smooth)


cov(pollute)
cor(pollute)
library(corrplot)
corrplot(cor(pollute))


#  make sure that you have down-loaded the tree library from CRAN
#

library(tree)
model<-tree(Pollution~.,data=pollute)
plot(model)
text(model)
?tree



# Question 2

# Here we should be careful not estimating the full model with
# all interactions due to the rather low number of samples
#
# Hence we will limit our attention to 2-way interactions
#
# After some initial estimation with the purpose of investigating
# 2-way interactions we consider the following model
#
?formula

model1<-lm(Pollution~Temp+Industry+Population+Wind+Rain+Wet.days+
           Wind:Rain+Wind:Wet.days+Population:Wind+Temp:Rain)
summary(model1)
par(mfrow=c(2,2))
plot(model1)
library(MASS)
boxcox(model1)
library(car)
qqPlot(model1,reps=10000)
qqPlot(model1,simulate=FALSE)
?qqPlot

model1a<-lm(log(Pollution) ~ Temp + Industry + Population +
            Population:Industry + Wind + Rain + Wet.days +
            Wind:Rain + Wind:Wet.days + Population:Wind +
            Temp:Rain)
summary(model1a)
par(mfrow=c(2,2))
plot(model1a)
drop1(model1,test="F")

qqPlot(model1a,reps=10000)
qqPlot(model1a,simulate=FALSE)


# The most insignificiant 2-way interaction is Temp:Rain
model2<-update(model1a,~. -Temp:Rain)
summary(model2)
drop1(model2,test="F")

# Now the most insigficiant 2-way interaction is Population:Wind

model3<-update(model2,~. -Industry:Population)

summary(model3)
drop1(model3,test="F")

model4<-update(model3,~. -Population:Wind)
summary(model4)
drop1(model4,test="F")
plot(model4)

# This is rather good. But what about the higher-order interactions?
# One reasonable solution would be to fit 3-way terms for the
# variables that already appear in 2-way interactions:

model5<-update(model4,~. + Wind:Rain:Wet.days)
summary(model5)

drop1(model5,test="F")


# This seems to be a reasonable model

# Let us try some model validation (... however for this you have
# not yet learned all the details ...)

plot(model4) 


# Other way around modelling

modela<-lm(Pollution~Industry)
summary(modela)

#

modelb<-update(modela,~. +Population)
summary(modelb)

#

modelc<-update(modelb,~. +Population:Industry)
summary(modelc)

#non-significant

# 

modeld<-update(modelb,~. +Wet.days)
summary(modeld)

# non-significant


modele<-update(modelb,~. +Temp)
summary(modele)
