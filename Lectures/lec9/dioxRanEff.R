######################################################################################
dat <- read.table("~/Kurser/kursus02424/2018/Assignments/dioxin.csv",sep=",",header=TRUE)
library(nlme)
fit1 <- gls(log(DIOX)~PLANT+factor(TIME)+LAB+O2COR+
                NEFFEKT+SO2+HCL,data=na.omit(dat),
             method="ML")

fit2 <- gls(log(DIOX)~PLANT+factor(TIME)+LAB+O2COR+
                NEFFEKT+SO2+HCL,data=na.omit(dat),
             weights=varIdent(form=~1|LAB),method="ML")


## Some obs are nearly equal...
## Can be identified by obs. nr.
head(dat)
i <- 1
I <- numeric(57)
for(j in 2:dim(dat)[1]){
    if((dat$OBSERV[j]-dat$OBSERV[j-1])!=1){i <- i+1}
    I[j] <- i
}

dat<- cbind(obs=I,dat)
## e.g
dat[dat$obs==15, ]

fit3 <- lme(log(DIOX)~PLANT+factor(TIME)+LAB+O2COR+
                NEFFEKT+SO2+HCL,data=na.omit(dat), random = ~ 1 | factor(obs),
            weights = varIdent(form= ~ 1 | LAB),method="ML")


anova(fit1,fit2,fit3)

##################################################
