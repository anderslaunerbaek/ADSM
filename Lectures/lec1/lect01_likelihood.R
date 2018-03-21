y <- 3 + rnorm(10,sd=2)
y

nll <- function(theta, y){
  -sum(dnorm(y,mean=theta[1],sd=sqrt(theta[2]), log=TRUE))
}

mean(y)
var(y)
nlminb(c(0,1),nll,y=y)

?nlminb
