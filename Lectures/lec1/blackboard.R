y <- 3 + rnorm(10, mean = 0, sd=2)
y

nll <- function(theta, y) {
  # negative log-likelihood
  -sum(dnorm(y, mean = theta[1], sd = sqrt(theta[2]), log = TRUE))
}

# choose optimizer
nlminb(start = c(0,1), objective = nll, y = y)

mean(y)

# the impirial variance are not the same as the estimated variance. There has been introduced some bias.
var(y)

