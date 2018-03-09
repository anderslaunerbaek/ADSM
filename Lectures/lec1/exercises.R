# exercise 1 ----
# Calculate the probability for each of the following events:
# a) A standard normally distributed variable is larger than 2.
pnorm(2, mean=0, sd=1, lower.tail = F)

# b) A normally distributed variable with mean 40 and variance equal to 9 is smaller than 34.
pnorm(34, mean=40, sd=9, lower.tail = F)

# c) Getting 9 successes out of 10 trials in a binomial experiment with p = 0.8.
pbinom(9, 10, 0.8)

# d) X > 6.2 in a Xi^2 distribution with 2 degrees of freedom.
pchisq(6.2, df=2)

# exercise 2 ----

x <- c(-1,0,1,2,3,4,5,6,7,8)
y <- c(1.4,4.7,5.1,8.3,9.0,14.5,14.0,13.4,19.2,18)
# fit model
lm(y ~ x)

# exercise 3 ----
# Use the following observations from a negative binomial distribution.
y <- c(13, 5, 28, 28, 15, 4, 13, 4, 10, 17, 11, 13, 12, 17, 3)
fun <- function(x) { (x[1] - 3)^2 + x[2]^2 }
fit <- optim(par = c(2, 2), fn = fun)
fit$par

# use loglike
nll <- function(theta, y) {
  # negative log-likelihood
  -sum(dnbinom(x = y, size = theta[1], prob = theta[2], log = TRUE))
}

nlminb(start = c(3,0.1), objective = nll, y = y)$par
