# ex 10
rm(list = ls())

library(dplyr)

dat <- read.table("./lec10/seed.dat",
                  header=TRUE, sep = ";") %>% 
    mutate(#plate = 1:n(),
           variety = as.integer(variety - 1),
           root = as.integer(root - 1))
# colnames(dat) <- c("extraxt", "seed", "r", "n","plate")

plot(y~n, data=dat, col=root, pch=19)
plot(y~n, data=dat, col=variety, pch=19)
 
X <- dat %>% select(-y) %>% as.matrix(.)
y <- dat %>% select(y) %>% as.matrix(.)

p <- dim(X)[2]

theta <- rep(1, p) %>% as.matrix(.)
Bi <- rnorm(n = N, 0, 1) %>% as.matrix(.)

pi <- X %*% theta + Bi

rbinom(N, dat$n, pi)


########








#
library(numDeriv)
l.LA <- function(th) {
    u.init <- rep(0, nlevels(Orange$Tree))
    obj <- function(u) nl(th, u, Orange$cir)
    est <- nlminb(u.init, obj)
    lval <- est$obj
    u <- est$par
    H <- hessian(obj, u)
    lval + 0.5 * log(det(H)) - length(u)/2 * log(2 * pi)
}

#
fit <- nlminb(c(300, 700, 200, 0, 0), l.LA)
