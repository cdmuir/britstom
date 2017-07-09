setwd("~/Google Drive/StomatalRatio/BES")
library(rstan)
library(MASS)
library(boot)
library(rstanarm)

# Simple linear regression (for understanding)
dat <- list(N = 1e2, K = 2, x = array(mvrnorm(1e2, c(0, 0), diag(2)), dim = c(100, 2)))
dat$y <- as.numeric(10 + dat$x %*% array(c(2, 3), dim = c(2, 1)) + rnorm(1e2))
m1 <- lm(dat$y ~ dat$x[, 1] + dat$x[, 2])
m1 <- stan(file = "Stan/linearRegression.stan", chains = 1, iter = 2e3, data = dat)

# Beta regression
# Two groups: group 1 ~ beta(phi = 0.25, lambda = 1); group 1 ~ beta(phi = 0.75, lambda = 1);
dat <- list(N = 1e3, K = 1)
dat$x <- array(c(rep(0:1, each = dat$N / 2)), dim = c(dat$N, dat$K))
phi <- 0.25 + dat$x %*% array(0.5, dim = c(dat$K, 1))
lambda <- 1
alpha <- lambda * phi
beta <- lambda * (1 - phi)
dat$y <- rbeta(dat$N, alpha, beta)
m1a <- lm(dat$y ~ dat$x)
summary(m1a)
m1b <- stan(file = "Stan/betaRegression.stan", chains = 1, iter = 2e4, data = dat)
summary(m1b)
plot(m1b, pars = c("phi[1]", "phi[501]", "lambda"))

# Truncated Beta regression
# One group: group 1 ~ beta(phi = 0.25, lambda = 1)
simPar <- list(N = 1e3, L = 0.01, U = 0.99, phi = 0.25, lambda = 1)
simPar$alpha <- simPar$lambda * simPar$phi
simPar$beta <- simPar$lambda * (1 - simPar$phi)
dat <- list(y = rbeta(simPar$N, simPar$alpha, simPar$beta))
dat$y_obs <- dat$y[dat$y > simPar$L & dat$y < simPar$U]
dat$N_obs <- length(dat$y_obs)
dat$N_censL <- length(which(dat$y <= simPar$L))
dat$N_censU <- length(which(dat$y >= simPar$U))
dat$L <- simPar$L
dat$U <- simPar$U

m2 <- stan(file = "Stan/truncBetaRegression.stan", chains = 1, iter = 2e3, data = dat)
summary(m2)
plot(m2, pars = c("phi", "lambda"))

# Normal mixture (for understanding)
simPar <- list(N = 3e2, mu = c(-5, 0, 5), sigma = c(1, 0.1, 2), theta = c(0.2, 0.2, 0.6))
dat <- list(y = rnorm(simPar$N, rep(simPar$mu, c(simPar$N * simPar$theta)),
                      rep(simPar$sigma, c(simPar$N * simPar$theta))))
dat$N <- simPar$N
dat$K <- 3

m3 <- stan(file = "Stan/normalMix.stan", chains = 1, iter = 2e3, data = dat)
summary(m3)

# Beta mixture
simPar <- list(N = 3e2, phi = c(0.05, 0.5, 0.95), lambda = c(10, 10, 10), theta = c(1/3, 1/3, 1/3))
simPar$alpha <- simPar$lambda * simPar$phi
simPar$beta <- simPar$lambda * (1 - simPar$phi)
dat <- list(y = rbeta(simPar$N, rep(simPar$alpha, c(simPar$N * simPar$theta)),
                      rep(simPar$beta, c(simPar$N * simPar$theta))))
dat$N <- simPar$N
dat$K <- 3

m4 <- stan(file = "Stan/betaMix.stan", chains = 1, iter = 2e4, data = dat)
summary(m4)


# Next:
# beta mixture
# truncated beta mixture
# truncated beta mixture regression





# Example of stan_betareg
N <- 200
x <- rnorm(N, 2, 1)
z <- rnorm(N, 2, 1)
mu <- binomial(link = "logit")$linkinv(1 + 0.2*x)
phi <- exp(1.5 + 0.4*z)
y <- rbeta(N, mu * phi, (1 - mu) * phi)
hist(y, col = "dark grey", border = FALSE, xlim = c(0,1))
fake_dat <- data.frame(y, x, z)

fit <- stan_betareg(y ~ x | z, data = fake_dat, 
                    link = "logit", link.phi = "log", 
                    chains = 1, iter = 2e3) # for speed
print(fit, digits = 2)
plot(fit)
pp_check(fit)
prior_summary(fit)







library(rethinking)
# lambda = 2
# phi = c(0.25, 0.75)
# alpha = c(0.5, 1.5)
# beta = c(1.5 0.5)
N <- 1e3
dat <- data.frame(X1 = rep(c(0, 1), each = N / 2),
                  X2 = rep(c(0, 1, 0, 1), each = N / 4))
dat$Y <- 1 + 2 * dat$X1 + 3 * dat$X2 + rnorm(N)
#dat$X1 <- as.factor(dat$X1)
#dat$X2 <- as.factor(dat$X2)
summary(fit <- lm(Y ~ X1 + X2, data = dat))
m1 <- map2stan(alist(
  Y ~ dnorm(mu, sigma),
  mu <- a + b1 * X1 + b2 * X2,
  a ~ dnorm(0, 10),
  b1 ~ dnorm(0, 10),
  b2 ~ dnorm(0, 10),
  sigma ~ dcauchy(0, 5)
), data = dat, iter = 2e3)

precis(m1, depth = 2)

m1 <- map2stan(alist(
  X ~ dbeta(alpha, beta),
  alpha <- lambda[Y] * phi,
  beta <- lambda[Y] * (1 - phi),
  phi ~ dbeta(1, 1),
  lambda[Y] ~ dlnorm(0.1, 1.5)
), data = dat, iter = 2e2)

phi <- 0.5; lambda <- 10

dat <- list(N = 1e3, X = rbeta(1e3, 1, 1))

m2 <- stan(file = "Stan/betaRegression1.stan",
           chains = 1, iter = 2e3, data = dat)
summary(m2)
plot(m2, pars = c("phi", "lambda", "alpha[1]", "beta[1]"))

###

N <- 1e2
dat <- list(Y = rbeta(N, rep(c(0.5, 1.5), each = N / 2), rep(c(1.5, 0.5), each = N / 2)), 
            X = rep(c(1, 2), each = N / 2), N = N, N_X = 2)
m3 <- stan(file = "Stan/betaRegression2.stan",
           chains = 1, iter = 2e3, data = dat)
summary(m3)
plot(m3, pars = c("phi[1]", "phi[2]", "lambda"))

###

N <- 1e2
dat <- list(X1 = rep(c(0, 1), each = N / 2), X2 = rep(c(0, 1, 0, 1), each = N / 4))
dat$Y <- rbeta(N, rep(c(0.5, 0.5, 0.5, 1.5), each = N / 4), rep(c(1.5, 1.5, 1.5, 0.5), each = N / 4))
dat$N <- N

m4 <- stan(file = "Stan/betaRegression3.stan",
           chains = 1, iter = 2e3, data = dat)
summary(m4)
plot(m4, pars = c("phi[1]", "phi[26]", "phi[26]", "phi[51]", "phi[76]", "lambda"))
stancode(m2)
