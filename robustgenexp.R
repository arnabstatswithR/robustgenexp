
#-----------------------
# distribution function of the Generalized exponential (GE) distribution

pge <- function(x, lambda, nu){(1 - exp(-lambda * x))^nu}

# density function of the Generalized exponential (GE) distribution

dge <- function(x, lambda, nu){
  lambda * nu * exp(-lambda * x) * (1 - exp(-lambda * x))^(nu - 1)}

# quantile function of the Generalized exponential (GE) distribution

qge <- function(p, lambda, nu){-log(1 - p^(1 / nu)) / lambda}

# simulation function of the Generalized exponential (GE) distribution

rge <- function(n, lambda, nu){-log(1 - runif(n)^(1 / nu)) / lambda}

#------------------------
# estimation approaches

# ML estimation

nll <- function(theta, X){
  lambda <- exp(theta[1])
  nu <- exp(theta[2])
  n <- length(X)
  out <- -n * log(lambda) - n * log(nu) + lambda * sum(X) - 
    (nu - 1) * sum(log(1 - exp(-lambda * X)))
  out}

mle_ge <- function(X){
  nll.out <- optim(c(0, 0), nll, X = X, hessian = F)
  exp(nll.out$par)}

# MM estimation

library(pracma)

sqdiff.mme <- function(theta, X){
  lambda <- exp(theta[1])
  nu <- exp(theta[2])
  out <- (sqrt(psi(1, 1) - psi(1, nu + 1)) / (psi(0, nu + 1) - psi(0, 1)) - sd(X) / mean(X))^2 +
    (lambda - (psi(0, nu + 1) - psi(0, 1)) / mean(X))^2
  out}

mme_ge <- function(X){
  sqdiff.out <- optim(c(0, 0), sqdiff.mme, X = X, hessian = F)
  exp(sqdiff.out$par)}

# Percentile estimation

sqdiff.pt <- function(theta, X){
  lambda <- exp(theta[1])
  nu <- exp(theta[2])
  n <- length(X)
  X.sort <- sort(X)
  out <- sum((X.sort - qge(1:n / (n + 1), lambda, nu))^2)
  out}

pt_ge <- function(X){
  sqdiff.out <- optim(c(0, 0), sqdiff.pt, X = X, hessian = F)
  exp(sqdiff.out$par)}

# LS estimation

sqdiff.ls <- function(theta, X){
  lambda <- exp(theta[1])
  nu <- exp(theta[2])
  n <- length(X)
  X.sort <- sort(X)
  out <- sum((1:n / (n + 1) - pge(X.sort, lambda, nu))^2)
  out}

ls_ge <- function(X){
  sqdiff.out <- optim(c(0, 0), sqdiff.ls, X = X, hessian = F)
  exp(sqdiff.out$par)}

# WLS estimation

sqdiff.wls <- function(theta, X){
  lambda <- exp(theta[1])
  nu <- exp(theta[2])
  n <- length(X)
  X.sort <- sort(X)
  probs <- 1:n / (n + 1)
  out <- sum((probs - pge(X.sort, lambda, nu))^2 * (n + 2) / (probs * (1 - probs)))
  out}

wls_ge <- function(X){
  sqdiff.out <- optim(c(0, 0), sqdiff.wls, X = X, hessian = F)
  exp(sqdiff.out$par)}

# LM estimation

nleq.nu <- function(nu){(psi(0, 2 * nu + 1) - psi(0, nu + 1)) / (psi(0, nu + 1) - psi(0, 1))}

lme_ge <- function(X){
  n <- length(X)
  rhs.nu <- (2 * sum((1:n - 1) * sort(X)) / (n * (n - 1)) - mean(X)) / mean(X)
  nu.hat <- optimize(function(nu){(nleq.nu(nu) - rhs.nu)^2}, 
                     lower = 0, upper = 10)$minimum
  lambda.hat <- (psi(0, nu.hat + 1) - psi(0, 1)) / mean(X)
  out <- c(lambda.hat, nu.hat)
  out}

#---------------------------------
# MDPDE
#---------------------------------

# Expression for V given in Equation (11)

V <- function(x, lambda, nu, alpha){
  out <- lambda^alpha * nu^(1 + alpha) * 
    beta(1 + alpha, (nu - 1) * (1 + alpha) + 1) - 
    (1 + 1 / alpha) * dge(x, lambda, nu)^alpha
  out}

# MDPDE main estimation function

dpd <- function(theta, X, alpha){
  lambda <- exp(theta[1])
  nu <- exp(theta[2])
  mean(V(X, lambda, nu, alpha))}

dpd_ge <- function(X, alpha){
  dpd.out <- optim(c(0, 0), dpd, X = X, alpha = alpha, hessian = F)
  exp(dpd.out$par)}

# Optimal tuning parameter finder

cmvdist.mdpde.ge <- function(X, alpha){
  X <- sort(X)
  n <- length(X)
  
  mdpde.ests <- t(sapply(1:n, function(i){dpd_ge(X[-i], alpha = alpha)}))
  mean((seq_len(n) / (n + 1) - 
          sapply(1:n, function(i){pge(X[i], mdpde.ests[i, 1], mdpde.ests[i, 2])}))^2)}

optim.alpha.ge <- function(X){
  out <- optimize(function(alpha){cmvdist.mdpde.ge(X, alpha)}, lower = 0, upper = 0.5)
  out <- list(optim.alpha = out$minimum, min.cvm.dist = out$objective)
  out}

# Asymptotic covariance matrix finder

library(pracma)

xi.vec <- function(lambda, nu, alpha){
  
  xi1 <- lambda^(alpha - 1) * nu^(alpha + 1) * 
    beta(1 + alpha, (1 + alpha) * (nu - 1) + 1) * 
    (1 + psi(0, alpha + 1) - psi(0, alpha + 2))
  
  xi2 <- lambda^alpha * nu^alpha * beta(1 + alpha, (1 + alpha) * (nu - 1) + 1) * 
    (1 + nu * (psi(0, (1 + alpha) * (nu - 1) + 1) - psi(0, (1 + alpha) * nu + 1)))
  
  c(xi1, xi2)}

J.mat <- function(lambda, nu, alpha){
  
  J11 <- lambda^(alpha - 2) * nu^(alpha+1) * beta(1 + alpha, (1 + alpha)*(nu - 1) + 1) * 
    (psi(1, 1 + alpha) - psi(1, nu * (1 + alpha) + 1) + 
       (1 + psi(0, 1 + alpha) - psi(0, nu * (1 + alpha) + 1))^2 - 
       2 * (psi(1, 2 + alpha) - psi(1, nu * (1 + alpha) + 1) + 
              psi(0, 2 + alpha) - psi(0, nu * (1 + alpha) + 1) + 
              (psi(0, 2 + alpha) - psi(0, nu * (1 + alpha) + 1))^2) + 
       (nu - 1) * (alpha + 2) / ((nu - 1) * (1 + alpha) - 1) * 
       (psi(1, 3 + alpha) - psi(1, nu * (1 + alpha) + 1) + 
          (psi(0, 3 + alpha) - psi(0, nu * (1 + alpha) + 1))^2))
  
  J22 <- lambda^alpha * nu^(alpha - 1) * beta(1 + alpha, (nu - 1) * (1 + alpha) + 1) * 
    ((1 + nu * (psi(0, (nu - 1) * (1 + alpha) + 1) - psi(0, nu * (1 + alpha) + 1)))^2 + 
       nu^2 * (psi(1, (nu - 1) * (1 + alpha) + 1) - psi(1, nu * (1 + alpha) + 1)))
  
  J12 <- J21 <- lambda^(alpha - 1) * nu^alpha * beta(1 + alpha, (1 + alpha) * (nu - 1) + 1) * 
    (1 + psi(0, alpha + 1) - psi(0, alpha + 2) + 
       nu * ((psi(0, (1 + alpha) * (nu - 1) + 1) - psi(0, nu * (1 + alpha) + 1)) * 
               (1 + psi(0, 1 + alpha) - psi(0, nu * (1 + alpha) + 1)) - 
               (psi(0, alpha + 2) - psi(0, nu * (1 + alpha) + 1)) * 
               (psi(0, (1 + alpha) * (nu - 1)) - psi(0, nu * (1 + alpha) + 1))))
  
  J <- matrix(NA, 2, 2)
  J[1, 1] <- J11
  J[2, 2] <- J22
  J[1, 2] <- J[2, 1] <- J12
  J}

covmat.theta.alpha <- function(lambda, nu, alpha){
  
  J <- J.mat(lambda, nu, alpha)
  xi <- xi.vec(lambda, nu, alpha)
  K <- J.mat(lambda, nu, 2 * alpha) - tcrossprod(xi)
  
  (solve(J) %*% (K %*% solve(J)))}

# Asymptotic relative efficiency

# nu.alpha <- vector of nu and alpha, independent of lambda

are.fn <- function(nu.alpha){
  out1 <- covmat.theta.alpha(1, nu.alpha[1], nu.alpha[2])
  out1 <- c(out1[1, 1], out1[2, 2])
  out2 <- covmat.theta.alpha(1, nu.alpha[1], 0)
  out2 <- c(out2[1, 1], out2[2, 2])
  out2 / out1}
