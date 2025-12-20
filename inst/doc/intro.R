## ----eval=FALSE---------------------------------------------------------------
# G.hat <- function(t, Dat) {
# newDat = list()
# newDat$x = Dat$x
# newDat$Delta = 1 - Dat$Delta
# newDat$Y = Dat$Y
# fit = survfit(Surv(Y, Delta) ~ 1, data = newDat)
# if (t <= fit$time[1]) {
# return(1)
# } else {
# q = max(which(fit$time <= t))
# return(fit$surv[q])
# }
# }

## ----eval=FALSE---------------------------------------------------------------
# library(survival)
# set.seed(123)
# n <- 100
# Dat <- list(
# x = matrix(rnorm(n * 3), n, 3),
# Y = rexp(n),
# Delta = rbinom(n, 1, 0.7)
# )
# G.hat(0.5, Dat)

## ----eval=FALSE---------------------------------------------------------------
# S.hat <- function(newX, newX_test, ns, ns_test, t, Datas, Datas_test, h) {
#   p <- ncol(newX)
#   K <- ceiling(p / h)
# 
#   S <- NULL
#   S_test <- NULL
# 
#   Dat1 <- Datas
#   Dat1$X <- newX  ## reordered/processed X
# 
#   Dat_test1 <- Datas_test
#   Dat_test1$X <- newX_test
# 
#   pre.dat <- Dat1
#   pre.dat$Y <- rep(t, ns)        ## fix time t
# 
#   pre.dat_test <- Dat_test1
#   pre.dat_test$Y <- rep(t, ns_test)
# 
#   for (k in 1:K) {
#     ind <- unique(c(1, ((k - 1) * h + 1):min(p, k * h)))
#     res <- coxph(Surv(Y, as.numeric(Delta)) ~ X[, ind], data = Dat1)
#     S.res <- stats::predict(res, newdata = pre.dat, type = "survival")
#     S_test.res <- stats::predict(res, newdata = pre.dat_test, type = "survival")
#     S <- cbind(S, S.res)
#     S_test <- cbind(S_test, S_test.res)
#   }
# 
#   list(S = S, S_test = S_test)
# }

## ----eval=FALSE---------------------------------------------------------------
# h <- 2
# newX <- Dat$x
# newX_test <- matrix(rnorm(50 * 3), 50, 3)
# 
# Dat_test <- list(
# x = newX_test,
# Y = rexp(50),
# Delta = rbinom(50, 1, 0.7)
# )
# 
# res <- S.hat(
# newX = newX,
# newX_test = newX_test,
# ns = n,
# ns_test = 50,
# t = 0.5,
# Datas = Dat,
# Datas_test = Dat_test
# )

## ----eval=FALSE---------------------------------------------------------------
# bs <- function(newX, newX_test, ns, ns_test, t, Dat, Dat_test, k) {
#   bs0 <- c()
#   res.S <- S.hat(newX, newX_test, ns, ns_test, t, Dat, Dat_test)
# 
#   for (i in 1:ns) {
#     if (Dat$Y[i] > t) {
#       bs0 <- c(bs0, (1 - res.S$S[i, k])^2 / G.hat(t, Dat))
#     } else if (as.numeric(Dat$Delta)[i] == 1) {
#       bs0 <- c(bs0, (0 - res.S$S[i, k])^2 / G.hat(Dat$Y[i], Dat))
#     }
#   }
# 
#   sum(bs0) / ns
# }

## ----eval=FALSE---------------------------------------------------------------
# bs(
# newX = newX,
# newX_test = newX_test,
# ns = n,
# ns_test = 50,
# t = 0.5,
# Dat = Dat,
# Dat_test = Dat_test,
# k = 1
# )

## ----eval=FALSE---------------------------------------------------------------
# Rcpp::cppFunction(
# List generate_subsetsC(IntegerVector set) {
# int n = set.size();
# int total = (1 << n);
# List subsets;
# 
# for (int mask = 1; mask < total; ++mask) {
# IntegerVector subset;
# for (int i = 0; i < n; ++i) {
# if (mask & (1 << i)) {
# subset.push_back(set[i]);
# }
# }
# subsets.push_back(subset);
# }
# return subsets;
# }
# )

## ----eval=FALSE---------------------------------------------------------------
# library(SA25204001)
# generate_subsetsC(1:3)

## ----eval=FALSE---------------------------------------------------------------
# Rcpp::cppFunction(
# List generate_nested_subsetsC(IntegerVector set) {
# int n = set.size();
# List subsets(n);
# for (int i = 0; i < n; ++i) {
# subsets[i] = set[Range(0, i)];
# }
# return subsets;
# }
# )

## ----eval=FALSE---------------------------------------------------------------
# library(SA24204150)
# generate_nested_subsetsC(1:5)

