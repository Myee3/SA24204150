#' @useDynLib SA24204150, .registration = TRUE
#' @importFrom Rcpp evalCpp
NULL


#' @title Kaplan–Meier estimator of the survival function with censoring
#' @description
#' Compute \eqn{\hat G(t)} using a Kaplan–Meier estimator for the censoring
#' distribution, which is used for inverse probability of censoring weighting (IPCW).
#'
#'
#' @param t A time point at which \eqn{\hat G(t)} is evaluated (numeric scalar).
#' @param Dat A data object containing at least \code{$x}, \code{$Delta}, and \code{$Y}.
#' \describe{
#'   \item{x}{Covariates (not used in the computation but kept for consistency).}
#'   \item{Delta}{Event indicator (1 = event, 0 = censored).}
#'   \item{Y}{Observed time.}
#' }
#'
#' @return A numeric scalar giving \eqn{\hat G(t)}.
#'
#' @details
#' If \code{t} is smaller than the first observed time in the fitted Kaplan–Meier
#' curve, the function returns 1.
#'
#' @examples
#' \dontrun{
#' library(survival)
#' Dat <- list(
#'   x = matrix(rnorm(100), ncol = 2),
#'   Delta = rbinom(50, 1, 0.7),
#'   Y = rexp(50)
#' )
#' G.hat(0.5, Dat)
#' }
#'
#' @importFrom survival survfit Surv
#' @export
G.hat <- function(t, Dat) {
  newDat <- list()
  newDat$x <- Dat$x
  newDat$Delta <- 1 - Dat$Delta  ## IPCW: treat censoring as event
  newDat$Y <- Dat$Y
  fit <- survfit(Surv(Y, Delta) ~ 1, data = newDat)
  
  if (t <= fit$time[1]) {
    return(1)
  } else {
    q <- max(which(fit$time <= t))
    return(fit$surv[q])
  }
}


#' @title Cox-model-based survival prediction across grouped covariates
#' @description
#' Fit a sequence of Cox proportional hazards models on grouped covariates and
#' return predicted survival probabilities at a fixed time \code{t} for both
#' training and test samples.
#'
#' The covariates are split into \code{K = ceiling(p / h)} groups, where \code{p}
#' is the number of columns in \code{newX} and \code{h} is a global group size
#' (must exist in the calling environment).
#'
#' @param newX A numeric matrix of training covariates (dimension \code{ns x p}).
#' @param newX_test A numeric matrix of test covariates (dimension \code{ns_test x p}).
#' @param ns Number of training samples (should equal \code{nrow(newX)}).
#' @param ns_test Number of test samples (should equal \code{nrow(newX_test)}).
#' @param t The prediction time point (numeric scalar). Survival probabilities are
#' predicted at time \code{t}.
#' @param Datas A data object for training containing at least \code{$Y}, \code{$Delta},
#' and an \code{$X} field will be created internally.
#' @param Datas_test A data object for testing containing at least \code{$Y}, \code{$Delta},
#' and an \code{$X} field will be created internally.
#' @param h Number of covariates in each candidate model.
#'
#' @return A list with two components:
#' \describe{
#'   \item{S}{A numeric matrix of predicted survival probabilities at time \code{t}
#'   for training data (dimension \code{ns x K}).}
#'   \item{S_test}{A numeric matrix of predicted survival probabilities at time \code{t}
#'   for test data (dimension \code{ns_test x K}).}
#' }
#'
#' @details
#' The index set for the \code{k}-th Cox model is
#' \code{ind = unique(c(1, ((k-1)*h+1):min(p, k*h)))}.
#' This always includes the first covariate.
#'
#' @examples
#' \dontrun{
#' library(survival)
#' set.seed(1)
#' h <- 5
#' ns <- 100; ns_test <- 50; p <- 12
#' newX <- matrix(rnorm(ns * p), ns, p)
#' newX_test <- matrix(rnorm(ns_test * p), ns_test, p)
#' Datas <- data.frame(Y = rexp(ns), Delta = rbinom(ns, 1, 0.7))
#' Datas_test <- data.frame(Y = rexp(ns_test), Delta = rbinom(ns_test, 1, 0.7))
#'
#' out <- S.hat(newX, newX_test, ns, ns_test, t = 0.5, Datas, Datas_test)
#' dim(out$S)
#' dim(out$S_test)
#' }
#'
#' @importFrom survival coxph Surv
#' @export
S.hat <- function(newX, newX_test, ns, ns_test, t, Datas, Datas_test, h) {
  p <- ncol(newX)
  K <- ceiling(p / h)
  
  S <- NULL
  S_test <- NULL
  
  Dat1 <- Datas
  Dat1$X <- newX  ## reordered/processed X
  
  Dat_test1 <- Datas_test
  Dat_test1$X <- newX_test
  
  pre.dat <- Dat1
  pre.dat$Y <- rep(t, ns)        ## fix time t
  
  pre.dat_test <- Dat_test1
  pre.dat_test$Y <- rep(t, ns_test)
  
  for (k in 1:K) {
    ind <- unique(c(1, ((k - 1) * h + 1):min(p, k * h)))
    res <- coxph(Surv(Y, as.numeric(Delta)) ~ X[, ind], data = Dat1)
    S.res <- stats::predict(res, newdata = pre.dat, type = "survival")
    S_test.res <- stats::predict(res, newdata = pre.dat_test, type = "survival")
    S <- cbind(S, S.res)
    S_test <- cbind(S_test, S_test.res)
  }
  
  list(S = S, S_test = S_test)
}


#' @title IPCW Brier score at a fixed time for a selected candidate Cox model
#' @description
#' Compute the inverse-probability-of-censoring-weighted (IPCW) Brier score at time
#' \code{t} for the \code{k}-th candidate Cox model produced by \code{\link{S.hat}}.
#'
#' @param newX A numeric matrix of training covariates (dimension \code{ns x p}).
#' @param newX_test A numeric matrix of test covariates (dimension \code{ns_test x p}).
#' @param ns Number of training samples.
#' @param ns_test Number of test samples.
#' @param t The time point at which the Brier score is evaluated (numeric scalar).
#' @param Dat Training data containing at least \code{$Y}, \code{$Delta}, and \code{$X}
#' (or will be created inside \code{\link{S.hat}}).
#' @param Dat_test Test data containing at least \code{$Y}, \code{$Delta}, and \code{$X}
#' (or will be created inside \code{\link{S.hat}}).
#' @param k An integer specifying which candidate model (column) to use from the
#' output of \code{\link{S.hat}}.
#'
#' @return A numeric scalar giving the IPCW Brier score at time \code{t}.
#'
#' @details
#' For each subject \code{i}:
#' \itemize{
#' \item If \code{Y_i > t}, the contribution is
#' \eqn{(1 - \hat S_k(t | X_i))^2 / \hat G(t)}.
#' \item If \code{Y_i <= t} and \code{Delta_i = 1} (event observed),
#' the contribution is
#' \eqn{(0 - \hat S_k(t | X_i))^2 / \hat G(Y_i)}.
#' \item Otherwise (censored before \code{t}), the observation is skipped.
#' }
#'
#' @examples
#' \dontrun{
#' library(survival)
#' set.seed(1)
#' h <- 5
#' ns <- 80; ns_test <- 40; p <- 10
#' newX <- matrix(rnorm(ns * p), ns, p)
#' newX_test <- matrix(rnorm(ns_test * p), ns_test, p)
#' Dat <- data.frame(Y = rexp(ns), Delta = rbinom(ns, 1, 0.7))
#' Dat_test <- data.frame(Y = rexp(ns_test), Delta = rbinom(ns_test, 1, 0.7))
#'
#' bs(newX, newX_test, ns, ns_test, t = 0.8, Dat, Dat_test, k = 1)
#' }
#'
#' @seealso \code{\link{G.hat}}, \code{\link{S.hat}}
#' @export
bs <- function(newX, newX_test, ns, ns_test, t, Dat, Dat_test, k) {
  bs0 <- c()
  res.S <- S.hat(newX, newX_test, ns, ns_test, t, Dat, Dat_test)
  
  for (i in 1:ns) {
    if (Dat$Y[i] > t) {
      bs0 <- c(bs0, (1 - res.S$S[i, k])^2 / G.hat(t, Dat))
    } else if (as.numeric(Dat$Delta)[i] == 1) {
      bs0 <- c(bs0, (0 - res.S$S[i, k])^2 / G.hat(Dat$Y[i], Dat))
    }
  }
  
  sum(bs0) / ns
}

