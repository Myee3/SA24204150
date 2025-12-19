## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)
library(ggplot2)

## ----warning=FALSE------------------------------------------------------------
set.seed(1111111)
n = 10000 #进行10000次模拟
t <- rep(0,10000)
for(i in 1:n){
  draw <- sample(54) #进行一次抽取
  for (j in 1:54) {
    if(draw[j] == j){
      t[i] = t[i] + 1 #进行一次成功抽取
    }
    else{t[i]=t[i]}
  }
}
mean(t)
var(t)

hist(t,
     probability = TRUE,
     main = "Histogram of Replications",
     xlab = "Value")

## -----------------------------------------------------------------------------
# Beta(n,1)分布
set.seed(111)
beta <- function(n,N){
  u <- runif(N,0,1)
  x <- u^(1/n)
  x
}
t <- beta(5,100) # n = 5
knitr::kable(head(t, 10), caption = "随机产生的数据")
plot.ecdf(t, main="Beta Distribution")
y <- seq(min(t), max(t), 0.01)
lines(y, pbeta(y, shape1 = 5,shape2 = 1), 
      col="red", type="h", lwd=2)

## -----------------------------------------------------------------------------
set.seed(111)
cauchy <- function(N){
  u <- runif(N,0,1)
  x <- tan(pi*(u-0.5))
  x
}
t <- cauchy(100)
knitr::kable(head(t, 10), caption = "随机产生的数据")

plot.ecdf(t, main="Cauchy Distribution")
y <- seq(min(t), max(t), 10)
lines(y, pcauchy(y), col="red", type="h", lwd=2)

## ----fig.width=10, fig.height=8-----------------------------------------------
set.seed(123)
n = 1000
generate_rl <- function(sigma) {
  u = runif(n)
  x = sigma*sqrt(-2*log(1-u))
  return(x)
}

par(mfrow = c(2, 2), mar = c(3.5, 3.5, 2, 1), mgp = c(2.1, 0.7, 0))
sigma = c(0.5, 1, 2, 10)   # 考虑不同取值
for (i in sigma) {
  x = generate_rl(i)
  hist(x, breaks = 50, freq = FALSE, 
       main = paste("Rayleigh Distribution (σ =", i, ")"),
       xlab = "x", ylab = "Density")
  x_true = seq(0, max(x), length.out = 1000)
  y_true = (x_true/i^2)*exp(-x_true^2/(2*i^2))
  lines(x_true, y_true, col = "red", lwd = 2)
  legend("topright", c("True"), col = c("red"), lwd = 2)
}

## -----------------------------------------------------------------------------
set.seed(123)
n = 1000
x = c(0:4)
p = c(0.1, 0.2, 0.2, 0.2, 0.3)
cum_p <- cumsum(p)

# 逆变换法
generate_x <- function(n, x, cum_p){
  u = runif(n)
  samples = numeric(n)
  for (i in 1:n) {
    samples[i] <- x[which(cum_p >= u[i])[1]]
  }
  return(samples)
}

samples_it = generate_x(n, x, cum_p)
ftable = table(samples_it) / n
comparison <- data.frame(
  x = x,
  Theoretical = p,
  Empirical = as.numeric(ftable[as.character(x)])
)
comparison

## -----------------------------------------------------------------------------
set.seed(123)
samples_r = sample(x, n, replace = TRUE, prob = p)
ftable_r = table(samples_r) / n
comparison_r <- data.frame(
  x = x,
  Theoretical = p,
  Empirical = as.numeric(ftable_r[as.character(x)])
)
comparison_r

## -----------------------------------------------------------------------------
set.seed(123)
n = 1000; i = 0; k = 0
generate_beta <- function(n, a, b){
  if (a>1 && b>1) {
    mode = (a-1)/(a+b-2)
    M = dbeta(mode, a, b)
  }else{
    M = max(dbeta(seq(0.01, 0.99, 0.01), a, b)) + 0.1
  }
  samples = numeric(n)
  while (i < n) {
    x = runif(1)
    k = k + 1
    u = runif(1)
    if (u <= dbeta(x, a, b)/M) {
      i = i + 1
      samples[i] = x
    }
  }
  return(samples)
}

samples = generate_beta(1000, 3, 2)

hist(samples, breaks = 20, freq = FALSE, 
     main = "Beta(3,2) Distribution",
     xlab = "x", ylab = "Density")
x_true = seq(0, 1, length.out = 1000)
y_true = dbeta(x_true, 3, 2)
lines(x_true, y_true, col = "red", lwd = 2)
legend("topleft", "Theoretical Distribution", col = "red", lwd = 2)

## -----------------------------------------------------------------------------
set.seed(123)

draw_hist <- function(p1){
  n = 1000
  X1 = rnorm(n, 0, 1) 
  X2 = rnorm(n, 3, 1)
  r = sample(c(0, 1), n, replace = TRUE, prob = c(1-p1, p1))
  Z = r * X1 + (1-r) * X2
  
  hist(Z, breaks = 20, freq = FALSE,
       main = "Normal Mixture",
       xlab = "x", ylab = "Density")
  x_true = seq(min(Z), max(Z), length.out = 1000)
  y_true = p1 * dnorm(x_true, 0, 1) + (1 - p1) * dnorm(x_true, 3, 1)
  lines(x_true, y_true, col = "red", lwd = 2)
}

draw_hist(0.75)

## -----------------------------------------------------------------------------
draw_hist(0.5)
draw_hist(0.25)

## -----------------------------------------------------------------------------
set.seed(123)
n = 1000
r = 4
beta = 2
lambda = rgamma(n, shape = r, rate = beta)
y = rexp(n, lambda)
hist(y, freq = TRUE, breaks = 100,
     main = "Exponential-Gamma mixture",
     xlab = "y", ylab = "Freq")

## -----------------------------------------------------------------------------
set.seed(123)
MC_Beta <- function(x, n = 1e5) {
  samples = rbeta(n, 3, 3)
  est = sapply(x, function(xi) mean(samples <= xi))
  return(est)
}

x = seq(0.1, 0.9, by = 0.1)
est_MC = MC_Beta(x)

true = pbeta(x, 3, 3)

## -----------------------------------------------------------------------------
compare = data.frame(
  x = x,
  MC_Estimate = round(est_MC, 4),
  True_pbeta = round(true, 4)
)
compare

## -----------------------------------------------------------------------------
n = 10000
set.seed(123)
U = runif(n)
Y1 = exp(U)
Y2 = exp(1 - U)

## -----------------------------------------------------------------------------
var_simple = var(Y1)
var_simple

## -----------------------------------------------------------------------------
cov_Y1Y2 = cov(Y1, Y2)
cov_Y1Y2

## -----------------------------------------------------------------------------
var_Y = var(Y1 + Y2)
var_Y

## -----------------------------------------------------------------------------
var_anti = var((Y1+Y2)/2)
reduction = 1 - var_anti / var_simple
cat("Variance reduction =", reduction, ".\n")

## -----------------------------------------------------------------------------
n = 10000
set.seed(123)

# f1
x1 = rexp(n, 1) + 1
w1 = (x1^2/sqrt(2*pi)*exp(-x1^2/2))/exp(-(x1-1))
I1 = mean(w1)
var1 = var(w1)

# f2
rtnorm <- function(n) {
  u = runif(n, pnorm(q=1, mean=1, sd=1), 1)
  qnorm(u, mean=1, sd=1)
}

x2 = rtnorm(n)
dens_f2 = dnorm(x2, 1, 1) / (1 - pnorm(1, 1, 1))
w2 = (x2^2 / sqrt(2*pi) * exp(-x2^2/2)) / dens_f2
I2 = mean(w2)
var2 = var(w2)
var2/var1

## -----------------------------------------------------------------------------
library(rbenchmark)
set.seed(123)
n = c(1e4, 2e4, 4e4, 6e4, 8e4)
an = numeric(5)
for (i in 1:5) {
  ni = n[i]
  x = sample(1:ni, ni,replace = FALSE) 
  bm = benchmark(sort(x), replications = 1000, columns = c("test", "replications", "elapsed"))
  an[i] = bm$elapsed 
}
an

## -----------------------------------------------------------------------------
tn = n * log(n)
fit = lm(an ~ tn)
summary(fit)

plot(x=tn, y=an, pch = 19, col = "blue",
     xlab = expression(t[n]),
     ylab = expression(a[n]),
     main = "an vs tn")
abline(fit, col = "red", lwd = 2)
legend("topleft", legend = c("Computation Time", "Regression line"), 
       col = c("blue", "red"), pch = c(19, NA), lty = c(NA, 1), lwd = 2)

## -----------------------------------------------------------------------------
rm(list = ls())

rf3 <- function(n) { u <- runif(n); -log(1 - u*(1 - exp(-1))) }

gen_data <- function(m, K, seed){
  set.seed(seed)
  m_j = rep(m/K, K)
  a = (0:(K-1))/K
  b = (1:K)/K
  X = vector("list", K)
  for(j in 1:K){
    u = runif(m_j[j])
    p = a[j]
    q = b[j]
    X[[j]] = -log(exp(-p)-u*(exp(-p)-exp(-q)))
  }
  list(X = X, a = a, b = b, m_j = m_j)
}
m = 1e4; K = 5; seed = 123
res = gen_data(m, K, seed)
save(res, file = "E:/code/gen_data.RData")

## -----------------------------------------------------------------------------
rm(list = ls())
load("E:/code/gen_data.RData")

stat_infer <- function(data){
  a = data$a; b = data$b; X = data$X; m_j = data$m_j; K = length(X)
  Pj = (exp(-a) - exp(-b))/(1 - exp(-1))         
  r_bar = s2 = numeric(K)
  for(j in 1:K){
    rj = (1 - exp(-1))/(1 + X[[j]]^2)
    r_bar[j] = mean(rj)
    s2[j] = var(rj)
  }
  theta_hat = sum(Pj * r_bar)
  var_hat = sum(Pj^2 * s2 / m_j)
  se_hat = sqrt(var_hat)
  cbind(theta_hat = theta_hat, se_hat = se_hat)
}

inf = as.data.frame(stat_infer(data = res))
write.csv(inf, file = "E:/code/newinf.csv")

## -----------------------------------------------------------------------------
rm(list = ls())
inf = read.csv("E:/code/newinf.csv")

print_results <- function(inf){
  out <- data.frame(
    method = c("Eg6.14", "Eg6.11"),
    theta_hat = c(inf$theta_hat, 0.5257801),
    mc_se     = c(inf$se_hat, 0.0970314)
  )
  out
}
print_results(inf) 

## -----------------------------------------------------------------------------
rm(list = ls())

gen_data <- function( n_seq, mu_seq, sigma, m, file){
  set.seed(123)
  dat <- list()
  for(n in n_seq){
    for(mu in mu_seq){
      dat[[paste0("n",n,"_mu",mu)]] <-
        replicate(m, rnorm(n, mu, sigma), simplify = FALSE)
    }
  }
  dat
}

res = gen_data(
  n_seq = c(10,20,30,40,50),
  mu_seq = seq(450,650,10),
  sigma = 100, m = 1000,
  file = "gen_data.rds"
)

save(res, file = "E:/code/gen_data1.RData")

## -----------------------------------------------------------------------------
rm(list = ls())
load("E:/code/gen_data1.RData")
calc_power <- function(res, mu0, alpha){
  dat <- res
  inf <- data.frame(n = numeric(), mu = numeric(), power = numeric())

  for (nm in names(dat)) {
    parts <- as.numeric(unlist(regmatches(nm, gregexpr("\\d+", nm))))
    n <- parts[1]; mu <- parts[2]
    pvals <- sapply(dat[[nm]], function(x)
      t.test(x, alternative = "greater", mu = mu0)$p.value)
    inf <- rbind(inf, data.frame(n = n, mu = mu, power = mean(pvals <= alpha)))
  }

  inf
}
inf = calc_power(res, 500, 0.05)
write.csv(inf, "E:/code/inf1.csv")

## -----------------------------------------------------------------------------
rm(list = ls())
inf = read.csv("E:/code/inf1.csv")
plot_power <- function(inf){
  cols = 1:length(unique(inf$n))
  plot(NA, xlim = range(inf$mu), ylim = c(0,1),
       xlab = expression(mu), ylab = "Empirical Power",
       main = "Empirical Power Curves")
  abline(v = 500, lty = 2, col = "gray40")

  for(i in seq_along(unique(inf$n))){
    n = unique(inf$n)[i]
    lines(inf$mu[inf$n==n], inf$power[inf$n==n],
          col = cols[i], lwd = 2)
  }
  legend("bottomright", legend = paste0("n=", unique(inf$n)),
         col = cols, lwd = 2, cex = .9)
}
png("E:/code/power_plot.png")
plot_power(inf)
dev.off()

plot_power(inf)

## -----------------------------------------------------------------------------
rm(list = ls())
gen_data <- function(n, m){
  set.seed(123)
  dat = replicate(m, rchisq(n, df = 2), simplify = FALSE)
  dat
}
dat = gen_data(n=20, m=1000)
save(dat, file = "E:/code/gen_data2.RData")

# ------------------------------------------

rm(list = ls())
load("E:/code/gen_data2.RData")
infer_coverage <- function(dat, mu_true, alpha){
  n = length(dat[[1]])
  cover = sapply(dat, function(x) {
    tcrit = qt(1-alpha/2, df = n-1)
    ci_low = mean(x)-tcrit*sd(x)/sqrt(n)
    ci_up  = mean(x)+tcrit*sd(x)/sqrt(n)
    mu_true >= ci_low & mu_true <= ci_up
  })
  res = data.frame(
    n = n, alpha = alpha,
    coverage = mean(cover),
    se = sqrt(mean(cover)*(1-mean(cover))/length(cover))
  )
  res
}

inf = infer_coverage(dat, 2, 0.05)
write.csv(inf, "E:/code/inf2.csv")

# -------------------------------------------

rm(list = ls())
inf = read.csv("E:/code/inf2.csv")
report_result <- function(inf){
  cat("\n Monte Carlo 估计结果 \n")
  print(inf)
}
report_result(inf)


## -----------------------------------------------------------------------------
rm(list = ls())
gen_data <- function(n, m){
  set.seed(123)
  dat = replicate(m, rchisq(n, df = 2), simplify = FALSE)
  dat
}
dat = gen_data(n=20, m=1000)
save(dat, file = "E:/code/gen_data2.RData")

# ------------------------------------------

rm(list = ls())
load("E:/code/gen_data2.RData")
infer_coverage <- function(dat, mu_true, alpha){
  n = length(dat[[1]])
  cover = sapply(dat, function(x) {
    crit = qnorm(1-alpha/2)
    ci_low = mean(x)-crit*sd(x)/sqrt(n)
    ci_up  = mean(x)+crit*sd(x)/sqrt(n)
    mu_true >= ci_low & mu_true <= ci_up
  })
  res = data.frame(
    n = n, alpha = alpha,
    coverage = mean(cover),
    se = sqrt(mean(cover)*(1-mean(cover))/length(cover))
  )
  res
}

inf = infer_coverage(dat, 2, 0.05)
write.csv(inf, "E:/code/inf2.csv")

# -------------------------------------------

rm(list = ls())
inf = read.csv("E:/code/inf2.csv")
report_result <- function(inf){
  cat("\n Monte Carlo 估计结果 \n")
  print(inf)
}
report_result(inf)


## -----------------------------------------------------------------------------
rm(list = ls())
gen_data <- function(n, m){
  set.seed(123)
  dat = replicate(m, runif(n, min = 0, max = 2), simplify = FALSE)
  dat
}
dat = gen_data(n=20, m=2000)
save(dat, file = "E:/code/gen_data4.RData")

# ------------------------------------------

rm(list = ls())
load("E:/code/gen_data4.RData")
infer_t1e <- function(dat, mu_true, alpha){
  n = length(dat[[1]])
  pvals = sapply(dat, function(x)
    t.test(x, alternative = "two.sided", mu = mu_true)$p.value)
  type1 = mean(pvals <= alpha)  
  res = data.frame(n = n, alpha = alpha, type1 = type1)
}

inf = infer_t1e(dat, 1, 0.05)
write.csv(inf, "E:/code/inf4.csv")

# -------------------------------------------

rm(list = ls())
inf = read.csv("E:/code/inf4.csv")
report_result <- function(inf){
  cat("\n Monte Carlo 估计结果 \n")
  print(inf)
}
report_result(inf)

## -----------------------------------------------------------------------------
rm(list = ls())
gen_data <- function(n, m){
  set.seed(123)
  dat = replicate(m, rchisq(n, df = 1), simplify = FALSE)
  dat
}
dat = gen_data(n=20, m=2000)
save(dat, file = "E:/code/gen_data3.RData")

# ------------------------------------------

rm(list = ls())
load("E:/code/gen_data3.RData")
infer_t1e <- function(dat, mu_true, alpha){
  n = length(dat[[1]])
  pvals = sapply(dat, function(x)
    t.test(x, alternative = "two.sided", mu = mu_true)$p.value)
  type1 = mean(pvals <= alpha)  
  res = data.frame(n = n, alpha = alpha, type1 = type1)
}

inf = infer_t1e(dat, 1, 0.05)
write.csv(inf, "E:/code/inf3.csv")

# -------------------------------------------

rm(list = ls())
inf = read.csv("E:/code/inf3.csv")
report_result <- function(inf){
  cat("\n Monte Carlo 估计结果 \n")
  print(inf)
}
report_result(inf)

## -----------------------------------------------------------------------------
rm(list = ls())
gen_data <- function(n, m){
  set.seed(123)
  dat = replicate(m, rexp(n, rate = 1), simplify = FALSE)
  dat
}
dat = gen_data(n=20, m=2000)
save(dat, file = "E:/code/gen_data5.RData")

# ------------------------------------------

rm(list = ls())
load("E:/code/gen_data5.RData")
infer_t1e <- function(dat, mu_true, alpha){
  n = length(dat[[1]])
  pvals = sapply(dat, function(x)
    t.test(x, alternative = "two.sided", mu = mu_true)$p.value)
  type1 = mean(pvals <= alpha)  
  res = data.frame(n = n, alpha = alpha, type1 = type1)
}

inf = infer_t1e(dat, 1, 0.05)
write.csv(inf, "E:/code/inf5.csv")

# -------------------------------------------

rm(list = ls())
inf = read.csv("E:/code/inf5.csv")
report_result <- function(inf){
  cat("\n Monte Carlo 估计结果 \n")
  print(inf)
}
report_result(inf)

## -----------------------------------------------------------------------------
infer_multi <- function(N = 1000, m0 = 950, alpha = 0.1, m = 10000, seed = 123, out_file = "multi_step1.rds"){
  set.seed(seed)
  fwer_b = fwer_bh = fdr_b = fdr_bh = tpr_b = tpr_bh = 0
  m1 = N - m0
  for(rep in 1:m){
    p = c(runif(m0), rbeta(m1, 0.1, 1)) 
    padj_b  = p.adjust(p, "bonferroni")
    padj_bh = p.adjust(p, "BH")
    for (padj in list(padj_b, padj_bh)) {
      rej <- padj <= alpha
      R = sum(rej); V = sum(rej[1:m0]); S = sum(rej[(m0+1):N])
      fwer = as.integer(V > 0)
      fdr = if (R > 0) V / R else 0
      tpr = S / m1
      if (identical(padj, padj_b)) { fwer_b = fwer_b + fwer; fdr_b = fdr_b + fdr; tpr_b = tpr_b + tpr }
      else                         { fwer_bh= fwer_bh+ fwer; fdr_bh= fdr_bh+ fdr; tpr_bh= tpr_bh+ tpr }
      }
    }
    res <- matrix(c(fwer_b, fwer_bh,  # FWER
                    fdr_b,  fdr_bh,   # FDR
                    tpr_b,  tpr_bh),  # TPR
                  nrow=3, byrow=TRUE) / m
    dimnames(res) <- list(c("FWER","FDR","TPR"),
                          c("Bonferroni","B-H"))
    saveRDS(list(table=res, params=par), out_file)
    invisible(res)
}

result_multi <- function(in_file = "multi_step1.rds"){
  obj = readRDS(in_file)
  tab = obj$table
  par = obj$params
  print(round(tab, 4))
}

infer_multi()
result_multi()

## -----------------------------------------------------------------------------
gen_aircond <- function(file = "aircond_step1.rds") {
  x = c(3,5,7,18,43,85,91,98,100,130,230,487)
  saveRDS(x, file); invisible(x)
}

infer_aircond <- function(in_file = "aircond_step1.rds",
                          B, seed,
                          out_file = "aircond_step2.rds") {
  x = readRDS(in_file)
  n = length(x)
  lam_hat = 1 / mean(x)
  set.seed(seed)
  lam_star = replicate(B, {1/mean(sample(x, n, replace = TRUE))})
  bias = mean(lam_star) - lam_hat  
  se = sd(lam_star)    
  res = list(lam_hat = lam_hat, bias = bias, se = se)
  saveRDS(res, out_file); invisible(res)
}

result_aircond <- function(in_file = "aircond_step2.rds") {
  r = readRDS(in_file)
  cat(sprintf("MLE Estimetor = %.6f\n", r$lam_hat))
  cat(sprintf("Bias (bootstrap) = %.6f\n", r$bias))
  cat(sprintf("SE   (bootstrap) = %.6f\n", r$se))
}

gen_aircond()
infer_aircond(B = 10000, seed = 1)
result_aircond("aircond_step2.rds")

## -----------------------------------------------------------------------------
library(bootstrap)
gen_scor <- function(file = "scor_step1.rds"){
  X = as.matrix(scor)
  saveRDS(X, file); invisible(X)
}


infer_theta <- function(in_file = "scor_step1.rds",
                        B, seed,
                        out_file = "scor_step2.rds"){
  X = readRDS(in_file)
  n = nrow(X)
  get_theta <- function(M){
    ev = eigen(cov(M), only.values = TRUE)$values
    max(ev) / sum(ev)
  }
  theta_hat = get_theta(X)

  set.seed(seed)
  thetas = replicate(B, {
    idx = sample.int(n, n, replace = TRUE)
    get_theta(X[idx, , drop = FALSE])
  })
  res = list(
    theta_hat = theta_hat,
    bias = mean(thetas) - theta_hat,
    se   = sd(thetas)
  )
  saveRDS(res, out_file); invisible(res)
}

result_theta <- function(in_file = "scor_step2.rds"){
  r = readRDS(in_file)
  cat(sprintf("Theta hat = %.4f\n", r$theta_hat))
  cat(sprintf("Bootstrap bias = %.4f\n", r$bias))
  cat(sprintf("Bootstrap SE = %.4f\n", r$se))
}

gen_scor()
infer_theta(B = 10000, seed = 123)
result_theta("scor_step2.rds")

## -----------------------------------------------------------------------------
library(bootstrap)
gen_scor <- function(file = "scor_step1.rds") {
  X = as.matrix(scor)
  saveRDS(X, file); invisible(X)
}

infer_jackknife <- function(in_file = "scor_step1.rds",
                            out_file = "scor_jack_step2.rds") {
  X = readRDS(in_file)
  n = nrow(X)
  get_theta <- function(M) {
    ev = eigen(cov(M), only.values = TRUE)$values
    max(ev)/sum(ev)
  }
  theta_hat = get_theta(X)         # 原始估计
  theta_jack = numeric(n)          # leave-one-out
  for (i in 1:n){
  theta_jack[i] = get_theta(X[-i, , drop = FALSE])
  }
  theta_bar = mean(theta_jack)
  bias_jack = (n-1)*(theta_bar-theta_hat)
  se_jack = sqrt((n-1)/n*sum((theta_jack-theta_bar)^2))

  res = list(bias_jack = bias_jack, se_jack = se_jack)
  saveRDS(res, out_file); invisible(res)
}

report_jackknife <- function(in_file = "scor_jack_step2.rds") {
  r = readRDS(in_file)
  cat(sprintf("Jackknife 偏差估计 = %.6f\n", r$bias_jack))
  cat(sprintf("Jackknife 标准误 = %.6f\n", r$se_jack))
}

gen_scor()
infer_jackknife()
report_jackknife()

## ----warning=FALSE------------------------------------------------------------
library(DAAG)
data("ironslag")

magnetic <- ironslag$magnetic
chemical <- ironslag$chemical
n <- length(magnetic)
N2 <- choose(n, 2)
e1 <- e2 <- e3 <- e4 <- numeric(N2)

pairs <- combn(n, 2)
for (k in 1:N2) {
  idx  <- pairs[, k]  
  y.tr <- magnetic[-idx]
  x.tr <- chemical[-idx]

  y.te <- magnetic[idx]
  x.te <- chemical[idx]

  ## Model 1: y ~ x
  J1 <- lm(y.tr ~ x.tr)
  yhat1 <- J1$coef[1] + J1$coef[2] * x.te
  e1[k] <- mean((y.te - yhat1)^2)

  ## Model 2: y ~ x + x^2
  J2 <- lm(y.tr ~ x.tr + I(x.tr^2))
  yhat2 <- J2$coef[1] + J2$coef[2] * x.te + J2$coef[3] * x.te^2
  e2[k] <- mean((y.te - yhat2)^2)

  ## Model 3: log(y) ~ x
  J3 <- lm(log(y.tr) ~ x.tr)
  logyhat3 <- J3$coef[1] + J3$coef[2] * x.te
  yhat3 <- exp(logyhat3)
  e3[k] <- mean((y.te - yhat3)^2)

  ## Model 4: log(y) ~ log(x)
  J4 <- lm(log(y.tr) ~ log(x.tr))
  logyhat4 <- J4$coef[1] + J4$coef[2] * log(x.te)
  yhat4 <- exp(logyhat4)
  e4[k] <- mean((y.te - yhat4)^2)
}

c(mean(e1^2), mean(e2^2), mean(e3^2), mean(e4^2))

## -----------------------------------------------------------------------------
set.seed(101)

mc_boot_ci <- function(m = 1e4, n = 50, B = 1000, mu = 0) {
  res <- list(
    normal = c(hit = 0, missL = 0, missR = 0),
    basic  = c(hit = 0, missL = 0, missR = 0),
    perc   = c(hit = 0, missL = 0, missR = 0)
  )

  for (i in 1:m) {
    x = rnorm(n, mean = mu, sd = 1)
    theta_hat = mean(x)
    boot_means = replicate(B, mean(sample(x, n, replace = TRUE)))
    se_boot = sd(boot_means)

    q = quantile(boot_means, c(0.025, 0.975))
    qL = q[1]; qU = q[2]

    ## (1) bootstrap-normal CI
    z = qnorm(0.975)
    ciN_L = theta_hat - z * se_boot
    ciN_U = theta_hat + z * se_boot

    ## (2) basic CI
    ciB_L = 2*theta_hat - qU
    ciB_U = 2*theta_hat - qL

    ## (3) percentile CI
    ciP_L = qL
    ciP_U = qU
    update <- function(ciL, ciU, name) {
      if (mu < ciL) {
        res[[name]]["missL"] <<- res[[name]]["missL"] + 1
      } else if (mu > ciU) {
        res[[name]]["missR"] <<- res[[name]]["missR"] + 1
      } else {
        res[[name]]["hit"]   <<- res[[name]]["hit"]   + 1
      }
    }
    update(ciN_L, ciN_U, "normal")
    update(ciB_L, ciB_U, "basic")
    update(ciP_L, ciP_U, "perc")
  }
  out <- sapply(res, function(v) v / m)
  rownames(out) <- c("coverage", "miss_left", "miss_right")
  return(out)
}

res = mc_boot_ci()
round(res, 4)

## -----------------------------------------------------------------------------
set.seed(101)

lambda_true = 2
n_vec = c(5, 10, 20)
B = 1000
m = 1000

res <- data.frame(
  n = n_vec,
  boot_bias_mean = NA,
  theory_bias = NA,
  boot_se_mean = NA,
  theory_se = NA
)

for (k in seq_along(n_vec)) {
  n = n_vec[k]
  boot_bias = numeric(m)
  boot_se = numeric(m)

  for (j in 1:m){
    x = rexp(n, lambda_true)
    lambda_hat = 1/mean(x)

    # bootstrap
    boot_lam = replicate(B, {
      xb = sample(x, n, replace = TRUE)
      1 / mean(xb)
    })

    boot_bias[j] = mean(boot_lam) - lambda_hat
    boot_se[j] = sd(boot_lam)
  }

  # Monte Carlo
  res$boot_bias_mean[k] = mean(boot_bias)
  res$boot_se_mean[k] = mean(boot_se)

  # 理论
  res$theory_bias[k] = lambda_true / (n - 1)
  res$theory_se[k] = lambda_true * n / ((n - 1) * sqrt(n - 2))
}

print(res)

## -----------------------------------------------------------------------------
attach(chickwts)
x = sort(weight[feed == "soybean"])
y = sort(weight[feed == "linseed"])
detach(chickwts)

cvm2 <- function(x, y){
  n = length(x)
  m = length(y)
  N = n + m
  rk = rank(c(x, y), ties.method = "average")
  r = rk[1:n]
  s = rk[(n+1):N]
  U  = n * sum((r - 1:n)^2) + m * sum((s - 1:m)^2)
  U/(n*m*(n+m)) - (4*n*m - 1)/(6*(n+m))
}

options(warn =-1)
W0 = cvm2(x, y)

## 置换检验
R = 999                  
z = c(x, y) 
n = length(x)
m = length(y)
K = 1:(n + m)
W = numeric(R)
for(i in 1:R){
  k  = sample(K, size = n, replace = FALSE) 
  x1 = z[k]
  y1 = z[-k]                                
  W[i] = cvm2(x1, y1)
}
p = mean(c(W0, W) >= W0)
options(warn = 0)
p

## -----------------------------------------------------------------------------
count5test = function(x, y) {
 X = x - mean(x)
 Y = y - mean(y)
 outx = sum(X > max(Y)) + sum(X < min(Y))
 outy = sum(Y > max(X)) + sum(Y < min(X))
 return(max(outx,outy))
}

perm_count5 <- function(x, y, R){
  z = c(x, y)
  n = length(x)
  T0 = count5test(x, y)
  Tstar = numeric(R)
  for (r in 1:R) {
    idx = sample.int(length(z), n)
    x1 = z[idx]
    y1 = z[-idx]
    Tstar[r] = count5test(x1, y1)
  }
  p = (1 + sum(Tstar >= T0)) / (R + 1)
  as.integer(p <= alpha)
}

set.seed(123)
n1 = 20; n2 = 30
m = 10000
alpha = 0.05
alphahat_correct1 = mean(replicate(m, {
  x = rnorm(n1)
  y = rnorm(n2)
  perm_count5(x, y, R = 999)
}))

n3 = 50
m = 10000
alpha = 0.05
alphahat_correct2 = mean(replicate(m, {
  x = rnorm(n1)
  y = rnorm(n3)
  perm_count5(x, y, R = 999)
}))

cat(sprintf("样本量为20/30下的经验 Type I 错误率: %.4f\n", alphahat_correct1))
cat(sprintf("样本量为20/50下的经验 Type I 错误率: %.4f\n", alphahat_correct2))

## ----fig.width=10, fig.height=8-----------------------------------------------
dlap <- function(x) 0.5 * exp(-abs(x)) 
rw_metropolis_laplace <- function(n, x0 = 25, sigma = 1) {
  x <- numeric(n); x[1] <- x0; acc <- 0L
  for (t in 2:n) {
    y = x[t-1] + rnorm(1, 0, sigma)        
    log_alpha = -abs(y) + abs(x[t-1])
    if (log(runif(1)) < log_alpha){ 
      x[t] = y
      acc = acc + 1L 
      }else x[t] <- x[t-1]
  }
  list(chain = x, acc_rate = acc/(n-1), sigma = sigma)
}

set.seed(123)
n_iter = 2000
sigmas = c(0.05, 0.5, 2, 16)

runs <- lapply(sigmas, function(s) rw_metropolis_laplace(n_iter, x0 = 25, sigma = s))
cat("Acceptance rates:\n")
for (r in runs) cat(sprintf("  sigma = %-4g : %6.2f%%\n", r$sigma, 100*r$acc_rate))

op <- par(mfrow = c(2,2), mar = c(3.5, 3.5, 2, 1), mgp = c(2.1, 0.7, 0))
for (r in runs) {
  xx = r$chain
  yl = range(xx) 
  plot(xx, type = "l", xlab = "", ylab = "X", main = "",
       ylim = yl)
  usr = par("usr")
  rect(xleft = usr[1], ybottom = -2, xright = usr[2], ytop = 2,
       col = "gray90", border = NA)
  lines(xx, type = "l")      
  abline(h = c(-2, 0, 2), col = "gray70", lwd = 1)
  mtext(bquote(sigma == .(r$sigma)), side = 1, line = 2) 
}
par(op)

## ----fig.width=10, fig.height=8-----------------------------------------------
dbetabinom <- function(x, n, a, b) {
  choose(n, x) * beta(x + a, n - x + b) / beta(a, b)
}

true <- function(n = 25, a = 2, b = 3) {
  list(n = n, a = a, b = b,
       Ey = a / (a + b),
       Vy = a*b / ((a + b)^2 * (a + b + 1)))
}

gibbs_xy <- function(truth, n_iter = 20000, burnin = 2000, y0 = 0.5, seed = 123){
  set.seed(seed)
  n = truth$n
  a = truth$a
  b = truth$b
  X = integer(n_iter)
  Y = numeric(n_iter)
  Y[1] = y0
  X[1] = rbinom(1, n, Y[1])
  for (t in 2:n_iter) {
    X[t] = rbinom(1, n, Y[t-1])
    Y[t] = rbeta(1, X[t] + a, n - X[t] + b)
  }
  list(X = X, Y = Y, burnin = burnin, truth = truth, seed = seed)
}

report_xy <- function(sim) {
  n = sim$truth$n 
  a = sim$truth$a
  b = sim$truth$b
  X = sim$X
  Y = sim$Y
  keep = (sim$burnin + 1):length(X)
  xk = X[keep]
  yk = Y[keep]


  op <- par(mfrow = c(2, 2), mar = c(3.5, 3.5, 2, 1), mgp = c(2.1, 0.7, 0))
  on.exit(par(op), add = TRUE)

  plot(X[1:2000], type = "l", xlab = "iteration", ylab = "X")
  plot(Y[1:2000], type = "l", xlab = "iteration", ylab = "Y")

  tabx = table(factor(xk, levels = 0:n))
  barplot(tabx / sum(tabx), space = 0, xlab = "x",
          main = "X marginal Dist")
  lines(0:n, dbetabinom(0:n, n, a, b), lwd = 2)

  hist(yk, breaks = 50, freq = FALSE, main = "Y marginal Dist",
       xlab = "y")
  curve(dbeta(x, a, b), add = TRUE, lwd = 2)
  par(op)
}

truth = true(n = 25, a = 2, b = 3)
sim = gibbs_xy(truth, n_iter = 20000, burnin = 2000, y0 = 0.5, seed = 11)
report_xy(sim)

## -----------------------------------------------------------------------------
Gelman.Rubin <- function(psi) {
  psi = as.matrix(psi)
  n = ncol(psi)
  k = nrow(psi)
  psi.means = rowMeans(psi)
  B = n * var(psi.means) 
  psi.w = apply(psi, 1, "var")
  W = mean(psi.w)          
  v.hat = W * (n - 1) / n + (B / n)
  r.hat = v.hat / W  
  return(r.hat)
}

run_gr_for_sigma <- function(sigma = 2, k = 4, n = 15000, b = 1000,
                             x0 = c(-10, -5, 5, 10), seed = 123){
  set.seed(seed)
  X = matrix(0, nrow = k, ncol = n)
  for (i in 1:k) {
    X[i, ] = rw_metropolis_laplace(n, x0 = x0[i], sigma = sigma)$chain
  }
  psi = t(apply(X, 1, cumsum))
  for (i in 1:nrow(psi)){ 
    psi[i, ] = psi[i, ] / (1:ncol(psi))
  }
  rhat = rep(NA_real_, n)
  for (j in (b+1):n){ 
    rhat[j] = Gelman.Rubin(psi[, 1:j])
  }
  hit = which(rhat < 1.2)[1]
  plot((b + 1):n, rhat[(b + 1):n], type = "l",
       xlab = "iteration (j)", ylab = expression(hat(R)),
       main = bquote("Gelman–Rubin  " ~ sigma == .(sigma)))
  abline(h = 1.2, lty = 2)
  if (!is.na(hit)) abline(v = hit, col = "gray40", lty = 3)

  invisible(list(X = X, psi = psi, rhat = rhat, hit = hit,
                 sigma = sigma, k = k, n = n, b = b))
}

## -----------------------------------------------------------------------------
sigmas = c(0.1, 0.5, 2, 16)
res_list = lapply(sigmas, function(s) run_gr_for_sigma(sigma = s, k = 4, n = 20000, b = 1000))

for (i in seq_along(sigmas)) {
  cat(sprintf("sigma = %-5g : first iter with R-hat < 1.2 is %s\n",
              sigmas[i],
              ifelse(is.na(res_list[[i]]$hit), "not reached", res_list[[i]]$hit)))
}

## ----fig.width=10, fig.height=8-----------------------------------------------
gr_for_gibbs <- function(truth,
                         k = 4, n_iter = 15000, burnin = 1000,
                         variable = c("Y","X"),
                         y0_vec = c(0.05, 0.25, 0.75, 0.95),
                         seeds = 1:4,
                         plot_it = TRUE) {
  variable = match.arg(variable)
  stopifnot(length(y0_vec) >= k, length(seeds) >= k)
  
  Mat = matrix(NA, nrow = k, ncol = n_iter)
  for (i in 1:k) {
    simi = gibbs_xy(truth, n_iter = n_iter, burnin = 0,
                     y0 = y0_vec[i], seed = seeds[i])
    Mat[i, ] <- if (variable == "Y") simi$Y else simi$X
  }

  psi = t(apply(Mat, 1, cumsum))
  for (i in 1:nrow(psi)) psi[i, ] = psi[i, ] / (1:ncol(psi))
  rhat = rep(NA_real_, n_iter)
  for (j in (burnin + 1):n_iter) rhat[j] = Gelman.Rubin(psi[, 1:j])
  hit = which(rhat < 1.2)[1]

  if (plot_it) {
    op = par(mfrow = c(2, 2), mar = c(3.5, 3.5, 2, 1), mgp = c(2.1, 0.7, 0))
    for (i in 1:k) {
      plot(psi[i, (burnin + 1):n_iter], type = "l",
           xlab = "iteration (j)", ylab = expression(psi),
           main = paste(variable, ": Chain", i))
    }
    par(op)

    plot((burnin + 1):n_iter, rhat[(burnin + 1):n_iter], type = "l",
         xlab = "iteration (j)", ylab = expression(hat(R)),
         main = paste0("Gelman–Rubin for ", variable))
    abline(h = 1.2, lty = 2)
    if (!is.na(hit)) abline(v = hit, col = "gray40", lty = 3)
  }

  list(psi = psi, rhat = rhat, hit = hit,
       n_iter = n_iter, burnin = burnin, variable = variable)
}

truth = true(n = 25, a = 2, b = 3)

y0s = seq(0.1, 0.9, length.out = 4)

grY = gr_for_gibbs(truth, k = 4, n_iter = 15000, burnin = 1000,
                   variable = "Y", y0_vec = y0s, seeds = 11:14)

grX = gr_for_gibbs(truth, k = 4, n_iter = 15000, burnin = 1000,
                   variable = "X", y0_vec = y0s, seeds = 21:24)

## -----------------------------------------------------------------------------
solve_alpha <- function(f0,
                        N = 1e6,
                        b1 = 1, b2 = 1, b3 = -1){
  interval = c(-20, 20)
  set.seed(123)
  x1 = rpois(N, 1)
  x2 = rexp(N, 1)
  x3 = rbinom(N, 1, 0.5)

  g <- function(alpha) {
    eta = alpha + b1 * x1 + b2 * x2 + b3 * x3
    p = 1 / (1 + exp(-eta))
    mean(p) - f0
  }

  solution = uniroot(g, interval = interval)
  alpha = solution$root

  invisible(list(alpha = alpha, solution = solution,
                 f0 = f0, b = c(b1 = b1, b2 = b2, b3 = b3),
                 N = N))
}

## -----------------------------------------------------------------------------
f0_values = c(0.1, 0.01, 0.001, 0.0001)
alphas = sapply(f0_values, function(f0)
  solve_alpha(f0 = f0, N = 1e6,
              b1 = 1, b2 = 1, b3 = -1)$alpha)

## -----------------------------------------------------------------------------
out_tab = data.frame(f0 = f0_values, alpha = alphas)
print(out_tab)

## -----------------------------------------------------------------------------
plot(-log(out_tab$f0), out_tab$alpha, type = "b", pch = 19,
     xlab = expression(-log(f[0])), ylab = "alpha")

## ----warning=FALSE------------------------------------------------------------
library(boot)

a = c(-4, -2, -9)
A1 = rbind(
  c(2, 1, 1),
  c(1, -1, 3)
)
b1 = c(2, 3)

sol = simplex(a = a, A1 = A1, b1 = b1, maxi = TRUE)
sol

## -----------------------------------------------------------------------------
u <- c(11, 8, 27, 13, 16, 0, 23, 10, 24, 2)
v <- c(12, 9, 28, 14, 17, 1, 24, 11, 25, 3)
n <- length(u)
lambda_init <- 0.1
tol <- 1e-8
max_iter <- 100

log_lik_func <- function(lambda, u, v) {
  sum(log(exp(-lambda * u) - exp(-lambda * v)))
}

newton_raphson <- function(u, v, lambda_start, tol, max_iter) {
  lambda <- lambda_start
  path <- c(lambda)
  for (i in 1:max_iter) {
    delta <- v - u
    term1 <- -u
    term2 <- delta / (exp(lambda * delta) - 1)
    score <- sum(term1 + term2)
    hessian <- -sum((delta^2 * exp(lambda * delta)) / (exp(lambda * delta) - 1)^2)
    lambda_new <- lambda - score / hessian
    path <- c(path, lambda_new)
    if (abs(lambda_new - lambda) < tol) {
      return(list(est = lambda_new, iter = i, path = path, converged = TRUE))
    }
    lambda <- lambda_new
  }
  return(list(est = lambda, iter = max_iter, path = path, converged = FALSE))
}

em_algorithm <- function(u, v, lambda_start, tol, max_iter) {
  lambda <- lambda_start
  path <- c(lambda)
  n <- length(u)
  for (i in 1:max_iter) {
    denom <- exp(-lambda * u) - exp(-lambda * v)
    num_part <- u * exp(-lambda * u) - v * exp(-lambda * v)
    expected_X <- (1 / lambda) + (num_part / denom)
    lambda_new <- n / sum(expected_X)
    path <- c(path, lambda_new)
    if (abs(lambda_new - lambda) < tol) {
      return(list(est = lambda_new, iter = i, path = path, converged = TRUE))
    }
    lambda <- lambda_new
  }
  return(list(est = lambda, iter = max_iter, path = path, converged = FALSE))
}

res_nr <- newton_raphson(u, v, lambda_init, tol, max_iter)
res_em <- em_algorithm(u, v, lambda_init, tol, max_iter)

cat(sprintf("Newton-Raphson (MLE): 估计值 = %.4f, 迭代次数 = %d\n", res_nr$est, res_nr$iter))
cat(sprintf("EM Algorithm (EM):    估计值 = %.4f, 迭代次数 = %d\n", res_em$est, res_em$iter))

diff_est <- abs(res_nr$est - res_em$est)
cat(sprintf("两者结果差异: %.4f \n", diff_est))

true_mle <- res_nr$est
nr_error <- abs(res_nr$path - true_mle) + 1e-20
em_error <- abs(res_em$path - true_mle) + 1e-20

max_len <- max(length(nr_error), length(em_error))
nr_plot <- c(nr_error, rep(NA, max_len - length(nr_error)))
em_plot <- c(em_error, rep(NA, max_len - length(em_error)))

plot(1:max_len, log10(em_plot), type = "o", col = "blue", pch = 16,
     xlab = "Iteration", ylab = "Log10(Error)",
     main = "Convergence Speed (Log Scale)")
lines(1:max_len, log10(nr_plot), type = "o", col = "red", pch = 17)
legend("topright", legend = c("EM Algo (Linear)", "Newton-Raphson (Quadratic)"), 
       col = c("blue", "red"), pch = c(16, 17), lty = 1)

## -----------------------------------------------------------------------------
vmap <- function(FUN, ..., FUN.VALUE) {
  args = list(...)
  n = length(args[[1]])

  vapply(
    X = seq_len(n),
    FUN = function(i) {
      ith_args = Map(`[`, args, i)
      do.call(FUN, ith_args)
    },
    FUN.VALUE = FUN.VALUE
  )
}

## -----------------------------------------------------------------------------
x = 1:5
y = 6:10

vmap(`+`, x, y, FUN.VALUE = numeric(1))

## -----------------------------------------------------------------------------
fast_chisq <- function(x, y) {
  stopifnot(is.numeric(x), is.numeric(y))
  stopifnot(length(x) == length(y))
  
  tab = table(x, y)
  rs = rowSums(tab)
  cs = colSums(tab)
  n = sum(tab)
  E = outer(rs, cs, "*") / n
  X2 = sum((tab - E)^2 / E)
  return(X2)
}

## -----------------------------------------------------------------------------
fast_table <- function(x, y) {
  ix = x - min(x) + 1
  iy = y - min(y) + 1
  nx = max(ix)
  ny = max(iy)
  tab = matrix(0L, nrow = nx, ncol = ny)
  for (i in seq_along(ix)) {
    tab[ix[i], iy[i]] <- tab[ix[i], iy[i]] + 1L
  }
  return(tab)
}

## -----------------------------------------------------------------------------
faster_chisq <- function(x, y) {
  tab = fast_table(x, y)
  rs = rowSums(tab)
  cs = colSums(tab)
  n = sum(tab)
  E = outer(rs, cs, "*") / n
  X2 = sum((tab - E)^2 / E)
  return(X2)
}

## -----------------------------------------------------------------------------
n_obs = c(125, 18, 20, 34)

logpost_theta <- function(theta, n = n_obs) {
  if (theta <= 0 || theta >= 1) return(-Inf)
  p1 = 0.5 + theta / 4
  p2 = (1 - theta) / 4
  p3 = (1 - theta) / 4
  p4 = theta / 4
  l = sum(n * log(c(p1, p2, p3, p4)))
  if (any(c(p1, p2, p3, p4) <= 0)) return(-Inf)
  return(l)
}

rw_theta_R <- function(N, theta0, w, n = n_obs) {
  theta = numeric(N)
  theta[1] = theta0
  accept = 0L
  for (i in 2:N) {
    y = theta[i - 1] + runif(1, -w, w)
    if (y <= 0 || y >= 1) { 
      theta[i] = theta[i - 1]
      next
    }
    logalpha = logpost_theta(y, n) - logpost_theta(theta[i - 1], n)
    if (log(runif(1)) < logalpha) {
      theta[i] = y
      accept = accept + 1L
    } else {
      theta[i] = theta[i - 1]
    }
  }
  list(theta = theta, acc_rate = accept / (N - 1))
}

## -----------------------------------------------------------------------------
library(Rcpp)
Rcpp::sourceCpp("E:/code/rw_theta.cpp")

## -----------------------------------------------------------------------------
set.seed(123)
N = 100000
w = 0.2
n = n_obs
th0 = 0.3

## R
out_R = rw_theta_R(N, th0, w, n)
theta_R = out_R$theta

## Rcpp
out_cpp = rw_theta_cpp(N, th0, w, n)
theta_cpp = out_cpp$theta

burn = 1000
theta_Rb = theta_R[-(1:burn)]
theta_cppb = theta_cpp[-(1:burn)]

qqplot(theta_Rb, theta_cppb,
       xlab = "R",
       ylab = "Rcpp",
       main = "QQ-plot: R vs Rcpp")
abline(0, 1, col = 2, lwd = 2)

## -----------------------------------------------------------------------------
library(microbenchmark)
set.seed(123)
N = 100000
bench <- microbenchmark(
  R = rw_theta_R(N, th0, w, n),
  Rcpp = rw_theta_cpp(N, th0, w, n),
  times = 20L
)
print(bench)

