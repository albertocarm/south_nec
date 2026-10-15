library(rstan)
library(survival)
library(loo)
library(ggplot2)
library(patchwork)

rstan_options(auto_write = TRUE)
options(mc.cores = 4)
dir.create("results", showWarnings = FALSE)

d <- read.csv("data/southnec.csv")
s_med <- median(d$chemo_start[d$chemo_type == "adjuvant" & d$chemo_start > 0], na.rm = TRUE)
KI_KNOTS <- c(30, 70, 90)
AGE_KNOTS <- c(49, 71, 83)
model <- stan_model("stan/cure_model.stan")

# second basis of a restricted cubic spline with three knots
rcs2 <- function(x, k) {
  (pmax(x - k[1], 0)^3 - pmax(x - k[2], 0)^3 * (k[3] - k[1]) / (k[3] - k[2]) +
     pmax(x - k[3], 0)^3 * (k[2] - k[1]) / (k[3] - k[2])) / (k[3] - k[1])^2
}

standata <- function(x, event, ki = "spline", age = "spline", interaction = FALSE, sd_clono = 2.5, sd_kin = 0.25) {
  ko <- !is.na(x$ki67)
  ao <- !is.na(x$age)
  kb <- rcs2(x$ki67, KI_KNOTS)
  ab <- rcs2(x$age, AGE_KNOTS)
  k_std <- c(mean(x$ki67[ko]), sd(x$ki67[ko]), mean(kb[ko]), sd(kb[ko]))
  a_std <- c(mean(x$age[ao]), sd(x$age[ao]), mean(ab[ao]), sd(ab[ao]))
  k1 <- ifelse(ko, (x$ki67 - k_std[1]) / k_std[2], 0)
  k2 <- ifelse(ko, (kb - k_std[3]) / k_std[4], 0)
  a1 <- ifelse(ao, (x$age - a_std[1]) / a_std[2], 0)
  a2 <- ifelse(ao, (ab - a_std[3]) / a_std[4], 0)
  X <- scale(cbind(x$stage == "II", x$stage == "III", x$colorectal))
  if (ki == "spline") {
    X <- cbind(X, k1, k2); Z <- cbind(k1, k2); ck <- c(4, 5, 1, 2, 0)
  } else {
    X <- cbind(X, k1); Z <- cbind(k1); ck <- c(4, 0, 1, 0, 0)
  }
  if (age == "spline") {
    X <- cbind(X, a1, a2); ca <- ncol(X) - 1:0
  } else {
    X <- cbind(X, a1); ca <- c(ncol(X), 0)
  }
  W <- matrix(1 / sd(x$chemo), nrow(x), 1)
  if (interaction) {
    W <- cbind(W, k1); ck[5] <- 2
  }
  s <- x$chemo_start
  path <- ifelse(is.na(s) | s >= x$time, 0, ifelse(s <= 0, 1, 2))
  list(N = nrow(x), K = ncol(X), Kc = ncol(W), L = ncol(Z), X = unname(X), W = W, Z = unname(Z), V = W[, 1, drop = FALSE],
       t = x$time, d = x[[event]], path = path, s = ifelse(path == 2, s, 0),
       sd_beta = as.array(rep(sd_clono, ncol(X))), sd_delta = as.array(rep(sd_clono, ncol(W))), sd_gamma = as.array(rep(sd_kin, ncol(Z))),
       sd_eta = sd_kin, sd_int = 2.5,
       Nk = sum(!ko), ik = as.array(which(!ko)), k_lo = as.array(x$ki67_min[!ko]), k_hi = as.array(rep(100, sum(!ko))),
       k_knots = KI_KNOTS, k_std = k_std, ck = ck,
       Na = sum(!ao), ia = as.array(which(!ao)), a_lo = as.array(rep(18, sum(!ao))), a_hi = as.array(rep(100, sum(!ao))),
       a_knots = AGE_KNOTS, a_std = a_std, ca = ca, sd_chemo = sd(x$chemo))
}

run <- function(x) {
  sampling(model, data = x, chains = 4, iter = 4000, warmup = 2000, seed = 42, refresh = 0,
           control = list(adapt_delta = 0.95, max_treedepth = 12), pars = c("Xf", "Wf", "Zf"), include = FALSE)
}

# chemotherapy ratios, cure probabilities and convergence
cure_summary <- function(f, x) {
  p <- extract(f)
  cl <- p$delta[, 1] / x$sd_chemo
  kn <- p$eta[, 1] / x$sd_chemo
  Fs <- 1 - exp(-(s_med / exp(p$gamma0))^p$alpha)
  c0 <- exp(-exp(p$beta0))
  c1 <- exp(-(exp(p$beta0) * Fs + exp(p$beta0 + cl) * (1 - Fs)))
  or <- (c1 / (1 - c1)) / (c0 / (1 - c0))
  s <- summary(f)$summary
  s <- s[is.finite(s[, "sd"]) & s[, "sd"] > 0, ]
  q <- function(v, tr = identity, k = 1, fmt = "%.2f") {
    v <- k * tr(quantile(v, c(0.5, 0.025, 0.975)))
    sprintf(paste0(fmt, " (", fmt, "-", fmt, ")"), v[1], v[2], v[3])
  }
  pd <- function(v) sprintf("%.1f%%", 100 * max(mean(v < 0), mean(v > 0)))
  data.frame(clonogenic_ratio = q(cl, exp), clonogenic_pd = pd(cl),
             kinetic_ratio = q(kn, exp), kinetic_pd = pd(kn),
             cure_fraction = q(rowMeans(p$cure), k = 100, fmt = "%.1f"),
             cure_untreated = q(c0, k = 100, fmt = "%.1f"), cure_chemo = q(c1, k = 100, fmt = "%.1f"),
             cure_or = q(or), cure_or_above_1 = sprintf("%.1f%%", 100 * mean(or > 1)),
             alpha = q(p$alpha), correlation = sprintf("%.2f", cor(p$delta[, 1], p$eta[, 1])),
             rhat = sprintf("%.3f", max(s[, "Rhat"])), n_eff = round(min(s[, "n_eff"])), divergent = get_num_divergent(f))
}
