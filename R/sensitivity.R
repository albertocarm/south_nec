source("R/setup.R")

# Priors (Supplementary Table A1)
pri <- data.frame(prior = paste0("P", 1:6), clonogenic = c(2.5, 1, 2.5, 0.25, 0.25, 5), kinetic = c(0.25, 1, 2.5, 0.25, 2.5, 5))
a1 <- NULL
for (i in 1:6) {
  x <- standata(d, "efs", sd_clono = pri$clonogenic[i], sd_kin = pri$kinetic[i])
  a1 <- rbind(a1, data.frame(pri[i, ], cure_summary(run(x), x)))
}
print(a1, row.names = FALSE)
write.csv(a1, "results/priors.csv", row.names = FALSE)

# Linear Ki-67 or age (Supplementary Table A10)
a10 <- NULL
for (e in c("efs", "rfs")) {
  for (f in list(c("spline", "spline"), c("linear", "spline"), c("spline", "linear"))) {
    x <- standata(d, e, ki = f[1], age = f[2])
    a10 <- rbind(a10, data.frame(endpoint = toupper(e), ki67 = f[1], age = f[2], cure_summary(run(x), x)))
  }
}
print(a10, row.names = FALSE)
write.csv(a10, "results/functional_form.csv", row.names = FALSE)

# Truncated follow-up (EFS)
tr <- NULL
for (tc in c(0.9 * max(d$time), 0.8 * max(d$time), 0.7 * max(d$time), 120)) {
  y <- d
  y$efs[y$time > tc] <- 0
  y$time <- pmin(y$time, tc)
  x <- standata(y, "efs")
  tr <- rbind(tr, data.frame(months = round(tc), events = sum(y$efs), cure_summary(run(x), x)))
}
print(tr, row.names = FALSE)
write.csv(tr, "results/truncation.csv", row.names = FALSE)
