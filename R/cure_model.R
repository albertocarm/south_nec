source("R/setup.R")

# Main models
x <- list(EFS = standata(d, "efs"), RFS = standata(d, "rfs"))
fit <- lapply(x, run)
res <- do.call(rbind, lapply(names(fit), function(e) data.frame(endpoint = e, cure_summary(fit[[e]], x[[e]]))))
print(t(res))
write.csv(res, "results/cure_models.csv", row.names = FALSE)

# Table 4 (EFS) and Supplementary Table A5 (RFS)
for (e in names(fit)) {
  p <- extract(fit[[e]])
  dk <- c(40 / x[[e]]$k_std[2], (rcs2(80, KI_KNOTS) - rcs2(40, KI_KNOTS)) / x[[e]]$k_std[4])
  da <- c(18 / x[[e]]$a_std[2], (rcs2(79, AGE_KNOTS) - rcs2(61, AGE_KNOTS)) / x[[e]]$a_std[4])
  v <- list("Stage II (vs I)" = p$beta[, 1] / sd(d$stage == "II"),
            "Stage III (vs I)" = p$beta[, 2] / sd(d$stage == "III"),
            "Site: colorectal (vs other)" = p$beta[, 3] / sd(d$colorectal),
            "Ki-67 (80% vs 40%)" = p$beta[, 4] * dk[1] + p$beta[, 5] * dk[2],
            "Age (79 vs 61 years)" = p$beta[, 6] * da[1] + p$beta[, 7] * da[2],
            "Perioperative chemotherapy" = p$delta[, 1] / x[[e]]$sd_chemo,
            "Kinetic: Ki-67 (80% vs 40%)" = p$gamma[, 1] * dk[1] + p$gamma[, 2] * dk[2],
            "Kinetic: perioperative chemotherapy" = p$eta[, 1] / x[[e]]$sd_chemo)
  tab <- data.frame(covariate = names(v),
                    coefficient = sprintf("%.3f", sapply(v, median)),
                    cri = sprintf("%.3f, %.3f", sapply(v, quantile, 0.025), sapply(v, quantile, 0.975)),
                    exp_coef = sprintf("%.2f", exp(sapply(v, median))),
                    pd = sprintf("%.1f%%", 100 * sapply(v, function(z) max(mean(z > 0), mean(z < 0)))))
  print(tab, row.names = FALSE)
  write.csv(tab, paste0("results/coefficients_", e, ".csv"), row.names = FALSE)
}

# Diagnostics
for (e in names(fit)) {
  p <- extract(fit[[e]])
  ll <- extract_log_lik(fit[[e]], merge_chains = FALSE)
  k <- loo(ll, r_eff = relative_eff(exp(ll)))$diagnostics$pareto_k
  cat(e, "| Ki-67 cross-correlation", round(max(abs(cor(p$beta[, 4:5], p$gamma[, 1:2]))), 2), "| Pareto k > 0.7:", sum(k > 0.7), "\n")
}

# Posterior predictive check (Supplementary Figure A3)
grid <- seq(0, 180, by = 3)
ppc <- NULL
km <- NULL
for (e in names(fit)) {
  p <- extract(fit[[e]], c("lt0", "lt1", "ll0", "ll1", "alpha"))
  n <- x[[e]]$N
  Tm <- matrix(grid, n, length(grid), byrow = TRUE)
  S <- matrix(x[[e]]$s, n, length(grid))
  P <- matrix(x[[e]]$path, n, length(grid))
  m <- sapply(seq_along(p$alpha), function(i) {
    a <- p$alpha[i]
    th0 <- exp(p$lt0[i, ]); th1 <- exp(p$lt1[i, ]); l0 <- exp(p$ll0[i, ]); l1 <- exp(p$ll1[i, ])
    H0 <- th0 * (1 - exp(-(Tm / l0)^a))
    H1 <- th1 * (1 - exp(-(Tm / l1)^a))
    Hs <- th0 * (1 - exp(-(S / l0)^a)) + th1 * exp(-(S / l0)^a) * (1 - exp(-pmax((Tm / l1)^a - (S / l1)^a, 0)))
    colMeans(exp(-ifelse(P == 1, H1, ifelse(P == 2 & Tm > S, Hs, H0))))
  })
  s <- survfit(Surv(d$time, d[[tolower(e)]]) ~ 1)
  ppc <- rbind(ppc, data.frame(endpoint = e, month = grid,
                               kaplan_meier = summary(s, times = grid, extend = TRUE)$surv,
                               model = apply(m, 1, median),
                               lower = apply(m, 1, quantile, 0.025),
                               upper = apply(m, 1, quantile, 0.975)))
  km <- rbind(km, data.frame(endpoint = e, month = c(0, s$time), surv = c(1, s$surv)))
}
write.csv(ppc, "results/ppc_curves.csv", row.names = FALSE)
print(ppc[ppc$month %in% c(12, 24, 60, 120), ], row.names = FALSE, digits = 2)
g <- ggplot(ppc, aes(month)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), fill = "#2E6F9E", alpha = 0.18) +
  geom_line(aes(y = model), colour = "#2E6F9E", linewidth = 0.9) +
  geom_step(data = km, aes(month, surv), linewidth = 0.6) +
  facet_wrap(~endpoint) + coord_cartesian(xlim = c(0, 180), ylim = c(0, 1)) +
  labs(x = "Months from surgery", y = "Survival probability") + theme_minimal(base_size = 11)
ggsave("results/figure_A3.pdf", g, width = 9.5, height = 4.4)

# Ki-67 x chemotherapy (Supplementary Table A6)
xi <- list(EFS = standata(d, "efs", ki = "linear", interaction = TRUE), RFS = standata(d, "rfs", ki = "linear", interaction = TRUE))
fi <- lapply(xi, run)
for (e in names(fi)) {
  p <- extract(fi[[e]])
  r <- sapply(c(50, 90), function(k) p$delta[, 1] / xi[[e]]$sd_chemo + p$delta[, 2] * (k - xi[[e]]$k_std[1]) / xi[[e]]$k_std[2])
  cat(e, sprintf("| interaction pd %.1f%%", 100 * max(mean(p$delta[, 2] > 0), mean(p$delta[, 2] < 0))),
      sprintf("| clonogenic ratio, Ki-67 %d%%: %.2f (%.2f-%.2f)", c(50, 90), exp(apply(r, 2, median)),
              exp(apply(r, 2, quantile, 0.025)), exp(apply(r, 2, quantile, 0.975))), "\n")
}

# Figure 2 (RFS)
ki <- seq(20, 100, length.out = 100)
p <- extract(fit$RFS)
z1 <- (ki - x$RFS$k_std[1]) / x$RFS$k_std[2]
z2 <- (rcs2(ki, KI_KNOTS) - x$RFS$k_std[3]) / x$RFS$k_std[4]
th <- exp(c(p$beta0) + outer(p$beta[, 4], z1) + outer(p$beta[, 5], z2))
mt <- exp(c(p$gamma0) + outer(p$gamma[, 1], z1) + outer(p$gamma[, 2], z2)) * log(2)^(1 / c(p$alpha))
abc <- rbind(data.frame(panel = "theta", ki67 = ki, t(apply(th, 2, quantile, c(0.5, 0.025, 0.975)))),
             data.frame(panel = "cure", ki67 = ki, 100 * t(apply(exp(-th), 2, quantile, c(0.5, 0.025, 0.975)))),
             data.frame(panel = "time", ki67 = ki, t(apply(mt, 2, quantile, c(0.5, 0.025, 0.975)))))
names(abc)[3:5] <- c("med", "lo", "hi")
p <- extract(fi$RFS)
z <- (ki - xi$RFS$k_std[1]) / xi$RFS$k_std[2]
lt0 <- c(p$beta0) + outer(p$beta[, 4], z)
lt1 <- lt0 + p$delta[, 1] / xi$RFS$sd_chemo + outer(p$delta[, 2], z)
Fs <- 1 - exp(-(s_med / exp(c(p$gamma0) + outer(p$gamma[, 1], z)))^c(p$alpha))
c0 <- exp(-exp(lt0))
c1 <- exp(-(exp(lt0) * Fs + exp(lt1) * (1 - Fs)))
pts <- round(seq(1, length(ki), length.out = 15))
de <- data.frame(ki67 = ki[pts], t(apply((c1 / (1 - c1)) / (c0 / (1 - c0)), 2, quantile, c(0.5, 0.025, 0.975)))[pts, ])
grp <- data.frame(ki67 = ki[pts], group = rep(c("Periop. chemotherapy", "No periop. chemotherapy"), each = length(pts)),
                  100 * rbind(t(apply(c1, 2, quantile, c(0.5, 0.025, 0.975)))[pts, ], t(apply(c0, 2, quantile, c(0.5, 0.025, 0.975)))[pts, ]))
names(de)[2:4] <- names(grp)[3:5] <- c("med", "lo", "hi")
grp$group <- factor(grp$group, c("Periop. chemotherapy", "No periop. chemotherapy"))
pd_int <- max(mean(p$delta[, 2] > 0), mean(p$delta[, 2] < 0))
thm <- theme_minimal(base_size = 13) +
  theme(plot.title = element_text(face = "bold", size = 12), plot.subtitle = element_text(colour = "grey45", size = 9),
        panel.grid.minor = element_blank(), legend.position = "bottom", legend.title = element_blank(), legend.text = element_text(size = 10))
panel <- function(k, line, fill, ylab, title) {
  ggplot(abc[abc$panel == k, ], aes(ki67)) + geom_ribbon(aes(ymin = lo, ymax = hi), fill = fill, alpha = 0.2) +
    geom_line(aes(y = med), colour = line, linewidth = 1.2) + labs(x = "Ki-67 (%)", y = ylab, title = title) + thm
}
pA <- panel("theta", "#1E40AF", "#3B82F6", expression(theta), "Clonogenic load")
pB <- panel("cure", "#047857", "#10B981", "P(never recur) %", "Cure rate") + ylim(0, NA)
pC <- panel("time", "#B45309", "#F59E0B", "Median failure time", "Relapse speed")
pD <- ggplot(de, aes(ki67, med)) + geom_hline(yintercept = 1, linetype = "dashed", colour = "grey60", linewidth = 0.4) +
  geom_linerange(aes(ymin = lo, ymax = hi), colour = "#047857", linewidth = 0.6, alpha = 0.55) +
  geom_point(colour = "#047857", fill = "#047857", shape = 21, size = 2.6, stroke = 0.3) + scale_y_log10() +
  labs(x = "Ki-67 (%)", y = "Cure Odds Ratio", title = "Cure OR: periop. chemotherapy vs none",
       subtitle = sprintf("95%% credible interval (log scale); chemotherapy started at %.1f months; interaction pd = %.1f%%", s_med, 100 * pd_int)) + thm
pE <- ggplot(grp, aes(ki67, med, colour = group)) +
  geom_linerange(aes(ymin = lo, ymax = hi), linewidth = 0.5, alpha = 0.55, position = position_dodge(width = 0.96)) +
  geom_point(size = 2.4, position = position_dodge(width = 0.96)) +
  scale_colour_manual(values = c("Periop. chemotherapy" = "#1E40AF", "No periop. chemotherapy" = "#B45309")) +
  ylim(0, NA) + labs(x = "Ki-67 (%)", y = "P(never recur) %", title = "Cure fraction by group") + thm
g <- (pA | pB | pC) / (pD | pE) + plot_annotation(tag_levels = "A") & theme(plot.tag = element_text(face = "bold", size = 14))
ggsave("results/figure_2.pdf", g, width = 14, height = 9)
