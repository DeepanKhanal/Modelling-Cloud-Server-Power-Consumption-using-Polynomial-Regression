###############################################################################
#  ST7089CEM  Introduction to Statistical Methods for Data Science
#  Coursework: Modelling Cloud Server Power Consumption using Polynomial
#              Regression
## ---------------------------------------------------------------------------
## 0. Set-up
## ---------------------------------------------------------------------------
rm(list = ls())
set.seed(7089)                      # makes the split and ABC reproducible
library(MASS)                       # kde2d() for the joint posterior contour

dir.create("figures", showWarnings = FALSE)
dir.create("outputs", showWarnings = FALSE)

# Helper: locate a data file under either naming convention
find_file <- function(candidates) {
  hit <- candidates[file.exists(candidates)]
  if (length(hit) == 0) stop("Data file not found. Looked for: ",
                             paste(candidates, collapse = ", "))
  hit[1]
}
f_x    <- find_file(c("X.csv", "x.csv", "x_-_Sheet1.csv"))
f_y    <- find_file(c("y.csv", "Y.csv", "y_-_Sheet1.csv"))
f_time <- find_file(c("time.csv", "time_-_Sheet1.csv"))

# Helper: open a PNG device with consistent settings
open_png <- function(name, w = 9, h = 6) {
  png(file.path("figures", name), width = w, height = h, units = "in",
      res = 200, type = "cairo")
}

## ---------------------------------------------------------------------------
## 1. Data import and quality checks
## ---------------------------------------------------------------------------
X_raw    <- read.csv(f_x)
y_raw    <- read.csv(f_y)
time_raw <- read.csv(f_time)

cat("Columns in X file:", names(X_raw), "\n")
cat("Columns in y file:", names(y_raw), "\n")
cat("Columns in time file:", names(time_raw), "\n")

x1 <- X_raw[[1]]      # CPU utilisation
x2 <- X_raw[[2]]      # memory usage
x3 <- X_raw[[3]]      # disk I/O
x4 <- X_raw[[4]]      # network traffic
y  <- y_raw[[1]]      # total power consumption
tt <- time_raw[[1]]   # sampling time (s)

signals <- data.frame(x1 = x1, x2 = x2, x3 = x3, x4 = x4, y = y)
sig_lab <- c(x1 = "x1: CPU utilisation", x2 = "x2: Memory usage",
             x3 = "x3: Disk I/O", x4 = "x4: Network traffic",
             y  = "y: Power consumption")

# Data-quality checks
n <- nrow(signals)
stopifnot(n == length(tt))                       # all files aligned
stopifnot(all(sapply(signals, is.numeric)))      # numeric types
dq <- data.frame(
  variable   = c(names(signals), "time"),
  n          = c(sapply(signals, length), length(tt)),
  n_missing  = c(sapply(signals, function(v) sum(is.na(v))), sum(is.na(tt))),
  n_nonfinite= c(sapply(signals, function(v) sum(!is.finite(v))),
                 sum(!is.finite(tt))),
  n_negative = c(sapply(signals, function(v) sum(v < 0)), sum(tt < 0))
)
print(dq)
write.csv(dq, "outputs/data_quality.csv", row.names = FALSE)

dt <- diff(tt)
cat(sprintf("Time: %.1f to %.1f s, step min %.3f / max %.3f s\n",
            min(tt), max(tt), min(dt), max(dt)))
cat("Duplicated time stamps:", sum(duplicated(tt)), "\n")

###############################################################################
## TASK 1 - PRELIMINARY DATA ANALYSIS
###############################################################################

## 1.0 Descriptive statistics ----------------------------------------------
skewness <- function(v) { m <- mean(v); mean((v - m)^3) / mean((v - m)^2)^1.5 }
kurtosis <- function(v) { m <- mean(v); mean((v - m)^4) / mean((v - m)^2)^2 }

desc <- data.frame(
  variable = names(signals),
  mean   = sapply(signals, mean),
  sd     = sapply(signals, sd),
  min    = sapply(signals, min),
  q25    = sapply(signals, quantile, 0.25),
  median = sapply(signals, median),
  q75    = sapply(signals, quantile, 0.75),
  max    = sapply(signals, max),
  skew   = sapply(signals, skewness),
  kurt   = sapply(signals, kurtosis)
)
print(desc, digits = 4)
write.csv(desc, "outputs/task1_descriptive_stats.csv", row.names = FALSE)

## 1.1 Time-series plots ----------------------------------------------------
cols <- c(x1 = "#1f77b4", x2 = "#2ca02c", x3 = "#ff7f0e",
          x4 = "#9467bd", y = "#d62728")

# (a) Full record (1000 s) for every signal
open_png("fig01_timeseries_full.png", 10, 11)
par(mfrow = c(5, 1), mar = c(3.5, 4.5, 2, 1), mgp = c(2.3, 0.7, 0))
for (v in names(signals)) {
  plot(tt, signals[[v]], type = "l", lwd = 0.3, col = cols[v],
       xlab = "Time (s)", ylab = v, main = sig_lab[v])
  # 10-second (100-sample) moving average to expose any trend/drift
  ma <- stats::filter(signals[[v]], rep(1 / 100, 100), sides = 2)
  lines(tt, ma, col = "black", lwd = 1.2)
}
dev.off()

# (b) Zoom on the first 20 s (200 samples) to show sample-to-sample behaviour
idx_zoom <- tt <= 20
open_png("fig02_timeseries_zoom.png", 10, 11)
par(mfrow = c(5, 1), mar = c(3.5, 4.5, 2, 1), mgp = c(2.3, 0.7, 0))
for (v in names(signals)) {
  plot(tt[idx_zoom], signals[[v]][idx_zoom], type = "o", pch = 16, cex = 0.4,
       lwd = 0.6, col = cols[v], xlab = "Time (s)", ylab = v,
       main = paste(sig_lab[v], "(first 20 s)"))
}
dev.off()

# (c) Autocorrelation - are consecutive samples (approximately) independent?
acf_lag1 <- sapply(signals, function(v) acf(v, lag.max = 1, plot = FALSE)$acf[2])
print(round(acf_lag1, 4))
# How many of the first 50 lags exceed the approximate 95% white-noise band?
acf_band <- 1.96 / sqrt(n)
acf_summary <- data.frame(
  signal = names(signals), lag1 = acf_lag1,
  n_outside_band = sapply(signals, function(v)
    sum(abs(acf(v, lag.max = 50, plot = FALSE)$acf[-1]) > acf_band)),
  max_abs_acf = sapply(signals, function(v)
    max(abs(acf(v, lag.max = 50, plot = FALSE)$acf[-1]))),
  band = acf_band)
print(acf_summary)
write.csv(acf_summary, "outputs/task1_acf_summary.csv", row.names = FALSE)
open_png("fig03_acf.png", 10, 6)
par(mfrow = c(2, 3), mar = c(4, 4, 3, 1))
for (v in names(signals)) acf(signals[[v]], lag.max = 50, main = sig_lab[v])
dev.off()

## 1.2 Distribution of each signal -----------------------------------------
open_png("fig04_histograms.png", 11, 7)
par(mfrow = c(2, 3), mar = c(4.2, 4.2, 3, 1))
for (v in names(signals)) {
  vv <- signals[[v]]
  h  <- hist(vv, breaks = 60, plot = FALSE)
  ytop <- max(h$density, density(vv)$y, dnorm(mean(vv), mean(vv), sd(vv)))
  hist(vv, breaks = 60, freq = FALSE, col = adjustcolor(cols[v], 0.35),
       border = "white", main = sig_lab[v], xlab = v, ylim = c(0, ytop * 1.05))
  lines(density(vv), lwd = 2, col = cols[v])                    # KDE
  curve(dnorm(x, mean(vv), sd(vv)), add = TRUE, lty = 2, lwd = 1.5)  # normal ref.
  rug(quantile(vv, c(0.25, 0.5, 0.75)), lwd = 2)
}
plot.new()
legend("center", c("Kernel density", "Normal with same mean/sd",
                   "Quartiles (rug)"),
       lty = c(1, 2, 1), lwd = c(2, 1.5, 2), bty = "n", cex = 1.2)
dev.off()

# y on a log scale - the raw y is dominated by a long right tail
open_png("fig05_y_log_hist.png", 8, 4.5)
par(mfrow = c(1, 2), mar = c(4.2, 4.2, 3, 1))
hist(log(y), breaks = 60, freq = FALSE, col = adjustcolor(cols["y"], 0.35),
     border = "white", main = "log(y)", xlab = "log(power)")
lines(density(log(y)), lwd = 2, col = cols["y"])
boxplot(signals[, 1:4], col = adjustcolor(cols[1:4], 0.4),
        main = "Inputs x1-x4 (box plots)", ylab = "value")
dev.off()

## 1.3 Correlation and scatter plots ----------------------------------------
cor_p <- cor(signals, method = "pearson")
cor_s <- cor(signals, method = "spearman")
print(round(cor_p, 3)); print(round(cor_s, 3))
write.csv(round(cor_p, 4), "outputs/task1_cor_pearson.csv")
write.csv(round(cor_s, 4), "outputs/task1_cor_spearman.csv")

# Significance of each input-output Pearson correlation
cor_tests <- do.call(rbind, lapply(c("x1", "x2", "x3", "x4"), function(v) {
  ct <- cor.test(signals[[v]], y)
  data.frame(input = v, r = unname(ct$estimate), t = unname(ct$statistic),
             df = unname(ct$parameter), p_value = ct$p.value,
             ci_low = ct$conf.int[1], ci_high = ct$conf.int[2],
             spearman = cor(signals[[v]], y, method = "spearman"))
}))
print(cor_tests)
write.csv(cor_tests, "outputs/task1_cor_tests.csv", row.names = FALSE)

# Heat-map of the Pearson correlation matrix
open_png("fig06_cor_heatmap.png", 7, 6)
par(mar = c(5, 5, 3, 6))
pal <- colorRampPalette(c("#2166ac", "white", "#b2182b"))(101)
image(1:5, 1:5, cor_p[, 5:1], col = pal, zlim = c(-1, 1), axes = FALSE,
      xlab = "", ylab = "", main = "Pearson correlation matrix")
axis(1, 1:5, colnames(cor_p)); axis(2, 1:5, rev(colnames(cor_p)), las = 1)
for (i in 1:5) for (j in 1:5)
  text(i, 6 - j, sprintf("%.2f", cor_p[i, j]), cex = 1.1)
box()
dev.off()

# Input-output scatter plots with a LOESS smoother
open_png("fig07_scatter_xy.png", 11, 8)
par(mfrow = c(2, 2), mar = c(4.2, 4.5, 3, 1))
for (v in c("x1", "x2", "x3", "x4")) {
  plot(signals[[v]], y, pch = 16, cex = 0.3, col = adjustcolor(cols[v], 0.25),
       xlab = sig_lab[v], ylab = "y: Power consumption",
       main = sprintf("y vs %s  (Pearson r = %.3f, Spearman = %.3f)",
                      v, cor_p[v, "y"], cor_s[v, "y"]))
  lo <- lowess(signals[[v]], y, f = 0.2)
  lines(lo, lwd = 2.5)
}
dev.off()

# Same scatter on log(y) to reveal structure hidden by the y tail
open_png("fig08_scatter_logy.png", 11, 8)
par(mfrow = c(2, 2), mar = c(4.2, 4.5, 3, 1))
for (v in c("x1", "x2", "x3", "x4")) {
  plot(signals[[v]], log(y), pch = 16, cex = 0.3,
       col = adjustcolor(cols[v], 0.25), xlab = sig_lab[v], ylab = "log(y)",
       main = paste("log(y) vs", v))
  lines(lowess(signals[[v]], log(y), f = 0.2), lwd = 2.5)
}
dev.off()

# Pairwise scatter among inputs (sub-sample of 2000 points for legibility)
sub <- sample(n, 2000)
open_png("fig09_pairs_inputs.png", 8, 8)
pairs(signals[sub, 1:4], pch = 16, cex = 0.3, col = adjustcolor("grey20", 0.4),
      labels = sig_lab[1:4], main = "Pairwise scatter of inputs (2000-point sample)")
dev.off()

## 1.4 Baseline linear model (first-order, all inputs) -----------------------
lin_fit <- lm(y ~ x1 + x2 + x3 + x4, data = signals)
print(summary(lin_fit))
lin_tab <- data.frame(term = names(coef(lin_fit)),
                      summary(lin_fit)$coefficients, row.names = NULL)
names(lin_tab) <- c("term", "estimate", "std_error", "t_value", "p_value")
write.csv(lin_tab, "outputs/task1_linear_model.csv", row.names = FALSE)
write.csv(data.frame(r_squared = summary(lin_fit)$r.squared,
                     adj_r_squared = summary(lin_fit)$adj.r.squared,
                     sigma = summary(lin_fit)$sigma),
          "outputs/task1_linear_model_fit.csv", row.names = FALSE)

open_png("fig10_linear_model_diag.png", 10, 4.5)
par(mfrow = c(1, 2), mar = c(4.2, 4.5, 3, 1))
plot(fitted(lin_fit), resid(lin_fit), pch = 16, cex = 0.3,
     col = adjustcolor("grey30", 0.3), xlab = "Fitted y", ylab = "Residual",
     main = "Linear model: residuals vs fitted")
abline(h = 0, col = "red", lty = 2)
lines(lowess(fitted(lin_fit), resid(lin_fit)), col = "blue", lwd = 2)
qqnorm(resid(lin_fit), pch = 16, cex = 0.3, main = "Linear model: normal Q-Q")
qqline(resid(lin_fit), col = "red", lwd = 2)
dev.off()

###############################################################################
## TASK 2 - REGRESSION
###############################################################################

## Design matrices for the five candidate models ---------------------------
## Column order follows the theta numbering in the brief; theta_bias is last.
## (Model 4 in the brief labels its third coefficient theta_4; it is the
##  model's third parameter and is called theta_3 here.)
make_design <- function(x1, x2, x3, x4) {
  ones <- rep(1, length(x1))
  list(
    Model1 = cbind(x4 = x4, x1_2 = x1^2, x1_3 = x1^3, x2_4 = x2^4, x1_4 = x1^4,
                   bias = ones),
    Model2 = cbind(x4 = x4, x1_3 = x1^3, x3_4 = x3^4, bias = ones),
    Model3 = cbind(x3_3 = x3^3, x3_4 = x3^4, bias = ones),
    Model4 = cbind(x2 = x2, x1_3 = x1^3, x3_4 = x3^4, bias = ones),
    Model5 = cbind(x4 = x4, x1_2 = x1^2, x1_3 = x1^3, x3_4 = x3^4, bias = ones)
  )
}
model_formula <- c(
  Model1 = "y = t1*x4 + t2*x1^2 + t3*x1^3 + t4*x2^4 + t5*x1^4 + t_bias",
  Model2 = "y = t1*x4 + t2*x1^3 + t3*x3^4 + t_bias",
  Model3 = "y = t1*x3^3 + t2*x3^4 + t_bias",
  Model4 = "y = t1*x2 + t2*x1^3 + t3*x3^4 + t_bias",
  Model5 = "y = t1*x4 + t2*x1^2 + t3*x1^3 + t4*x3^4 + t_bias")
designs <- make_design(x1, x2, x3, x4)

## Task 2.1 - least-squares estimates  theta_hat = (X'X)^-1 X'y -------------
ls_fit <- function(Xm, yv) {
  XtX <- crossprod(Xm)                       # X'X
  theta <- solve(XtX, crossprod(Xm, yv))     # solves (X'X) theta = X'y
  drop(theta)
}
thetas <- lapply(designs, ls_fit, yv = y)

# Cross-check against lm() (QR decomposition) - must agree to rounding
for (m in names(designs)) {
  Xm <- designs[[m]]
  chk <- coef(lm(y ~ Xm - 1))
  stopifnot(max(abs(chk - thetas[[m]])) < 1e-6 * max(1, max(abs(chk))))
}
cat("Normal-equation estimates agree with lm() for all models.\n")

## Tasks 2.2-2.4 - RSS, log-likelihood, AIC, BIC ----------------------------
model_stats <- function(Xm, yv, theta) {
  nn    <- length(yv)
  k     <- ncol(Xm)                          # number of estimated thetas
  res   <- drop(yv - Xm %*% theta)
  rss   <- sum(res^2)                        # Task 2.2
  s2    <- rss / (nn - 1)                    # sigma^2 as defined in the brief
  loglik<- -nn / 2 * log(2 * pi) - nn / 2 * log(s2) - rss / (2 * s2)  # Task 2.3
  aic   <- 2 * k - 2 * loglik                # Task 2.4
  bic   <- k * log(nn) - 2 * loglik
  r2    <- 1 - rss / sum((yv - mean(yv))^2)
  list(res = res, k = k, rss = rss, sigma2 = s2, loglik = loglik,
       aic = aic, bic = bic, r2 = r2, rmse = sqrt(rss / nn))
}
stats_all <- mapply(model_stats, designs, MoreArgs = list(yv = y),
                    theta = thetas, SIMPLIFY = FALSE)

# Standard errors / t-values for every parameter (for interpretation)
param_table <- do.call(rbind, lapply(names(designs), function(m) {
  Xm <- designs[[m]]; st <- stats_all[[m]]
  s2_unb <- st$rss / (n - st$k)
  se <- sqrt(diag(s2_unb * solve(crossprod(Xm))))
  tv <- thetas[[m]] / se
  data.frame(model = m,
             parameter = c(paste0("theta", seq_len(st$k - 1)), "theta_bias"),
             regressor = colnames(Xm), estimate = thetas[[m]], std_error = se,
             t_value = tv, p_value = 2 * pt(-abs(tv), df = n - st$k),
             row.names = NULL)
}))
print(param_table, digits = 5)
write.csv(param_table, "outputs/task2_1_parameters.csv", row.names = FALSE)

summary_table <- data.frame(
  model   = names(designs),
  formula = model_formula,
  k       = sapply(stats_all, `[[`, "k"),
  RSS     = sapply(stats_all, `[[`, "rss"),
  sigma2  = sapply(stats_all, `[[`, "sigma2"),
  logLik  = sapply(stats_all, `[[`, "loglik"),
  AIC     = sapply(stats_all, `[[`, "aic"),
  BIC     = sapply(stats_all, `[[`, "bic"),
  R2      = sapply(stats_all, `[[`, "r2"),
  RMSE    = sapply(stats_all, `[[`, "rmse"),
  row.names = NULL)
summary_table$dAIC <- summary_table$AIC - min(summary_table$AIC)
summary_table$dBIC <- summary_table$BIC - min(summary_table$BIC)
summary_table$AIC_weight <- exp(-summary_table$dAIC / 2) /
                            sum(exp(-summary_table$dAIC / 2))
print(summary_table, digits = 6)
write.csv(summary_table, "outputs/task2_2to4_model_comparison.csv",
          row.names = FALSE)

# Bar charts of AIC and BIC
open_png("fig11_aic_bic.png", 10, 4.5)
par(mfrow = c(1, 2), mar = c(4.5, 5.5, 3, 1))
barplot(summary_table$AIC, names.arg = summary_table$model, col = "#4c72b0",
        main = "AIC by candidate model", ylab = "", las = 1, cex.names = 0.9)
mtext("AIC", side = 2, line = 4.3)
barplot(summary_table$BIC, names.arg = summary_table$model, col = "#dd8452",
        main = "BIC by candidate model", ylab = "", las = 1, cex.names = 0.9)
mtext("BIC", side = 2, line = 4.3)
dev.off()

## Task 2.5 - residual distributions and Q-Q plots ---------------------------
jarque_bera <- function(r) {
  nn <- length(r); S <- skewness(r); K <- kurtosis(r)
  JB <- nn / 6 * (S^2 + (K - 3)^2 / 4)
  c(skew = S, kurt = K, JB = JB, p = pchisq(JB, df = 2, lower.tail = FALSE))
}
idx_sw <- sample(n, 5000)          # shapiro.test() accepts at most 5000 values
resid_tests <- do.call(rbind, lapply(names(designs), function(m) {
  r  <- stats_all[[m]]$res
  jb <- jarque_bera(r)
  sw <- shapiro.test(r[idx_sw])
  data.frame(model = m, mean = mean(r), sd = sd(r), skewness = jb["skew"],
             kurtosis = jb["kurt"], JB = jb["JB"], JB_p = jb["p"],
             SW_W = unname(sw$statistic), SW_p = sw$p.value, row.names = NULL)
}))
print(resid_tests, digits = 4)
write.csv(resid_tests, "outputs/task2_5_residual_tests.csv", row.names = FALSE)

open_png("fig12_residual_hist.png", 12, 6.5)
par(mfrow = c(2, 3), mar = c(4.2, 4.2, 3, 1))
for (m in names(designs)) {
  r <- stats_all[[m]]$res
  hist(r, breaks = 80, freq = FALSE, col = "grey80", border = "white",
       main = paste(m, "residuals"), xlab = "Residual")
  curve(dnorm(x, mean(r), sd(r)), add = TRUE, col = "red", lwd = 2)
}
plot.new()
legend("center", c("Residual histogram", "Normal density (same mean, sd)"),
       fill = c("grey80", NA), border = c("grey50", NA),
       lty = c(NA, 1), col = c(NA, "red"), lwd = c(NA, 2), bty = "n", cex = 1.2)
dev.off()

open_png("fig13_residual_qq.png", 12, 6.5)
par(mfrow = c(2, 3), mar = c(4.2, 4.2, 3, 1))
for (m in names(designs)) {
  r <- stats_all[[m]]$res
  qqnorm(r, pch = 16, cex = 0.3, main = paste(m, "normal Q-Q"))
  qqline(r, col = "red", lwd = 2)
}
plot.new()
dev.off()

# Residuals vs fitted values and vs time for the two best models
open_png("fig14_resid_vs_fitted.png", 11, 7.5)
par(mfrow = c(2, 2), mar = c(4.2, 4.5, 3, 1))
for (m in c("Model2", "Model5")) {
  fit <- drop(designs[[m]] %*% thetas[[m]]); r <- stats_all[[m]]$res
  plot(fit, r, pch = 16, cex = 0.3, col = adjustcolor("grey20", 0.3),
       log = "x", xlab = "Fitted y (log axis)", ylab = "Residual",
       main = paste(m, ": residuals vs fitted"))
  abline(h = 0, col = "red", lty = 2)
  plot(tt, r, type = "l", lwd = 0.3, col = "grey30", xlab = "Time (s)",
       ylab = "Residual", main = paste(m, ": residuals vs time"))
  abline(h = 0, col = "red", lty = 2)
}
dev.off()
resid_acf1 <- sapply(names(designs), function(m)
  acf(stats_all[[m]]$res, lag.max = 1, plot = FALSE)$acf[2])
print(round(resid_acf1, 4))
write.csv(data.frame(model = names(resid_acf1), lag1_acf = resid_acf1),
          "outputs/task2_5_residual_acf.csv", row.names = FALSE)

## Task 2.6 - supporting test: Model 2 is nested in Model 5 -------------------
## Model 5 = Model 2 + theta*x1^2.  Partial F-test and likelihood-ratio test.
rss2 <- stats_all$Model2$rss; rss5 <- stats_all$Model5$rss
df5  <- n - stats_all$Model5$k
F_nested <- ((rss2 - rss5) / 1) / (rss5 / df5)
p_F      <- pf(F_nested, 1, df5, lower.tail = FALSE)
LR       <- n * log(rss2 / rss5)                  # ML-based LR statistic
p_LR     <- pchisq(LR, df = 1, lower.tail = FALSE)
nested <- data.frame(F = F_nested, df1 = 1, df2 = df5, p_F = p_F,
                     LR = LR, p_LR = p_LR)
print(nested)
write.csv(nested, "outputs/task2_6_nested_test.csv", row.names = FALSE)

# Polynomial terms of the same input are strongly collinear, which inflates
# the standard errors of Model 5 (and Model 1) relative to Model 2.
poly_cor <- data.frame(
  pair = c("x1^2 vs x1^3", "x1^3 vs x1^4", "x1^2 vs x1^4"),
  r = c(cor(x1^2, x1^3), cor(x1^3, x1^4), cor(x1^2, x1^4)))
se_m2 <- param_table$std_error[param_table$model == "Model2"]
se_m5 <- param_table$std_error[param_table$model == "Model5"]
vif_like <- data.frame(term = c("x4", "x1^3", "x3^4", "bias"),
                       se_model2 = se_m2, se_model5 = se_m5[c(1, 3, 4, 5)],
                       se_ratio = se_m5[c(1, 3, 4, 5)] / se_m2)
print(poly_cor); print(vif_like)
write.csv(poly_cor, "outputs/task2_6_poly_collinearity.csv", row.names = FALSE)
write.csv(vif_like, "outputs/task2_6_se_inflation.csv", row.names = FALSE)

best_model <- summary_table$model[which.min(summary_table$BIC)]

# Share of the variance of y explained by each fitted term of the selected
# model (inputs are independent, so the term variances are close to additive)
Xb_sel <- designs[[best_model]]; th_sel <- thetas[[best_model]]
term_var <- sapply(seq_len(ncol(Xb_sel) - 1), function(j) var(th_sel[j] * Xb_sel[, j]))
var_share <- data.frame(term = colnames(Xb_sel)[-ncol(Xb_sel)],
                        variance = term_var, share_of_var_y = term_var / var(y))
var_share <- rbind(var_share, data.frame(term = "residual",
                   variance = var(stats_all[[best_model]]$res),
                   share_of_var_y = var(stats_all[[best_model]]$res) / var(y)))
print(var_share)
write.csv(var_share, "outputs/task2_6_variance_shares.csv", row.names = FALSE)
cat("Selected model (min AIC and BIC):", best_model, "\n")
stopifnot(best_model == summary_table$model[which.min(summary_table$AIC)])

# Component-plus-residual (partial residual) plots for the selected model:
# each panel shows theta_j * term_j + residual against the original input,
# so the fitted shape of every term can be checked against the data.
Xb <- designs[[best_model]]; thb <- thetas[[best_model]]
rb <- stats_all[[best_model]]$res
term_src <- c(x4 = "x4", x1_3 = "x1", x3_4 = "x3", x1_2 = "x1", x2 = "x2",
              x3_3 = "x3", x2_4 = "x2", x1_4 = "x1")
nterm <- ncol(Xb) - 1
open_png("fig14b_partial_residuals.png", 4 * nterm, 4.3)
par(mfrow = c(1, nterm), mar = c(4.2, 4.5, 3, 1))
for (j in seq_len(nterm)) {
  src <- term_src[colnames(Xb)[j]]
  comp <- thb[j] * Xb[, j]
  o <- order(signals[[src]])
  plot(signals[[src]], comp + rb, pch = 16, cex = 0.25,
       col = adjustcolor("grey30", 0.3), xlab = sig_lab[src],
       ylab = "Component + residual",
       main = sprintf("%.4g * %s", thb[j], colnames(Xb)[j]))
  lines(signals[[src]][o], comp[o], col = "red", lwd = 2)
}
dev.off()

## Task 2.7 - 70/30 train/test split, prediction and 95% intervals ----------
idx_train <- sort(sample(n, round(0.7 * n)))
idx_test  <- setdiff(seq_len(n), idx_train)
X_all  <- designs[[best_model]]
X_tr   <- X_all[idx_train, , drop = FALSE]; y_tr <- y[idx_train]
X_te   <- X_all[idx_test,  , drop = FALSE]; y_te <- y[idx_test]

theta_tr <- ls_fit(X_tr, y_tr)                                  # 2.7(1)
n_tr <- length(y_tr); k_b <- ncol(X_tr)
rss_tr   <- sum((y_tr - X_tr %*% theta_tr)^2)
s2_tr    <- rss_tr / (n_tr - k_b)        # unbiased residual variance
XtX_inv  <- solve(crossprod(X_tr))
cov_theta<- s2_tr * XtX_inv              # Cov(theta_hat)

y_hat_te <- drop(X_te %*% theta_tr)                              # 2.7(2)
se_mean  <- sqrt(rowSums((X_te %*% XtX_inv) * X_te) * s2_tr)     # SE of E[y|x]
se_pred  <- sqrt(se_mean^2 + s2_tr)                              # SE of new y
t_crit   <- qt(0.975, df = n_tr - k_b)
ci_lo <- y_hat_te - t_crit * se_mean;  ci_hi <- y_hat_te + t_crit * se_mean
pi_lo <- y_hat_te - t_crit * se_pred;  pi_hi <- y_hat_te + t_crit * se_pred

res_te <- y_te - y_hat_te
test_metrics <- data.frame(
  n_train = n_tr, n_test = length(y_te),
  RMSE = sqrt(mean(res_te^2)), MAE = mean(abs(res_te)),
  R2 = 1 - sum(res_te^2) / sum((y_te - mean(y_te))^2),
  train_RMSE = sqrt(rss_tr / n_tr), sigma_hat = sqrt(s2_tr), t_crit = t_crit,
  mean_CI_halfwidth = mean(t_crit * se_mean),
  max_CI_halfwidth  = max(t_crit * se_mean),
  mean_PI_halfwidth = mean(t_crit * se_pred),
  coverage_CI = mean(y_te >= ci_lo & y_te <= ci_hi),
  coverage_PI = mean(y_te >= pi_lo & y_te <= pi_hi))
print(test_metrics)
write.csv(test_metrics, "outputs/task2_7_test_metrics.csv", row.names = FALSE)

param_tr <- data.frame(parameter = c(paste0("theta", 1:(k_b - 1)), "theta_bias"),
                       regressor = colnames(X_tr), estimate_train = theta_tr,
                       std_error = sqrt(diag(cov_theta)),
                       ci95_low = theta_tr - t_crit * sqrt(diag(cov_theta)),
                       ci95_high = theta_tr + t_crit * sqrt(diag(cov_theta)),
                       estimate_full = thetas[[best_model]], row.names = NULL)
print(param_tr, digits = 6)
write.csv(param_tr, "outputs/task2_7_train_parameters.csv", row.names = FALSE)

pred_out <- data.frame(index = idx_test, time = tt[idx_test], y_test = y_te,
                       y_hat = y_hat_te, ci_low = ci_lo, ci_high = ci_hi,
                       pi_low = pi_lo, pi_high = pi_hi)
write.csv(pred_out, "outputs/task2_7_test_predictions.csv", row.names = FALSE)

# Figure: predicted vs observed for all test points
open_png("fig15_pred_vs_obs.png", 10, 4.8)
par(mfrow = c(1, 2), mar = c(4.2, 4.5, 3, 1))
plot(y_hat_te, y_te, pch = 16, cex = 0.35, col = adjustcolor("#1f77b4", 0.4),
     xlab = "Predicted y", ylab = "Observed y (test)",
     main = sprintf("%s on test data (R2 = %.4f)", best_model, test_metrics$R2))
abline(0, 1, col = "red", lwd = 2)
hist(res_te, breaks = 60, freq = FALSE, col = "grey80", border = "white",
     main = "Test residuals", xlab = "Observed - predicted")
curve(dnorm(x, 0, sqrt(s2_tr)), add = TRUE, col = "red", lwd = 2)
dev.off()

# Figure: 60 consecutive test samples with 95% CI error bars, PI band and data
sel <- 1:60
xs  <- tt[idx_test][sel]
open_png("fig16_prediction_ci.png", 11, 8)
par(mfrow = c(2, 1), mar = c(4.2, 4.8, 3, 1))
plot(xs, y_te[sel], type = "n", ylim = range(pi_lo[sel], pi_hi[sel], y_te[sel]),
     xlab = "Time (s)", ylab = "Power consumption y",
     main = "Test set (first 60 samples): prediction, 95% CI and test data")
polygon(c(xs, rev(xs)), c(pi_lo[sel], rev(pi_hi[sel])),
        col = adjustcolor("orange", 0.25), border = NA)
lines(xs, y_hat_te[sel], col = "#1f77b4", lwd = 1.2)
suppressWarnings(  # bars shorter than the arrow head trigger a harmless warning
  arrows(xs, ci_lo[sel], xs, ci_hi[sel], angle = 90, code = 3, length = 0.02,
         col = "#1f77b4"))
points(xs, y_hat_te[sel], pch = 16, cex = 0.7, col = "#1f77b4")
points(xs, y_te[sel], pch = 4, cex = 0.9, col = "red", lwd = 1.5)
legend("topleft", c("Model prediction", "95% CI (error bars)",
                    "95% prediction interval", "Test data"),
       col = c("#1f77b4", "#1f77b4", adjustcolor("orange", 0.5), "red"),
       pch = c(16, NA, 15, 4), lty = c(1, 1, NA, NA), pt.cex = c(1, 1, 2, 1),
       bty = "n", cex = 0.85)
# Same samples, centred on the prediction so the interval widths are visible
plot(xs, y_te[sel] - y_hat_te[sel], type = "n",
     ylim = range(pi_lo[sel] - y_hat_te[sel], pi_hi[sel] - y_hat_te[sel],
                  y_te[sel] - y_hat_te[sel]),
     xlab = "Time (s)", ylab = "Deviation from prediction",
     main = "Same samples, centred on the prediction (y - y_hat)")
polygon(c(xs, rev(xs)), c(pi_lo[sel] - y_hat_te[sel], rev(pi_hi[sel] - y_hat_te[sel])),
        col = adjustcolor("orange", 0.25), border = NA)
arrows(xs, ci_lo[sel] - y_hat_te[sel], xs, ci_hi[sel] - y_hat_te[sel],
       angle = 90, code = 3, length = 0.03, col = "#1f77b4", lwd = 1.5)
abline(h = 0, col = "#1f77b4")
points(xs, y_te[sel] - y_hat_te[sel], pch = 4, col = "red", lwd = 1.5)
dev.off()

###############################################################################
## TASK 3 - APPROXIMATE BAYESIAN COMPUTATION (rejection ABC)
###############################################################################

## 3.1 Choose the two parameters with largest |theta_hat| (from Task 2.1) ----
theta_best <- thetas[[best_model]]
pnames     <- c(paste0("theta", 1:(length(theta_best) - 1)), "theta_bias")
names(theta_best) <- pnames
ord   <- order(abs(theta_best), decreasing = TRUE)
abc_i <- ord[1:2]                         # indices of the 2 free parameters
fix_i <- setdiff(seq_along(theta_best), abc_i)
cat("ABC parameters:", pnames[abc_i], "=", theta_best[abc_i], "\n")
cat("Fixed parameters:", pnames[fix_i], "=", theta_best[fix_i], "\n")

X_best <- designs[[best_model]]
se_best <- param_table$std_error[param_table$model == best_model]
rss_hat <- stats_all[[best_model]]$rss
s2_hat  <- rss_hat / (n - ncol(X_best))

## 3.2 Uniform priors centred on the LS estimates -----------------------------
## Half-width = 10 standard errors: wide relative to the likelihood (so the
## data, not the prior, determine the posterior) yet narrow enough that the
## rejection sampler accepts a usable number of draws.
prior_halfwidth <- 10 * se_best[abc_i]
prior_lo <- theta_best[abc_i] - prior_halfwidth
prior_hi <- theta_best[abc_i] + prior_halfwidth
prior_tab <- data.frame(parameter = pnames[abc_i],
                        regressor = colnames(X_best)[abc_i],
                        ls_estimate = theta_best[abc_i], std_error = se_best[abc_i],
                        prior_low = prior_lo, prior_high = prior_hi)
print(prior_tab, digits = 6)
write.csv(prior_tab, "outputs/task3_priors.csv", row.names = FALSE)

## 3.3 Rejection ABC ---------------------------------------------------------
## For each prior draw: simulate the model output y_sim = X theta (other
## parameters fixed), compute the distance d = RSS(y, y_sim) and accept if
## d <= epsilon.  epsilon = RSS_hat + chi2_{2,0.95} * sigma2_hat, i.e. the
## RSS increase that corresponds to a 95% joint region for 2 parameters.
N_abc <- 200000
draws <- cbind(runif(N_abc, prior_lo[1], prior_hi[1]),
               runif(N_abc, prior_lo[2], prior_hi[2]))
colnames(draws) <- pnames[abc_i]

y_fixed <- drop(X_best[, fix_i, drop = FALSE] %*% theta_best[fix_i])
X_free  <- X_best[, abc_i, drop = FALSE]

abc_rss <- numeric(N_abc)
chunk <- 1000
for (s in seq(1, N_abc, by = chunk)) {         # chunked to limit memory use
  j <- s:min(s + chunk - 1, N_abc)
  y_sim <- y_fixed + X_free %*% t(draws[j, , drop = FALSE])   # n x |j| matrix
  abc_rss[j] <- colSums((y - y_sim)^2)
}

eps_levels <- c(strict = qchisq(0.50, 2), main = qchisq(0.95, 2),
                loose = qchisq(0.99, 2))
eps <- rss_hat + eps_levels * s2_hat
accept <- abc_rss <= eps["main"]
post   <- draws[accept, , drop = FALSE]
cat(sprintf("Accepted %d of %d draws (%.2f%%)\n", sum(accept), N_abc,
            100 * mean(accept)))

## Analytical reference: under a flat prior and Gaussian noise the posterior
## of the two free parameters (others fixed) is normal with covariance
## sigma^2 (X_free' X_free)^-1, centred on the conditional LS estimate.
cond_cov <- s2_hat * solve(crossprod(X_free))
cond_mean <- drop(solve(crossprod(X_free), crossprod(X_free, y - y_fixed)))

post_summary <- data.frame(
  parameter   = pnames[abc_i],
  ls_estimate = theta_best[abc_i],
  post_mean   = colMeans(post),
  post_median = apply(post, 2, median),
  post_sd     = apply(post, 2, sd),
  cri95_low   = apply(post, 2, quantile, 0.025),
  cri95_high  = apply(post, 2, quantile, 0.975),
  analytic_sd_conditional = sqrt(diag(cond_cov)),
  ls_se_full_model = se_best[abc_i], row.names = NULL)
post_cor <- cor(post)[1, 2]
print(post_summary, digits = 6)
cat("Posterior correlation:", post_cor, "\n")
cat("Analytical conditional correlation:",
    cond_cov[1, 2] / sqrt(prod(diag(cond_cov))), "\n")
write.csv(post_summary, "outputs/task3_posterior_summary.csv", row.names = FALSE)

abc_meta <- data.frame(N_draws = N_abc, accepted = sum(accept),
                       acceptance_rate = mean(accept), rss_hat = rss_hat,
                       sigma2_hat = s2_hat, epsilon = eps["main"],
                       posterior_correlation = post_cor,
                       analytic_correlation = cond_cov[1, 2] /
                         sqrt(prod(diag(cond_cov))))
write.csv(abc_meta, "outputs/task3_abc_settings.csv", row.names = FALSE)

# Sensitivity of the posterior to the tolerance epsilon
sens <- do.call(rbind, lapply(names(eps), function(e) {
  a <- abc_rss <= eps[e]; p <- draws[a, , drop = FALSE]
  data.frame(tolerance = e, chi2_quantile = eps_levels[e],
             epsilon = eps[e], accepted = sum(a),
             mean_1 = mean(p[, 1]), sd_1 = sd(p[, 1]),
             mean_2 = mean(p[, 2]), sd_2 = sd(p[, 2]),
             correlation = cor(p)[1, 2], row.names = NULL)
}))
print(sens, digits = 6)
write.csv(sens, "outputs/task3_tolerance_sensitivity.csv", row.names = FALSE)

## 3.4 Plots: marginal and joint posteriors ---------------------------------
lab1 <- sprintf("%s (%s)", pnames[abc_i[1]], colnames(X_best)[abc_i[1]])
lab2 <- sprintf("%s (%s)", pnames[abc_i[2]], colnames(X_best)[abc_i[2]])

# Analytical Gaussian posterior (flat prior, other parameters fixed) for
# comparison with the ABC result
an_sd <- sqrt(diag(cond_cov))

open_png("fig17_abc_marginals.png", 11, 8.5)
par(mfrow = c(2, 2), mar = c(4.5, 4.5, 3, 1))
# Row 1: posterior against the full prior range
for (q in 1:2) {
  hist(post[, q], breaks = 40, freq = FALSE, col = c("#6baed6", "#74c476")[q],
       border = NA, xlab = c(lab1, lab2)[q], xlim = c(prior_lo[q], prior_hi[q]),
       main = paste(pnames[abc_i[q]], ": posterior vs uniform prior"))
  segments(prior_lo[q], 1 / (prior_hi[q] - prior_lo[q]),
           prior_hi[q], 1 / (prior_hi[q] - prior_lo[q]), lwd = 2, lty = 3)
  abline(v = theta_best[abc_i[q]], col = "red", lwd = 2, lty = 2)
  legend("topright", c("ABC posterior", "Uniform prior", "LS estimate"),
         pch = c(15, NA, NA), pt.cex = 2, lty = c(0, 3, 2), lwd = c(1, 2, 2),
         col = c(c("#6baed6", "#74c476")[q], "black", "red"),
         bty = "n", cex = 0.8)
}
# Row 2: zoom on the posterior
for (q in 1:2) {
  hist(post[, q], breaks = 40, freq = FALSE, col = c("#6baed6", "#74c476")[q],
       border = "white", xlab = c(lab1, lab2)[q],
       main = paste(pnames[abc_i[q]], ": posterior (zoomed)"),
       ylim = c(0, 1.25 * max(dnorm(0, 0, an_sd[q]),
                              hist(post[, q], breaks = 40, plot = FALSE)$density)))
  lines(density(post[, q]), lwd = 2)
  curve(dnorm(x, cond_mean[q], an_sd[q]), add = TRUE, lwd = 2, lty = 2,
        col = "darkorange")
  abline(v = theta_best[abc_i[q]], col = "red", lwd = 2, lty = 2)
  abline(v = quantile(post[, q], c(0.025, 0.975)), col = "purple", lty = 4)
  legend("topright", c("ABC posterior (KDE)", "Analytical Gaussian posterior",
                       "LS estimate", "95% credible interval"),
         lty = c(1, 2, 2, 4), col = c("black", "darkorange", "red", "purple"),
         lwd = c(2, 2, 2, 1), bty = "n", cex = 0.75)
}
dev.off()

# Theoretical acceptance region: {theta : RSS(theta) <= epsilon} is the ellipse
# (theta - theta_hat)' X'X (theta - theta_hat) <= epsilon - RSS_hat
ell <- function(level_rss) {
  ang <- seq(0, 2 * pi, length.out = 300)
  L <- t(chol(solve(crossprod(X_free)) * level_rss))
  pts <- t(L %*% rbind(cos(ang), sin(ang)))
  sweep(pts, 2, theta_best[abc_i], "+")
}
open_png("fig18_abc_joint.png", 11, 5.2)
par(mfrow = c(1, 2), mar = c(4.5, 4.5, 3, 1))
plot(draws[!accept, ][1:20000, ], pch = 16, cex = 0.2, col = "grey85",
     xlab = lab1, ylab = lab2, main = "Prior draws (grey) and accepted (blue)")
points(post, pch = 16, cex = 0.25, col = adjustcolor("#08519c", 0.5))
lines(ell(eps["main"] - rss_hat), col = "darkorange", lwd = 2)
points(theta_best[abc_i[1]], theta_best[abc_i[2]], pch = 4, cex = 2,
       lwd = 3, col = "red")
legend("topright", c("Accepted draws", "Theoretical acceptance ellipse",
                     "LS estimate"), pch = c(16, NA, 4), lty = c(NA, 1, NA),
       col = c("#08519c", "darkorange", "red"), lwd = c(NA, 2, 3), bty = "n",
       cex = 0.8)
bw <- c(bandwidth.nrd(post[, 1]), bandwidth.nrd(post[, 2])) * 1.3
kd <- kde2d(post[, 1], post[, 2], n = 120, h = bw)
image(kd, col = hcl.colors(30, "YlOrRd", rev = TRUE), xlab = lab1, ylab = lab2,
      main = "Joint posterior density (2-D KDE)")
contour(kd, add = TRUE, drawlabels = FALSE, col = "grey30", nlevels = 8)
points(theta_best[abc_i[1]], theta_best[abc_i[2]], pch = 4, cex = 2, lwd = 3,
       col = "blue")
legend("topright", "LS estimate", pch = 4, col = "blue", lty = 0, pt.lwd = 3, bty = "n")
dev.off()

## ---------------------------------------------------------------------------
## Session information (for reproducibility)
## ---------------------------------------------------------------------------
writeLines(capture.output(sessionInfo()), "outputs/session_info.txt")
cat("\nAnalysis complete. Figures in ./figures, tables in ./outputs\n")
