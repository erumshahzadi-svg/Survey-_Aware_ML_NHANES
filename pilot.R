# =============================================================================
# pilot.R  --  Survey-weighted and class-weighted logistic training, one scenario
# Prevalence 3%, n = 2,000, stratified two-stage cluster design. Base R only.
#
# HOW TO USE (RStudio):
#   source("pilot.R")                       # loads the functions, runs nothing
#   pop <- make_population(SEED)            # STEP 1  (about 10 seconds)
#   demo_one_replicate(pop)                 # STEPS 2-5 once, with printed output
#   res <- run_all(reps = 500)              # the whole simulation (about 10 min)
# HOW TO USE (terminal):   Rscript pilot.R 500
# =============================================================================

## ---- STEP 0. Settings ------------------------------------------------------
H <- 50; J <- 40; M <- 100        # 50 strata, 40 PSUs per stratum, 100 people per PSU
P <- 5                            # number of predictors
DRAWS <- 2; PER <- 20             # per stratum: 2 PSU draws, 20 people per draw -> n = 2000
LAMBDAS <- c(1, 2, 5, 10, 20)     # class weights (lambda = 1 is Method B)
ALPHA <- 0.05                     # FPR limit used to tune the threshold
PREV <- 0.03; TARGET_AUC <- 0.75  # population prevalence and signal strength
BASE <- c(1, 0.6, 0.4, 0.2, 0)    # relative size of the 5 true coefficients
SEED <- 2026

# AUROC from ranks (Mann-Whitney); avoids needing extra packages
auc <- function(y, p) {
  r <- rank(p); n1 <- sum(y == 1); n0 <- sum(y == 0)
  (sum(r[y == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}

## ---- STEP 1. Build the finite population (the "truth") ----------------------
# We create 200,000 people once. Because we know everyone's outcome, we can
# later compute the TRUE sensitivity, FPR and calibration of any fitted model.
make_population <- function(seed = SEED) {
  set.seed(seed)
  n_psu <- H * J; N <- n_psu * M
  psu   <- rep(seq_len(n_psu), each = M)            # PSU label of each person
  delta <- rnorm(n_psu, 0, 0.3)                     # PSU shift in x1 (clustering in predictors)
  u     <- rnorm(n_psu, 0, 0.5)                     # PSU random effect in the outcome
  X <- matrix(rnorm(N * P), N, P)
  X[, 2] <- 0.4 * X[, 1] + sqrt(1 - 0.16) * X[, 2]  # mild correlation x1-x2
  X[, 1] <- X[, 1] + delta[psu]
  v <- runif(N)                                     # fixed uniforms -> fixed outcomes
  eta <- function(s) as.vector(X %*% (s * BASE)) + u[psu]
  # intercept b0 so that mean probability = PREV
  solve_b0 <- function(e) uniroot(function(b) mean(plogis(b + e)) - PREV, c(-20, 5))$root
  # scale s so that the true probabilities have AUROC = TARGET_AUC
  auc_gap <- function(s) {
    e <- eta(s); p <- plogis(solve_b0(e) + e); auc(as.integer(v < p), p) - TARGET_AUC
  }
  s  <- uniroot(auc_gap, c(0.2, 3), tol = 1e-3)$root
  e  <- eta(s); b0 <- solve_b0(e); p_true <- plogis(b0 + e)
  y  <- as.integer(v < p_true)                      # the population outcome
  # size measure per PSU: bigger when mean x1 is bigger -> selection depends on a predictor
  size <- exp(0.8 * colMeans(matrix(X[, 1], nrow = M)) + 0.3 * rnorm(n_psu))
  cat(sprintf("Population: N=%d, prevalence=%.4f, AUROC of true probabilities=%.3f\n",
              N, mean(y), auc(y, p_true)))
  list(X = X, y = y, size = size, auc_true = auc(y, p_true))
}

## ---- STEP 2. Draw one survey sample ----------------------------------------
# In each stratum: draw 2 PSUs with probability proportional to size (with
# replacement), then 20 people by simple random sampling inside each draw.
# Weight of a sampled person = 1 / (probability of being selected)
#                            = 1 / (DRAWS * p_PSU * PER / M).
draw_sample <- function(pop) {
  idx <- vector("list", H * DRAWS); w <- vector("list", H * DRAWS); k <- 0
  for (h in seq_len(H)) {
    js <- ((h - 1) * J + 1):(h * J)                   # PSUs of stratum h
    pr <- pop$size[js] / sum(pop$size[js])            # draw probabilities
    for (d in sample.int(J, DRAWS, replace = TRUE, prob = pr)) {
      k <- k + 1
      idx[[k]] <- (js[d] - 1) * M + sample.int(M, PER)   # rows of the 20 sampled people
      w[[k]]   <- rep(1 / (DRAWS * pr[d] * PER / M), PER)
    }
  }
  list(idx = unlist(idx), w = unlist(w))
}

## ---- STEP 3. Fit the models (A, B, C) ---------------------------------------
# A: unweighted logistic regression.
# B: survey-weighted (weights w, rescaled to mean 1).
# C(lambda): weights w * lambda^Y -- each positive case counts lambda times more.
fit_logit <- function(X, y, sw) {
  f <- glm.fit(cbind(1, X), y, weights = sw, family = quasibinomial())
  list(b0 = f$coefficients[1], b = f$coefficients[-1], converged = f$converged)
}
predict_prob <- function(fit, X) plogis(fit$b0 + as.vector(X %*% fit$b))

## ---- STEP 4. Helper functions for evaluation ---------------------------------
# Calibration: regress the TRUE outcome on logit(predicted prob).
# Perfect calibration = intercept 0, slope 1. Returns c(intercept, slope).
calib <- function(y, p) {
  z <- qlogis(pmin(pmax(p, 1e-12), 1 - 1e-12)); A <- cbind(1, z); b <- c(0, 0)
  for (i in 1:30) {                                   # Newton-Raphson (IRLS)
    mu <- plogis(drop(A %*% b)); W <- mu * (1 - mu)
    step <- drop(solve(crossprod(A, A * W), crossprod(A, y - mu))); b <- b + step
    if (max(abs(step)) < 1e-9) break
  }
  b
}
# Tuned threshold: the smallest cut-off such that the WEIGHTED false-positive
# rate among sampled negatives is <= ALPHA. Predict "positive" if score > t.
tuned_threshold <- function(score, y, w) {
  sn <- score[y == 0]; wn <- w[y == 0]; o <- order(sn)
  cw <- cumsum(wn[o]) / sum(wn)
  sn[o][which(cw >= 1 - ALPHA)[1]]
}

## ---- STEP 5. Evaluate one fitted model against the population truth ------------
evaluate <- function(fit, pop, smp, ys, nm, lam) {
  p <- predict_prob(fit, pop$X); y <- pop$y; pos <- y == 1; neg <- !pos
  a_b <- unname(calib(y, p))
  t_or <- quantile(p[neg], 1 - ALPHA)                 # oracle threshold: FPR = 5% exactly
  t_tu <- tuned_threshold(p[smp$idx], ys, smp$w)      # threshold a real analyst could pick
  data.frame(method = nm, lam = lam, nonconv = as.integer(!fit$converged),
             sens_t05 = mean(p[pos] >= 0.5), fpr_t05 = mean(p[neg] >= 0.5),
             auc = auc(y, p), calib_int = a_b[1], calib_slope = a_b[2],
             sens_oracle = mean(p[pos] > t_or),
             sens_tuned = mean(p[pos] > t_tu), fpr_tuned = mean(p[neg] > t_tu))
}

one_replicate <- function(pop, rep_id, smp = NULL) {
  if (is.null(smp)) smp <- draw_sample(pop)
  Xs <- pop$X[smp$idx, ]; ys <- pop$y[smp$idx]
  wn <- smp$w / mean(smp$w)                           # rescale weights to mean 1
  out <- list(evaluate(fit_logit(Xs, ys, rep(1, length(ys))), pop, smp, ys, "A", NA))
  for (l in LAMBDAS) {
    sw <- wn * ifelse(ys == 1, l, 1)                  # survey weight x class weight
    out[[length(out) + 1]] <- evaluate(fit_logit(Xs, ys, sw), pop, smp, ys,
                                       ifelse(l == 1, "B", paste0("C", l)), l)
  }
  res <- do.call(rbind, out); res$rep <- rep_id
  diag <- data.frame(rep = rep_id, n_pos = sum(ys), weight_cv = sd(smp$w) / mean(smp$w),
                     kish_ess = sum(smp$w)^2 / sum(smp$w^2))
  list(res = res, diag = diag)
}

## ---- Demo: look inside ONE replicate -----------------------------------------
demo_one_replicate <- function(pop) {
  set.seed(1)
  smp <- draw_sample(pop); ys <- pop$y[smp$idx]
  cat("\nSample size n =", length(smp$idx), "| positives in sample =", sum(ys), "\n")
  cat("Weights: min =", round(min(smp$w), 1), " max =", round(max(smp$w), 1),
      " CV =", round(sd(smp$w) / mean(smp$w), 2), "\n")
  r <- one_replicate(pop, 1, smp)                     # fit all models on THIS sample
  tab <- r$res[, c("sens_t05", "fpr_t05", "calib_int", "calib_slope", "sens_tuned", "fpr_tuned")]
  rownames(tab) <- r$res$method
  cat("\nTrue performance of each model (population):\n"); print(round(tab, 3))
  invisible(r)
}

## Summary table: mean and Monte Carlo standard error of each metric, by method
summarise_results <- function(df) {
  cols <- c("sens_t05", "fpr_t05", "calib_int", "calib_slope", "sens_oracle", "sens_tuned", "fpr_tuned")
  do.call(rbind, lapply(split(df, factor(df$method, levels = unique(df$method))), function(d) {
    out <- data.frame(method = d$method[1])
    for (cc in cols) { out[[paste0(cc, "_mean")]] <- mean(d[[cc]]); out[[paste0(cc, "_se")]] <- sd(d[[cc]]) / sqrt(nrow(d)) }
    out
  }))
}

## ---- STEP 6. Run all replicates and summarise -------------------------------
run_all <- function(reps = 500, seed = SEED, outdir = "results_R") {
  dir.create(outdir, showWarnings = FALSE)
  pop <- make_population(seed); set.seed(seed + 1)
  all <- vector("list", reps); dg <- vector("list", reps); t0 <- Sys.time()
  for (r in seq_len(reps)) {
    o <- one_replicate(pop, r); all[[r]] <- o$res; dg[[r]] <- o$diag
    if (r %% 50 == 0) cat("replicate", r, "of", reps, "\n")
  }
  df <- do.call(rbind, all); dd <- do.call(rbind, dg)
  write.csv(df, file.path(outdir, "replicates.csv"), row.names = FALSE)
  summ <- summarise_results(df)
  write.csv(summ, file.path(outdir, "summary.csv"), row.names = FALSE)
  cat(sprintf("\nMean positives per sample: %.1f | mean weight CV: %.2f | mean Kish ESS: %.0f | non-converged fits: %d | %.1f min\n",
              mean(dd$n_pos), mean(dd$weight_cv), mean(dd$kish_ess), sum(df$nonconv),
              as.numeric(difftime(Sys.time(), t0, units = "mins"))))
  plot_results(df, outdir); invisible(list(df = df, summary = summ))
}

## ---- STEP 7. Figure ---------------------------------------------------------
plot_results <- function(df, outdir = "results_R") {
  R <- length(unique(df$rep)); ord <- c("B", "C2", "C5", "C10", "C20"); lam <- LAMBDAS
  ms <- function(col, m) { x <- df[[col]][df$method == m]; c(mean(x), 1.96 * sd(x) / sqrt(R)) }
  cur <- function(col) t(sapply(ord, function(m) ms(col, m)))   # rows = lambda; cols = mean, half-width
  panel <- function(col, main, ylab, ref = NULL, ylim = NULL) {
    v <- cur(col); a <- ms(col, "A")[1]
    if (is.null(ylim)) ylim <- range(c(v[, 1] - v[, 2], v[, 1] + v[, 2], a, ref))
    plot(lam, v[, 1], log = "x", type = "b", pch = 16, col = "#0072B2", xaxt = "n", ylim = ylim,
         xlab = "class weight lambda", ylab = ylab, main = main, cex.main = 0.9)
    axis(1, at = lam, labels = lam); ok <- v[, 2] > 1e-8; if (any(ok)) suppressWarnings(arrows(lam[ok], v[ok, 1] - v[ok, 2], lam[ok], v[ok, 1] + v[ok, 2], angle = 90, code = 3, length = 0.03, col = "#0072B2"))
    abline(h = a, lty = 2, col = "grey40"); if (!is.null(ref)) lines(lam, ref, lty = 3, col = "#D55E00", lwd = 2)
    grid()
  }
  png(file.path(outdir, "pilot_figure_R.png"), width = 2000, height = 1480, res = 200)
  par(mfrow = c(2, 2), mar = c(4, 4, 3, 1), oma = c(0, 0, 2, 0))
  panel("sens_t05", "(a) True sensitivity at threshold 0.5", "sensitivity")
  panel("fpr_t05", "(b) True FPR at threshold 0.5", "FPR")
  panel("calib_int", "(c) Calibration intercept (dotted: B intercept - log lambda)", "intercept",
        ref = ms("calib_int", "B")[1] - log(lam))
  panel("sens_tuned", "(d) True sensitivity at FPR ~ 5%, tuned threshold", "sensitivity", ylim = c(0, 0.3))
  mtext(sprintf("Pilot: prevalence 3%%, n = 2,000, %d replicates (mean +/- 1.96 MC SE; dashed = A)", R), outer = TRUE, line = 0.3, cex = 0.8)
  dev.off()
}

if (sys.nframe() == 0L) {                   # runs only under Rscript, not when sourced
  a <- commandArgs(trailingOnly = TRUE); run_all(reps = if (length(a)) as.integer(a[1]) else 500)
}
