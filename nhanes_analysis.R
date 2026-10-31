# ------------------------------------------------------------------
# NHANES application: survey-weighted vs class-weighted logistic models
# for undiagnosed diabetes (a rare outcome, about 2-3% of US adults).
# Base R + 'foreign' only.  Run:  source("nhanes_analysis.R")
# Set USE_SYNTHETIC <- TRUE to test the code without downloading data.
# ------------------------------------------------------------------
USE_SYNTHETIC <- FALSE
LAMBDAS <- c(1, 2, 5, 10, 20)   # lambda = 1 is Method B (survey-weighted)
ALPHA   <- 0.05                 # false-positive-rate limit for tuned threshold
K       <- 5                    # folds (whole PSUs are kept together)
REPEATS <- 10                   # repeated cross-validation
NBOOT   <- 500                  # PSU bootstrap replicates for uncertainty
set.seed(2027)
dir.create("results", showWarnings = FALSE)

# ---------------- 1. Data ------------------------------------------
get_xpt <- function(cycle_letter, year_start, file) {
  url <- sprintf("https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/%d/DataFiles/%s_%s.xpt",
                 year_start, file, cycle_letter)
  f <- file.path("data", sprintf("%s_%s.xpt", file, cycle_letter))
  dir.create("data", showWarnings = FALSE)
  if (!file.exists(f)) download.file(url, f, mode = "wb", quiet = TRUE)
  foreign::read.xport(f)
}
load_cycle <- function(letter, year_start) {
  d <- get_xpt(letter, year_start, "DEMO")[, c("SEQN","RIAGENDR","RIDAGEYR","RIDRETH3",
        "INDFMPIR","SDMVPSU","SDMVSTRA","WTMEC2YR")]
  for (nm in c("BMX","BPX","GHB","DIQ")) d <- merge(d, get_xpt(letter, year_start, nm), by = "SEQN", all.x = TRUE)
  d <- d[, c("SEQN","RIAGENDR","RIDAGEYR","RIDRETH3","INDFMPIR","SDMVPSU","SDMVSTRA",
             "WTMEC2YR","BMXBMI","BMXWAIST","BPXSY1","BPXDI1","LBXGH","DIQ010")]
  d$cycle <- letter
  d
}
make_synthetic <- function(n = 9000) {
  d <- data.frame(SEQN = 1:n, RIAGENDR = sample(1:2, n, TRUE), RIDAGEYR = sample(20:80, n, TRUE),
        RIDRETH3 = sample(c(1,2,3,4,6,7), n, TRUE), INDFMPIR = round(runif(n, 0, 5), 2),
        SDMVSTRA = sample(119:150, n, TRUE), SDMVPSU = sample(1:2, n, TRUE))
  d$WTMEC2YR <- exp(rnorm(n, 10, 0.7))
  d$BMXBMI <- rnorm(n, 22 + 0.06 * d$RIDAGEYR, 5); d$BMXWAIST <- 60 + 1.6 * d$BMXBMI + rnorm(n, 0, 6)
  d$BPXSY1 <- 100 + 0.4 * d$RIDAGEYR + rnorm(n, 0, 12); d$BPXDI1 <- 60 + rnorm(n, 0, 9)
  lp <- -9.2 + 0.05 * d$RIDAGEYR + 0.09 * d$BMXBMI + 0.01 * d$BPXSY1
  d$LBXGH <- ifelse(runif(n) < plogis(lp), 6.8, 5.4); d$DIQ010 <- sample(c(1,2,3), n, TRUE, c(.1,.85,.05))
  d$LBXGH[sample(n, 400)] <- NA; d$BPXSY1[sample(n, 700)] <- NA; d
}
raw <- if (USE_SYNTHETIC) make_synthetic() else rbind(load_cycle("I", 2015), load_cycle("J", 2017))

# Adults 20+, examined, with HbA1c, and no reported diagnosis of diabetes.
# Outcome: undiagnosed diabetes = HbA1c >= 6.5% with no diagnosis.
d <- raw[raw$RIDAGEYR >= 20 & !is.na(raw$LBXGH) & raw$DIQ010 %in% c(2, 3) & raw$WTMEC2YR > 0, ]
d$y <- as.integer(d$LBXGH >= 6.5)
d$w <- d$WTMEC2YR / (if (USE_SYNTHETIC) 1 else 2)          # two cycles combined
d$cluster <- paste(d$SDMVSTRA, d$SDMVPSU)
# One documented missing-data rule: median imputation for predictors (sensitivity check later)
xs <- c("RIDAGEYR","BMXBMI","BMXWAIST","BPXSY1","BPXDI1","INDFMPIR")
for (v in xs) d[[v]][is.na(d[[v]])] <- median(d[[v]], na.rm = TRUE)
X <- cbind(scale(d[, xs]), female = as.numeric(d$RIAGENDR == 2),
           model.matrix(~ factor(RIDRETH3), d)[, -1])
X <- cbind(1, X); colnames(X)[1] <- "intercept"

cat("Adults analysed:", nrow(d), "| undiagnosed cases:", sum(d$y),
    "| weighted prevalence:", round(sum(d$w * d$y) / sum(d$w), 4),
    "| Kish effective n:", round(sum(d$w)^2 / sum(d$w^2)), "\n")

# ---------------- 2. Functions -------------------------------------
fit <- function(X, y, w, lambda) {            # weight = survey weight * lambda^y
  ww <- w * ifelse(y == 1, lambda, 1); ww <- ww / mean(ww)
  suppressWarnings(glm.fit(X, y, weights = ww, family = quasibinomial())$coefficients)
}
prob <- function(b, X) plogis(pmin(pmax(as.vector(X %*% b), -30), 30))
wquant <- function(x, w, q) { o <- order(x); cw <- cumsum(w[o]) / sum(w); x[o][which(cw >= q)[1]] }
tuned_t <- function(p, y, w, alpha) wquant(p[y == 0], w[y == 0], 1 - alpha)
wauc <- function(p, y, w) {
  o <- order(p); p <- p[o]; y <- y[o]; w <- w[o]
  negcum <- cumsum(w * (y == 0)); neg_before <- negcum - w * (y == 0)
  # ties: treat as half credit by averaging within tied groups
  grp <- match(p, unique(p)); gneg <- tapply(w * (y == 0), grp, sum)
  gstart <- tapply(neg_before, grp, min)
  credit <- gstart[grp] + 0.5 * gneg[grp]
  sum(w * (y == 1) * credit) / (sum(w * (y == 1)) * sum(w * (y == 0)))
}
metrics <- function(p, y, w, t) {
  pr <- as.numeric(p >= t)
  c(sens = sum(w * y * pr) / sum(w * y), fpr = sum(w * (1 - y) * pr) / sum(w * (1 - y)))
}
calib <- function(p, y, w) {
  z <- qlogis(pmin(pmax(p, 1e-6), 1 - 1e-6)); ww <- w / mean(w)
  suppressWarnings(coef(glm.fit(cbind(1, z), y, weights = ww, family = quasibinomial())))
}

# ---------------- 3. Cross-validation by whole PSUs ----------------
one_repeat <- function(seed) {
  set.seed(seed)
  cl <- unique(d$cluster); fold_of <- setNames(sample(rep_len(1:K, length(cl))), cl)
  fold <- fold_of[d$cluster]
  P <- matrix(NA_real_, nrow(d), length(LAMBDAS)); Tt <- P    # predictions, tuned thresholds
  for (k in 1:K) {
    tr <- fold != k; te <- fold == k
    for (j in seq_along(LAMBDAS)) {
      b <- fit(X[tr, ], d$y[tr], d$w[tr], LAMBDAS[j])
      P[te, j] <- prob(b, X[te, ])
      Tt[te, j] <- tuned_t(prob(b, X[tr, ]), d$y[tr], d$w[tr], ALPHA)   # t chosen on training data only
    }
  }
  list(P = P, Tt = Tt)
}
summ <- function(P, Tt, idx = seq_len(nrow(d))) {
  y <- d$y[idx]; w <- d$w[idx]
  t(sapply(seq_along(LAMBDAS), function(j) {
    m05 <- metrics(P[idx, j], y, w, 0.5)
    # threshold varies by fold; use pooled-fold mean for a single tuned value
    mt <- metrics(P[idx, j], y, w, mean(Tt[idx, j]))
    c(lambda = LAMBDAS[j], sens05 = m05[["sens"]], fpr05 = m05[["fpr"]],
      sens_tuned = mt[["sens"]], fpr_tuned = mt[["fpr"]],
      auroc = wauc(P[idx, j], y, w), calib_int = calib(P[idx, j], y, w)[[1]],
      calib_slope = calib(P[idx, j], y, w)[[2]])
  }))
}
cat("Running", REPEATS, "x", K, "fold cross-validation...\n")
reps <- lapply(1:REPEATS, function(r) one_repeat(1000 + r))
tabs <- lapply(reps, function(r) summ(r$P, r$Tt))
res <- Reduce(`+`, tabs) / REPEATS
res <- as.data.frame(res); res$method <- ifelse(res$lambda == 1, "B (survey-weighted)", paste0("C, lambda=", res$lambda))
write.csv(res, "results/nhanes_main_table.csv", row.names = FALSE)
# unweighted method A for reference
fitA <- function(tr) fit(X[tr, ], d$y[tr], rep(1, sum(tr)), 1)

# ---------------- 4. PSU bootstrap (paired difference C20 - B) -----
# Uses all repeats: in each bootstrap draw, resample PSUs within strata, compute the
# difference for every repeat, and average over repeats.
cat("PSU bootstrap for the paired difference at the tuned threshold...\n")
idx_by_psu <- split(seq_len(nrow(d)), d$cluster)
psu_by_strat <- split(unique(d$cluster), sapply(strsplit(unique(d$cluster), " "), `[`, 1))
jB <- which(LAMBDAS == 1); jC <- which(LAMBDAS == 20)
pair_diff <- function(idx) {
  mean(sapply(reps, function(r) {
    y <- d$y[idx]; w <- d$w[idx]
    metrics(r$P[idx, jC], y, w, mean(r$Tt[, jC]))[["sens"]] -
      metrics(r$P[idx, jB], y, w, mean(r$Tt[, jB]))[["sens"]]
  }))
}
diff_est <- pair_diff(seq_len(nrow(d)))
boot <- replicate(NBOOT, {
  pick <- unlist(lapply(psu_by_strat, function(u) sample(u, length(u), replace = TRUE)))
  pair_diff(unlist(idx_by_psu[pick]))
})
ci <- quantile(boot, c(0.025, 0.975), na.rm = TRUE)
out <- c(paired_diff_sens_tuned = diff_est, ci_low = ci[[1]], ci_high = ci[[2]])
write.csv(data.frame(t(out)), "results/nhanes_paired_difference.csv", row.names = FALSE)

# ---------------- 5. Figure and printed summary --------------------
png("results/nhanes_figure.png", 1400, 450, res = 130); par(mfrow = c(1, 3), mar = c(4.2, 4.2, 2.5, 1))
plot(res$lambda, res$sens05, type = "b", log = "x", ylim = c(0, 1), xlab = "class weight lambda (log scale)",
     ylab = "sensitivity / FPR", main = "(a) at threshold 0.5", pch = 16)
lines(res$lambda, res$fpr05, type = "b", lty = 2, pch = 1); legend("topleft", c("sensitivity", "FPR"), lty = 1:2, pch = c(16, 1), bty = "n")
plot(res$lambda, res$calib_int, type = "b", log = "x", xlab = "class weight lambda (log scale)",
     ylab = "calibration intercept", main = "(b) calibration", pch = 16)
lines(res$lambda, res$calib_int[1] - log(res$lambda), lty = 3)
plot(res$lambda, res$sens_tuned, type = "b", log = "x", ylim = c(0, max(res$sens_tuned) * 1.4),
     xlab = "class weight lambda (log scale)", ylab = "sensitivity at FPR ~5%", main = "(c) threshold tuned", pch = 16)
dev.off()
options(width = 160); cat("\n"); print(cbind(method = res$method, round(res[, c("sens05","fpr05","sens_tuned","fpr_tuned","auroc","calib_int","calib_slope")], 3)), row.names = FALSE)
cat("\nPaired difference in sensitivity (C lambda=20 minus B, tuned threshold): ",
    round(diff_est, 4), " 95% PSU-bootstrap CI [", round(ci[[1]], 4), ",", round(ci[[2]], 4), "]\n")
writeLines(capture.output(sessionInfo()), "results/sessionInfo.txt")
