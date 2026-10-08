## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width = 7,
  fig.height = 4.5,
  dpi = 150,
  out.width = "100%"
)

## ----library------------------------------------------------------------------
library(proxymix)

## ----engines------------------------------------------------------------------
has_ggplot2 <- requireNamespace("ggplot2", quietly = TRUE)

## ----stored-results, include = FALSE------------------------------------------
## The comparison table reads stored simulation results. They must come
## from the same major.minor version of proxymix as this build.
res <- readRDS("results/missing_data_mnar.rds")
major_minor <- function(v) paste(unlist(package_version(v))[1:2],
                                 collapse = ".")
if (major_minor(res$proxymix_version) !=
    major_minor(as.character(packageVersion("proxymix")))) {
  stop("results/missing_data_mnar.rds was built under proxymix ",
       res$proxymix_version, ", but this is proxymix ",
       packageVersion("proxymix"), ". Rerun the simulation and ",
       "data-raw/vignette_results/missing_data_mnar.R.", call. = FALSE)
}

## Small numbers are written as plain decimals rather than in the
## scientific notation that knitr's inline hook would otherwise use.
fixed <- function(v, digits) {
  format(round(v, digits), nsmall = digits, scientific = FALSE)
}

## ----mechanism-table, echo = FALSE--------------------------------------------
knitr::kable(
  data.frame(
    mechanism = c("missing at random", "censored", "missing not at random"),
    known = c("nothing beyond the rest of its row",
              "it lies beyond a known limit",
              "its chance of going missing depended on its size"),
    call = c("`mar()`", "`censored()`", "`mnar()`"),
    settled = c("yes, once the mixture is assumed",
                "yes, the limit is known",
                "only through the assumed shape of `y`"),
    stringsAsFactors = FALSE
  ),
  col.names = c("Mechanism", "What is known about a missing value",
                "Function", "Can the observed data fix the imputation model?"),
  caption = "The three mechanisms that `gmm_impute()` accepts."
)

## ----dgp----------------------------------------------------------------------
set.seed(20260622)
n <- 600L
comp <- sample(1:2, n, replace = TRUE)
mu <- rbind(c(0, 0), c(1.5, 0.5))
chol_r <- chol(matrix(c(1, 0.6, 0.6, 1), 2L))
z_full <- matrix(rnorm(2 * n), n, 2L) %*% chol_r + mu[comp, ]
colnames(z_full) <- c("x1", "y")
truth <- mean(z_full[, 2L])

beta_true <- 0.7
miss <- runif(n) < plogis(-0.5 + beta_true * z_full[, 2L])
dat <- z_full
dat[miss, "y"] <- NA

## ----recover------------------------------------------------------------------
m_draws <- 5L
mar_fit <- gmm_impute(dat, N = 2L, m = m_draws, mechanism = mar(),
                      seed = 1L, max_iter = 500L)
mnar_fit <- gmm_impute(dat, N = 2L, m = m_draws,
                       mechanism = mnar("y", beta = beta_true), seed = 1L,
                       max_iter = 500L)
mar_est <- proxy_pool(mar_fit, "y")$estimate
mnar_est <- proxy_pool(mnar_fit, "y", method = "rubin")$estimate

## ----recover-table, echo = FALSE----------------------------------------------
rec_tbl <- data.frame(
  data_used = c("complete data, before deletion", "rows not deleted",
                "imputed, missing at random",
                "imputed, missing not at random (slope 0.7)"),
  estimate = c(truth, mean(dat[!miss, "y"]), mar_est, mnar_est),
  stringsAsFactors = FALSE
)
rec_tbl$diff <- rec_tbl$estimate - truth
knitr::kable(
  rec_tbl, digits = 3L,
  col.names = c("Data used", "Mean of y", "Difference from complete data"),
  caption = paste0(
    "The mean of y from the complete data, from the rows not deleted, and ",
    "pooled over each imputation. ", sum(miss), " of ", n, " values of y ",
    "were deleted."
  )
)

## ----sweep--------------------------------------------------------------------
sweep <- proxy_mnar_sensitivity(dat, "y",
                                beta_grid = seq(0, 1.2, by = 0.3),
                                N = 2L, m = m_draws, seed = 1L)
covers <- sweep$conf.low <= truth & sweep$conf.high >= truth
beta_first <- sweep$beta[min(which(covers))]
beta_best <- sweep$beta[which.max(sweep$loglik)]
ll_gain <- max(sweep$loglik) - sweep$loglik[1L]

## ----sweep-table, echo = FALSE------------------------------------------------
sweep_tbl <- as.data.frame(sweep)[, c("beta", "estimate", "conf.low",
                                      "conf.high", "loglik", "converged")]
sweep_tbl$converged <- ifelse(sweep_tbl$converged, "yes", "no")
knitr::kable(
  sweep_tbl,
  digits = c(1L, 3L, 3L, 3L, 1L, 0L),
  col.names = c("Assumed slope", "Pooled mean of y", "CI lower",
                "CI upper", "Log-likelihood", "Converged"),
  caption = paste0(
    "The pooled mean of y at each assumed slope. The mean of y in the ",
    "complete data is ", round(truth, 3), "."
  )
)

## ----fig-sweep, eval = has_ggplot2, echo = has_ggplot2, fig.cap = "The pooled mean of y and its 95% confidence interval at each assumed slope. The dashed line is the mean of the complete data. The dotted line marks the slope that generated the data.", fig.alt = "Pooled mean of y against the assumed slope, rising from left to right, with a shaded confidence band, a point at each assumed slope, a dashed horizontal line at the complete-data mean, and a dotted vertical line at the generating slope of 0.7."----
sweep_df <- as.data.frame(sweep)
ggplot2::ggplot(sweep_df, ggplot2::aes(beta, estimate)) +
  ggplot2::geom_ribbon(
    ggplot2::aes(ymin = conf.low, ymax = conf.high),
    fill = "#56B4E9", alpha = 0.3
  ) +
  ggplot2::geom_line(colour = "#0072B2", linewidth = 0.9) +
  ggplot2::geom_point(colour = "#0072B2", size = 2) +
  ggplot2::geom_hline(yintercept = truth, linetype = "dashed",
                      colour = "#000000") +
  ggplot2::geom_vline(xintercept = beta_true, linetype = "dotted",
                      colour = "#D55E00") +
  ggplot2::annotate("text", x = min(sweep_df$beta), y = truth,
                    label = "complete-data mean", hjust = 0, vjust = -0.6,
                    size = 3.2) +
  ggplot2::annotate("text", x = beta_true, y = min(sweep_df$conf.low),
                    label = "generating slope", hjust = 1.05, vjust = 0,
                    size = 3.2, colour = "#D55E00") +
  ggplot2::labs(
    x = "assumed slope",
    y = "pooled mean of y",
    title = "The pooled mean of y under each assumed slope"
  ) +
  ggplot2::theme_minimal(base_size = 11)

## ----fig-sweep-skip, eval = !has_ggplot2, echo = FALSE, results = "asis"------
# cat("ggplot2 is not installed on this build, so the sensitivity figure",
#     "is skipped. The table above gives the same values.\n")

## ----censor-------------------------------------------------------------------
thr <- 0.3
cmiss <- z_full[, 2L] < thr
cdat <- z_full
cdat[cmiss, "y"] <- NA
cfit <- gmm_impute(cdat, N = 2L, m = m_draws,
                   mechanism = censored("y", upper = thr), seed = 1L)
cens_est <- proxy_pool(cfit, "y", method = "rubin")$estimate
half_est <- mean(ifelse(cmiss, thr / 2, z_full[, 2L]))
below_half <- mean(z_full[cmiss, 2L] < thr / 2)

## ----censor-table, echo = FALSE-----------------------------------------------
cens_tbl <- data.frame(
  data_used = c("complete data, before censoring", "rows not censored",
                "hidden values set to half the limit",
                "imputed, censored below the limit"),
  estimate = c(truth, mean(cdat[!cmiss, "y"]), half_est, cens_est),
  stringsAsFactors = FALSE
)
cens_tbl$diff <- cens_tbl$estimate - truth
knitr::kable(
  cens_tbl, digits = 3L,
  col.names = c("Data used", "Mean of y", "Difference from complete data"),
  caption = paste0(
    "The mean of y when values below ", thr, " are hidden. ", sum(cmiss),
    " of ", n, " values were hidden."
  )
)

## ----compare-facts, include = FALSE-------------------------------------------
sim_value <- function(design, estimand, method, what) {
  s1 <- res$sim_tab$design == design & res$sim_tab$estimand == estimand &
    res$sim_tab$method == method
  res$sim_tab[[what]][s1]
}
mv <- function(method, what) sim_value("mnar", "mean", method, what)
cv <- function(method, what) sim_value("censored", "slope", method, what)
mar_methods <- c("proxymix, slope 0", "mice, shift 0", "Amelia")
mar_bias <- vapply(mar_methods, mv, numeric(1L), what = "bias")
mar_cov <- vapply(mar_methods, mv, numeric(1L), what = "coverage")
and_list <- function(v) {
  paste(paste(v[-length(v)], collapse = ", "), "and", v[length(v)])
}

## ----compare-table, echo = FALSE----------------------------------------------
cmp_rows <- data.frame(
  design = c(rep("mnar", 6L), rep("censored", 4L)),
  method = c("complete data", "available cases", "Amelia",
             "proxymix, slope 0.7", "mice, shift 0.4",
             "sampleSelection heckit",
             "complete data", "zeros at face value", "AER tobit",
             "proxymix"),
  label = c("Mean of y: complete data, before deletion",
            "Mean of y: rows not deleted",
            "Mean of y: Amelia (missing at random)",
            "Mean of y: proxymix, slope 0.7",
            "Mean of y: mice, shift 0.4",
            "Mean of y: Heckman two-step, no exclusion restriction",
            "Slope: complete data, before censoring",
            "Slope: zeros taken at face value",
            "Slope: Tobit model (AER, survival)",
            "Slope: proxymix, censored imputation"),
  stringsAsFactors = FALSE
)
cmp_rows$estimand <- ifelse(cmp_rows$design == "mnar", "mean", "slope")
cmp_tbl <- data.frame(
  label = cmp_rows$label,
  bias = mapply(sim_value, cmp_rows$design, cmp_rows$estimand,
                cmp_rows$method, "bias"),
  rmse = mapply(sim_value, cmp_rows$design, cmp_rows$estimand,
                cmp_rows$method, "rmse"),
  coverage = mapply(sim_value, cmp_rows$design, cmp_rows$estimand,
                    cmp_rows$method, "coverage"),
  width = mapply(sim_value, cmp_rows$design, cmp_rows$estimand,
                 cmp_rows$method, "width"),
  stringsAsFactors = FALSE
)
knitr::kable(
  cmp_tbl, digits = 3L, row.names = FALSE,
  align = c("l", "r", "r", "r", "r"),
  col.names = c("Target and method", "Bias", "Error", "Coverage",
                "Interval width"),
  caption = paste0(
    "Results over ", res$n_rep, " simulated datasets per design. Bias is ",
    "the average difference from the true value, and error is the root ",
    "mean squared error. The mean of y is ", res$truth[["mnar_mean"]],
    " in the population, and the true slope is ",
    res$truth[["cens_slope"]], ". The proxymix and mice rows are the grid ",
    "values closest to the truth. With ", res$n_rep, " datasets, a ",
    "coverage near 0.95 has a simulation standard error of about ",
    fixed(sqrt(0.95 * 0.05 / res$n_rep), 3L), ". The Heckman coverage is ",
    "over the ", res$n_rep - res$n_undefined, " datasets in which its ",
    "estimated variance was positive."
  )
)

## ----compare-code, eval = FALSE-----------------------------------------------
# library(proxymix)
# library(mice)
# library(Amelia)
# library(AER)
# library(survival)
# library(sampleSelection)
# 
# # one dataset of 300 rows in which larger values of y are more often deleted
# set.seed(1L)
# n <- 300L
# comp <- sample(1:2, n, replace = TRUE)
# mu <- rbind(c(0, 0), c(1.5, 0.5))
# chol_r <- chol(matrix(c(1, 0.6, 0.6, 1), 2L))
# z <- matrix(rnorm(2 * n), n, 2L) %*% chol_r + mu[comp, ]
# full <- data.frame(x1 = z[, 1L], y = z[, 2L])
# obs <- full
# obs$y[runif(n) < plogis(-0.5 + 0.7 * full$y)] <- NA
# 
# # proxymix: the pooled mean of y at four assumed slopes
# proxy_mnar_sensitivity(obs, "y", beta_grid = c(0, 0.35, 0.7, 1.05),
#                        m = 10L, seed = 1L)
# 
# # mice: add a fixed shift to every imputed y, then pool the mean
# lapply(c(0, 0.2, 0.4, 0.6), function(delta) {
#   post <- make.post(obs)
#   post["y"] <- paste0("imp[[j]][, i] <- imp[[j]][, i] + ", delta)
#   imp <- mice(obs, m = 10L, method = "norm", post = post, seed = 1L,
#               printFlag = FALSE)
#   summary(pool(with(imp, lm(y ~ 1))), conf.int = TRUE)
# })
# 
# # Amelia: ten completed datasets, pooled by mice::pool()
# fits <- lapply(amelia(obs, m = 10L, p2s = 0L)$imputations,
#                function(d) lm(y ~ 1, data = d))
# summary(pool(fits), conf.int = TRUE)
# 
# # sampleSelection: Heckman two-step with x1 in both parts; the mean of y
# # is the fitted model for y at the mean of x1, with an approximate interval
# obs$seen <- !is.na(obs$y)
# fit <- heckit(seen ~ x1, y ~ x1, data = obs)
# x_bar <- c(1, mean(obs$x1))
# est <- sum(x_bar * coef(fit)[3:4])
# se <- sqrt(as.numeric(t(x_bar) %*% vcov(fit)[3:4, 3:4] %*% x_bar))
# c(estimate = est, conf.low = est - qnorm(0.975) * se,
#   conf.high = est + qnorm(0.975) * se)
# 
# # one dataset of 300 rows with y recorded as 0 whenever it falls below 0
# set.seed(1L)
# a_cens <- -sqrt(2) * qnorm(0.3)
# x <- rnorm(n)
# y_star <- a_cens + x + rnorm(n)
# cens <- data.frame(x = x, y = pmax(y_star, 0))
# 
# # proxymix: set the zeros to missing, draw them below 0, pool the regression
# holes <- cens
# holes$y[cens$y == 0] <- NA
# imp <- gmm_impute(holes, m = 10L, mechanism = censored("y", upper = 0),
#                   seed = 1L)
# summary(pool(lapply(complete(as_mids(imp), "all"),
#                     function(d) lm(y ~ x, data = d))), conf.int = TRUE)
# 
# # the Tobit model, fitted by AER and by survival
# coef(tobit(y ~ x, data = cens))
# coef(survreg(Surv(y, y > 0, type = "left") ~ x, data = cens,
#              dist = "gaussian"))

## ----session-info, collapse = FALSE, class.output = "session-info"------------
sessionInfo()

