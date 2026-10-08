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
has_mice <- requireNamespace("mice", quietly = TRUE)

## ----stored-results, include = FALSE------------------------------------------
## The comparison table reads stored simulation results. They must come
## from the same major.minor version of proxymix as this build.
res <- readRDS("results/missing_data.rds")
major_minor <- function(v) paste(unlist(package_version(v))[1:2],
                                 collapse = ".")
if (major_minor(res$proxymix_version) !=
    major_minor(as.character(packageVersion("proxymix")))) {
  stop("results/missing_data.rds was built under proxymix ",
       res$proxymix_version, ", but this is proxymix ",
       packageVersion("proxymix"), ". Rerun the simulation and ",
       "data-raw/vignette_results/missing_data.R.", call. = FALSE)
}

## Small numbers are written as plain decimals rather than in the
## scientific notation that knitr's inline hook would otherwise use.
fixed <- function(v, digits) {
  format(round(v, digits), nsmall = digits, scientific = FALSE)
}

## ----data---------------------------------------------------------------------
set.seed(20260620)
n <- 600L
lab <- sample(c(-1, 1), n, replace = TRUE)
x1 <- 2 * lab + rnorm(n, 0, 0.6)
x2 <- 2 * lab + 0.5 * (x1 - 2 * lab) + rnorm(n, 0, 0.6)
truth <- cbind(x1 = x1, x2 = x2)

x_holes <- truth
missing <- runif(n) < plogis(0.6 * x1)
x_holes[missing, "x2"] <- NA
frac_missing <- mean(missing)

## ----gap-truth----------------------------------------------------------------
in_gap <- function(v) mean(abs(v) < 1)
gap_truth <- in_gap(truth[missing, "x2"])

## ----impute-------------------------------------------------------------------
imp <- gmm_impute(x_holes, N = 2L, m = 20L, seed = 1L)
imp

## ----complete-----------------------------------------------------------------
done <- gmm_complete(imp, 1L)
anyNA(done)

## ----single-------------------------------------------------------------------
imp1 <- gmm_impute(x_holes, N = 1L, m = 20L, seed = 1L)
done1 <- gmm_complete(imp1, 1L)

## ----gap-table, echo = FALSE--------------------------------------------------
gap_tbl <- data.frame(
  source = c("deleted values (truth)", "mixture imputation (N = 2)",
             "single-normal imputation (N = 1)"),
  gap_share = c(gap_truth, in_gap(done[missing, "x2"]),
                in_gap(done1[missing, "x2"])),
  stringsAsFactors = FALSE
)
knitr::kable(
  gap_tbl, digits = 3L,
  col.names = c("Values", "Share with $\\lvert x_2 \\rvert < 1$"),
  caption = paste0(
    "Share of values in the empty middle, over the ", sum(missing),
    " deleted entries. The imputation rows use the first completed ",
    "dataset of each imputation."
  )
)

## ----fig-modes, eval = has_ggplot2, echo = has_ggplot2, fig.height = 3.8, fig.cap = "Density of the deleted values of $x_2$ and of the values each imputation put in their place, from the first completed dataset of each. The shaded band is the empty middle, $|x_2| < 1$.", fig.alt = "Three density curves over x2: the deleted values, the mixture imputation and the single-normal imputation. All three have two peaks. The single-normal curve has lower peaks and more density in the shaded band between the two groups."----
dens_df <- function(v, label) {
  d <- density(v)
  data.frame(x2 = d$x, density = d$y, source = label,
             stringsAsFactors = FALSE)
}
plot_df <- rbind(
  dens_df(truth[missing, "x2"], "deleted values (truth)"),
  dens_df(done[missing, "x2"], "mixture (N = 2)"),
  dens_df(done1[missing, "x2"], "single normal (N = 1)")
)
plot_df$source <- factor(
  plot_df$source,
  levels = c("deleted values (truth)", "mixture (N = 2)",
             "single normal (N = 1)")
)
ggplot2::ggplot(plot_df, ggplot2::aes(x2, density, colour = source)) +
  ggplot2::geom_line(linewidth = 0.9) +
  ggplot2::annotate("rect", xmin = -1, xmax = 1, ymin = -Inf, ymax = Inf,
                    fill = "grey60", alpha = 0.18) +
  ggplot2::scale_colour_manual(
    name = NULL,
    values = c("deleted values (truth)" = "#000000",
               "mixture (N = 2)" = "#009E73",
               "single normal (N = 1)" = "#D55E00")
  ) +
  ggplot2::labs(
    x = expression(x[2]), y = "density",
    title = "Deleted values and the values imputed in their place",
    subtitle = expression("Shaded band: the empty middle, " *
                            group("|", x[2], "|") < 1)
  ) +
  ggplot2::theme_minimal(base_size = 11) +
  ggplot2::theme(legend.position = "top")

## ----fig-modes-skip, eval = !has_ggplot2, echo = FALSE, results = "asis"------
# cat("ggplot2 is not installed on this build, so the imputation-density",
#     "figure is skipped.\n")

## ----pool-mean----------------------------------------------------------------
pooled <- proxy_pool(imp, "x2")
pooled1 <- proxy_pool(imp1, "x2")
fmi_mix <- proxy_fmi(imp, "x2")
mean_truth <- mean(truth[, "x2"])
err_mix <- abs(pooled$estimate - mean_truth)
err_one <- abs(pooled1$estimate - mean_truth)

## ----pool-table, echo = FALSE-------------------------------------------------
pool_tbl <- data.frame(
  model = c("complete data, before deletion", "rows not deleted",
            "mixture imputation (N = 2)",
            "single-normal imputation (N = 1)"),
  estimate = c(mean(truth[, "x2"]), mean(x_holes[!missing, "x2"]),
               pooled$estimate, pooled1$estimate),
  std_error = c(NA_real_, NA_real_, pooled$std.error, pooled1$std.error),
  conf_low = c(NA_real_, NA_real_, pooled$conf.low, pooled1$conf.low),
  conf_high = c(NA_real_, NA_real_, pooled$conf.high, pooled1$conf.high),
  stringsAsFactors = FALSE
)
pool_tbl$abs_error <- abs(pool_tbl$estimate - mean(truth[, "x2"]))
old_opt <- options(knitr.kable.NA = "")
pool_out <- knitr::kable(
  pool_tbl, digits = 4L,
  col.names = c("Data used", "Mean of $x_2$", "SE",
                "CI lower", "CI upper", "Distance from complete data"),
  caption = paste(
    "The mean of $x_2$ from the complete data, from the rows that were",
    "not deleted, and pooled over each imputation."
  )
)
options(old_opt)
pool_out

## ----pool-mice, eval = has_mice, echo = has_mice------------------------------
mice_fit <- mice::pool(with(as_mids(imp), lm(x2 ~ x1)))

## ----pool-mice-table, eval = has_mice, echo = FALSE---------------------------
mice_tbl <- summary(mice_fit)
mice_tbl$p.value <- ifelse(mice_tbl$p.value < 0.001, "< 0.001",
                           format(round(mice_tbl$p.value, 3L), nsmall = 3L))
knitr::kable(
  mice_tbl, digits = 3L, align = c("l", "r", "r", "r", "r", "r"),
  caption = paste0(
    "The regression of $x_2$ on $x_1$, fitted in each of the ", imp@m,
    " completed datasets of the mixture imputation and pooled by ",
    "`mice::pool()`."
  )
)

## ----pool-mice-skip, eval = !has_mice, echo = FALSE, results = "asis"---------
# cat("mice is not installed on this build, so the pooled regression is",
#     "skipped. The pooled mean above does not need mice.\n")

## ----compare-facts, include = FALSE-------------------------------------------
sim_value <- function(design, estimand, method, what) {
  s1 <- res$sim_tab$design == design & res$sim_tab$estimand == estimand &
    res$sim_tab$method == method
  res$sim_tab[[what]][s1]
}
cov_b <- function(method) fixed(sim_value("B", "slope", method, "coverage"), 3)
cov_a <- function(method) fixed(sim_value("A", "slope", method, "coverage"), 3)

## ----compare-table, echo = FALSE----------------------------------------------
methods <- c("complete data", "proxymix", "mice", "Amelia")
cmp_tbl <- data.frame(
  method = c("complete data, before deletion", "proxymix", "mice",
             "Amelia"),
  cov_a = vapply(methods, function(s1) {
    sim_value("A", "slope", s1, "coverage")
  }, numeric(1L)),
  cov_b = vapply(methods, function(s1) {
    sim_value("B", "slope", s1, "coverage")
  }, numeric(1L)),
  rmse_b = vapply(methods, function(s1) {
    sim_value("B", "slope", s1, "rmse")
  }, numeric(1L)),
  width_b = vapply(methods, function(s1) {
    sim_value("B", "slope", s1, "width")
  }, numeric(1L)),
  stringsAsFactors = FALSE
)
knitr::kable(
  cmp_tbl, digits = 3L, row.names = FALSE,
  align = c("l", "r", "r", "r", "r"),
  col.names = c("Data used", "Coverage, one cloud", "Coverage, two groups",
                "Error, two groups", "Interval width, two groups"),
  caption = paste0(
    "Slope of $x_2$ on $x_1$ over ", res$n_rep, " simulated datasets per ",
    "design. Coverage is the share of 95% intervals that contained the ",
    "true slope. Error is the root mean squared error of the estimate. ",
    "With ", res$n_rep, " datasets, a coverage near 0.95 has a simulation ",
    "standard error of about ", fixed(sqrt(0.95 * 0.05 / res$n_rep), 2L),
    "."
  )
)

## ----compare-code, eval = FALSE-----------------------------------------------
# library(proxymix)
# library(mice)
# library(Amelia)
# 
# # one dataset of 500 rows with two groups
# set.seed(1L)
# n <- 500L
# grp <- runif(n) < 0.5
# rho <- ifelse(grp, -0.3, 0.6)
# z1 <- rnorm(n)
# full <- data.frame(
#   x1 = ifelse(grp, 1.5, -1.5) + z1,
#   x2 = ifelse(grp, 2, -1.5) + rho * z1 + sqrt(1 - rho^2) * rnorm(n)
# )
# 
# # delete x2 with a probability that rises with x1
# obs <- full
# obs$x2[runif(n) < plogis(0.4 + 0.6 * full$x1)] <- NA
# 
# # 20 completed datasets from each package
# sets <- list(
#   proxymix = complete(as_mids(gmm_impute(obs, m = 20L, seed = 1L)), "all"),
#   mice = complete(mice(obs, m = 20L, seed = 1L, printFlag = FALSE), "all"),
#   Amelia = amelia(obs, m = 20L, p2s = 0L)$imputations
# )
# 
# # the same regression in every completed dataset, pooled by mice::pool()
# lapply(sets, function(s) {
#   fits <- lapply(s, function(d) lm(x2 ~ x1, data = d))
#   summary(pool(fits), conf.int = TRUE)
# })

## ----gap-all, include = FALSE-------------------------------------------------
## share in the empty middle, averaged over every completed dataset
gap_over <- function(im) {
  mean(vapply(seq_len(im@m), function(i1) {
    in_gap(gmm_complete(im, i1)[missing, "x2"])
  }, numeric(1L)))
}
gap_all <- c(mixture = gap_over(imp), single = gap_over(imp1))
## line and spread along which each imputation model draws x2 given x1,
## from the covariance matrix of each fitted component
model_line <- function(fit) {
  covs <- fit@covariances
  list(slope = vapply(covs, function(v) v[2L, 1L] / v[1L, 1L], numeric(1L)),
       sd = vapply(covs, function(v) sqrt(v[2L, 2L] - v[2L, 1L]^2 / v[1L, 1L]),
                   numeric(1L)))
}
model_mix <- model_line(imp@point_fit)
model_one <- model_line(imp1@point_fit)
fit_within <- lm(x2 ~ x1 + factor(lab), data = as.data.frame(truth))
slope_within <- coef(fit_within)[["x1"]]
sd_within <- sd(resid(fit_within))
## pooled estimate of the analyst's regression = mean over completions
pooled_slope <- function(im) {
  mean(vapply(seq_len(im@m), function(i1) {
    coef(lm(x2 ~ x1, data = as.data.frame(gmm_complete(im, i1))))[["x1"]]
  }, numeric(1L)))
}
fit_slope <- c(mixture = pooled_slope(imp), single = pooled_slope(imp1))
slope_complete <- coef(lm(x2 ~ x1, data = as.data.frame(truth)))[["x1"]]

## ----session-info, collapse = FALSE, class.output = "session-info"------------
sessionInfo()

