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
res <- readRDS("results/end_of_sample.rds")
major_minor <- function(v) paste(unlist(package_version(v))[1:2],
                                 collapse = ".")
if (major_minor(res$proxymix_version) !=
    major_minor(as.character(packageVersion("proxymix")))) {
  stop("results/end_of_sample.rds was built under proxymix ",
       res$proxymix_version, ", but this is proxymix ",
       packageVersion("proxymix"), ". Rerun the simulation and ",
       "data-raw/vignette_results/end_of_sample.R.", call. = FALSE)
}

## Small numbers are written as plain decimals rather than in the
## scientific notation that knitr's inline hook would otherwise use.
fixed <- function(v, digits) {
  format(round(v, digits), nsmall = digits, scientific = FALSE)
}

## ----seed---------------------------------------------------------------------
set.seed(20260621)

## ----model--------------------------------------------------------------------
prior <- gmm(weights = 1, means = list(0), covariances = list(matrix(10)))
dynamics <- list(A = matrix(1), Q = matrix(0.04))  # the level's random walk
measurement <- list(C = matrix(1), R = matrix(1))  # the observation noise

n <- 120L
level <- cumsum(c(0, rnorm(n - 1L, 0, sqrt(0.04))))
y_stable <- level + rnorm(n, 0, 1)

## ----break--------------------------------------------------------------------
y_break <- y_stable
y_break[n] <- y_break[n] + 5

test_stable <- gmm_eos_test(
  prior, dynamics, measurement, y_stable, m = 1L, method = "andrews"
)
test_break <- gmm_eos_test(
  prior, dynamics, measurement, y_break, m = 1L, method = "andrews"
)

## ----tests-kable, echo = FALSE------------------------------------------------
eos_row <- function(x) {
  c(
    round(x$statistic, 3L), round(x$p_value, 4L), x$reject,
    x$method, x$m
  )
}
knitr::kable(
  data.frame(
    field = c("statistic", "p-value", "reject at 0.05", "calibration",
              "window m"),
    stable = eos_row(test_stable),
    broken = eos_row(test_break)
  ),
  col.names = c("Result", "Stable series", "Broken final value"),
  caption = paste(
    "The end-of-sample test on the same series before and after its",
    "final value is moved up by five standard deviations of the",
    "observation noise."
  )
)

## ----filter-path--------------------------------------------------------------
filtered <- gmm_filter(
  prior, dynamics, measurement, y_stable, ridge_eps = 0
)
# ridge_eps adds a tiny amount to each covariance for numerical stability;
# 0 keeps the filter's recursion exact
path <- filtered$summary
# forecast of y_t and its standard deviation, from the filter at t - 1
pred_mean <- c(prior@means[[1L]], path$mean_1[-n])
pred_sd <- sqrt(c(prior@covariances[[1L]][1L, 1L], path$sd_1[-n]^2) +
                  dynamics$Q[1L, 1L] + measurement$R[1L, 1L])
outside <- sum(abs(y_stable - pred_mean)[-1L] > 2 * pred_sd[-1L])

## ----fig-series, eval = has_ggplot2, echo = has_ggplot2, fig.height = 4, fig.cap = "The simulated series, the unobserved level that generated it, the level estimated by the filter, and a band of two forecast standard deviations either side of each forecast. The band starts at the second time step. The broken final value lies far outside the band.", fig.alt = "A time series of noisy observations in grey with the true and estimated level overlaid, a shaded forecast band around them, and a single isolated point at the right-hand end lying well above the band."----
series_df <- data.frame(
  t = seq_len(n),
  observed = y_stable,
  truth = level,
  filtered = path$mean_1,
  lo = pred_mean - 2 * pred_sd,
  hi = pred_mean + 2 * pred_sd
)
break_df <- data.frame(t = n, y = y_break[n])
ggplot2::ggplot(series_df, ggplot2::aes(t)) +
  ggplot2::geom_ribbon(
    data = series_df[-1L, ], ggplot2::aes(ymin = lo, ymax = hi),
    fill = "#0072B2", alpha = 0.18
  ) +
  ggplot2::geom_point(
    ggplot2::aes(y = observed, colour = "observed (stable)"),
    size = 1.1, alpha = 0.7
  ) +
  ggplot2::geom_line(
    ggplot2::aes(y = truth, colour = "true level"), linewidth = 0.8
  ) +
  ggplot2::geom_line(
    ggplot2::aes(y = filtered, colour = "estimated level"), linewidth = 0.8
  ) +
  ggplot2::geom_point(
    data = break_df, ggplot2::aes(t, y, colour = "broken final value"),
    size = 2.6
  ) +
  ggplot2::scale_colour_manual(
    name = NULL,
    values = c(
      "observed (stable)" = "grey60",
      "true level" = "#009E73",
      "estimated level" = "#0072B2",
      "broken final value" = "#D55E00"
    )
  ) +
  ggplot2::labs(
    x = "time step", y = "observation",
    title = "A break in the last observation"
  ) +
  ggplot2::theme_minimal(base_size = 11) +
  ggplot2::theme(legend.position = "top")

## ----fig-series-skip, eval = !has_ggplot2, echo = FALSE, results = "asis"-----
# cat("ggplot2 is not installed on this build, so the series figure is",
#     "skipped.\n")

## ----fig-blocks, eval = has_ggplot2, echo = has_ggplot2, fig.height = 4, fig.cap = "The cumulative distribution of the earlier block statistics in each series, with the final block's statistic as a dashed line. The p-value is $(1 + k) / (n - 2m + 2)$, where $k$ is the number of earlier blocks at or beyond the dashed line. The stable series' final block lies inside the range of earlier blocks. The broken series' final block lies beyond every earlier block.", fig.alt = "Two panels, each showing a rising step function of the cumulative share of earlier block statistics against a logarithmic axis, with a dashed vertical line. In the left panel the line falls inside the curve. In the right panel it lies far beyond the curve's right-hand end."----
blocks_df <- rbind(
  data.frame(
    block = test_stable$in_sample_blocks, series = "stable series"
  ),
  data.frame(
    block = test_break$in_sample_blocks, series = "broken final value"
  )
)
blocks_df$series <- factor(
  blocks_df$series, levels = c("stable series", "broken final value")
)
rule_df <- data.frame(
  series = factor(
    c("stable series", "broken final value"),
    levels = c("stable series", "broken final value")
  ),
  statistic = c(test_stable$statistic, test_break$statistic)
)
ggplot2::ggplot(blocks_df, ggplot2::aes(block)) +
  ggplot2::stat_ecdf(geom = "step", linewidth = 0.8, colour = "#0072B2") +
  ggplot2::geom_vline(
    data = rule_df, ggplot2::aes(xintercept = statistic),
    colour = "#D55E00", linetype = "dashed", linewidth = 0.8
  ) +
  ggplot2::facet_wrap(~ series) +
  ggplot2::scale_x_log10(labels = function(v) {
    format(v, scientific = FALSE, drop0trailing = TRUE, trim = TRUE)
  }) +
  ggplot2::labs(
    x = "block statistic (log scale)",
    y = "cumulative share of earlier blocks",
    title = "The final block against the earlier blocks of the same series"
  ) +
  ggplot2::theme_minimal(base_size = 11)

## ----fig-blocks-skip, eval = !has_ggplot2, echo = FALSE, results = "asis"-----
# cat("ggplot2 is not installed on this build, so the block figure is",
#     "skipped.\n")

## ----size-study---------------------------------------------------------------
simulate_null <- function(n_obs) {
  lv <- cumsum(c(0, rnorm(n_obs - 1L, 0, sqrt(0.04))))
  lv + rnorm(n_obs, 0, 1)
}
empirical_size <- function(n_obs, method, n_rep) {
  p <- vapply(seq_len(n_rep), function(i1) {
    gmm_eos_test(
      prior, dynamics, measurement, simulate_null(n_obs),
      m = 1L, method = method
    )$p_value
  }, numeric(1L))
  size <- mean(p < 0.05)
  c(size = size, se = sqrt(size * (1 - size) / n_rep))
}

n_rep <- 250L
size_chisq_30 <- empirical_size(30L, "chisq", n_rep)
size_andrews_30 <- empirical_size(30L, "andrews", n_rep)

## ----grid-size----------------------------------------------------------------
grid_size <- function(n_obs, m = 1L, alpha = 0.05) {
  n_block <- n_obs - 2L * m + 1L
  p_grid <- (1L + seq.int(0L, n_block)) / (1L + n_block)
  mean(p_grid < alpha)
}
grid_30 <- grid_size(30L)
grid_120 <- grid_size(120L)

## ----size-checks, include = FALSE---------------------------------------------
## the prose below says the chi-square size does not differ detectably
## from 0.05 and the subsampling size agrees with its grid value
stopifnot(abs(size_chisq_30[["size"]] - 0.05) < 2 * size_chisq_30[["se"]],
          abs(size_andrews_30[["size"]] - grid_30) <
            2 * size_andrews_30[["se"]])

## ----size-kable, echo = FALSE-------------------------------------------------
old_opts <- options(knitr.kable.NA = "--")
knitr::kable(
  data.frame(
    calibration = c("chi-square", "subsampling", "subsampling",
                    "subsampling"),
    n_obs = c(30L, 30L, 30L, 120L),
    how = c("simulated", "simulated", "from the p-value grid",
            "from the p-value grid"),
    size = round(c(size_chisq_30[["size"]], size_andrews_30[["size"]],
                   grid_30, grid_120), 4L),
    se = round(c(size_chisq_30[["se"]], size_andrews_30[["se"]],
                 NA_real_, NA_real_), 4L)
  ),
  col.names = c("Calibration", "Series length", "Obtained", "Size",
                "Standard error"),
  caption = paste0(
    "Size at a nominal 0.05 with the model's parameters known, for m = 1. ",
    "Simulated sizes use ", n_rep, " stable series each. Sizes from the ",
    "p-value grid are exact when the final block's rank is equally likely ",
    "to be any rank, so they have no standard error."
  )
)
options(old_opts)

## ----compare-facts, include = FALSE-------------------------------------------
sim_value <- function(m, scenario, method) {
  s1 <- res$sim_tab$m == m & res$sim_tab$scenario == scenario &
    res$sim_tab$method == method
  res$sim_tab$rate[s1]
}
size_at <- function(method, m) fixed(sim_value(m, "no shift", method), 3)
power_at <- function(method, m) fixed(sim_value(m, "shift", method), 3)
mcse <- sqrt(0.05 * 0.95 / res$n_rep)
ms_of <- function(pkg) fixed(1000 * res$time_secs[[pkg]], 0)

## the prose below rests on these orderings in the stored results
size_1 <- vapply(unique(res$sim_tab$method), function(k) {
  sim_value(1L, "no shift", k)
}, numeric(1L))
stopifnot(identical(names(size_1)[!is.na(size_1) & abs(size_1 - 0.05) <= mcse],
                    "proxymix subsampling"))
stopifnot(abs(sim_value(5L, "no shift", "strucchange sup-F") - 0.05) <
            min(abs(sim_value(5L, "no shift", "proxymix chi-square") - 0.05),
                abs(sim_value(5L, "no shift", "proxymix subsampling") - 0.05)))
stopifnot(sim_value(1L, "shift", "proxymix subsampling") <
            min(sim_value(1L, "shift", "proxymix chi-square"),
                sim_value(1L, "shift", "KFAS forecast interval"),
                sim_value(1L, "shift", "dlm forecast interval")))
stopifnot(names(which.max(res$time_secs)) == "KFAS",
          res$time_secs[["strucchange"]] < 0.001,
          res$time_secs[["changepoint"]] < 0.001)

## ----compare-table, echo = FALSE----------------------------------------------
method_order <- c("proxymix chi-square", "proxymix subsampling",
                  "strucchange sup-F", "strucchange OLS-CUSUM",
                  "changepoint cpt.mean", "KFAS forecast interval",
                  "dlm forecast interval", "t-test on forecast residuals")
method_label <- c("proxymix, chi-square (default)", "proxymix, subsampling",
                  "strucchange, sup-F", "strucchange, OLS-CUSUM",
                  "changepoint, cpt.mean()", "KFAS, forecast interval",
                  "dlm, forecast interval", "t-test on forecast errors")
rate_col <- function(scenario, m) {
  v <- vapply(method_order, function(k) sim_value(m, scenario, k),
              numeric(1L))
  ifelse(is.na(v), "--", format(round(v, 3L), nsmall = 3L))
}
m_shown <- c(1L, 5L)
cmp_tbl <- data.frame(method = method_label, stringsAsFactors = FALSE)
for (m1 in m_shown) {
  cmp_tbl[[paste0("size_", m1)]] <- rate_col("no shift", m1)
  cmp_tbl[[paste0("power_", m1)]] <- rate_col("shift", m1)
}
knitr::kable(
  cmp_tbl, row.names = FALSE,
  align = c("l", rep("r", 2L * length(m_shown))),
  col.names = c("Check", as.vector(rbind(paste0("Size, m&nbsp;=&nbsp;", m_shown),
                                         paste0("Power, m&nbsp;=&nbsp;", m_shown)))),
  caption = paste0(
    "Rejection rates at a nominal 0.05 over ", res$n_rep, " simulated ",
    "series of ", res$n, " values, with parameters estimated from each ",
    "series. Size is the rate without a change and should be near 0.05. ",
    "Power is the rate when the last m values were raised by ", res$delta,
    " noise standard deviations. A size near 0.05 has a simulation ",
    "standard error of about ", fixed(mcse, 3), ". A dash marks a check ",
    "that cannot be computed at that m. The extended version of this ",
    "article also reports m = 3."
  )
)

## ----compare-code, eval = FALSE-----------------------------------------------
# library(proxymix)
# library(strucchange)
# library(changepoint)
# library(KFAS)
# library(dlm)
# 
# # one series of 100 values whose last m values are raised by 3
# set.seed(1L)
# n <- 100L
# m <- 3L
# y <- as.numeric(arima.sim(list(ar = 0.6), n = n))
# y[(n - m + 1L):n] <- y[(n - m + 1L):n] + 3
# y_fit <- y[seq_len(n - m)]
# last <- (n - m + 1L):n
# 
# # proxymix, with the autoregressive model fitted by arima()
# ar <- arima(y_fit, order = c(1L, 0L, 0L))
# a1 <- coef(ar)[["ar1"]]
# mu <- coef(ar)[["intercept"]]
# s2 <- ar$sigma2
# prior <- gmm(weights = 1, means = list(mu),
#              covariances = list(matrix(s2 / (1 - a1^2))))
# dynamics <- list(A = matrix(a1), b = mu * (1 - a1), Q = matrix(s2))
# measurement <- list(C = matrix(1), R = matrix(0))
# p_chisq <- gmm_eos_test(prior, dynamics, measurement, y, m = m,
#                         method = "chisq")$p_value
# p_andrews <- gmm_eos_test(prior, dynamics, measurement, y, m = m,
#                           method = "andrews")$p_value
# resid <- y[last] - predict(ar, n.ahead = m)$pred
# p_t <- if (m > 1L) t.test(resid)$p.value else NA_real_
# 
# # strucchange needs at least three observations after a break in an AR(1)
# reg <- data.frame(y = y[-1L], y_lag = y[-n])
# n_reg <- nrow(reg)
# p_supf <- if (m >= 3L) {
#   fs <- Fstats(y ~ y_lag, data = reg, from = n_reg - m, to = n_reg - 3L)
#   unname(sctest(fs)$p.value)
# } else NA_real_
# p_cusum <- unname(sctest(efp(y ~ y_lag, data = reg,
#                              type = "OLS-CUSUM"))$p.value)
# 
# cp <- cpts(cpt.mean(y / sd(y_fit)))
# p_cpt <- if (any(cp >= n - m)) 0 else 1
# 
# y_bar <- mean(y_fit)
# kfas_update <- function(pars, model) {
#   part <- SSMarima(ar = 0.999 * tanh(pars[1L]), Q = exp(pars[2L]))
#   model["T", "arima"] <- part$T
#   model["R", "arima"] <- part$R
#   model["Q", "arima"] <- part$Q
#   model["P1", "arima"] <- part$P1
#   model
# }
# kfas_fit <- fitSSM(
#   SSModel(I(y_fit - y_bar) ~ -1 + SSMarima(ar = 0.5, Q = 1), H = 0),
#   inits = c(0, 0), updatefn = kfas_update, method = "BFGS"
# )
# kfas_par <- kfas_fit$optim.out$par
# kfas_model <- SSModel(
#   I(y - y_bar) ~ -1 + SSMarima(ar = 0.999 * tanh(kfas_par[1L]),
#                                 Q = exp(kfas_par[2L])),
#   H = 0
# )
# kfas_pred <- predict(kfas_model, interval = "prediction", level = 0.95,
#                      filtered = TRUE)
# kfas_z <- (y[last] - y_bar - kfas_pred[last, "fit"]) /
#   ((kfas_pred[last, "upr"] - kfas_pred[last, "lwr"]) / (2 * qnorm(0.975)))
# p_kfas <- min(1, m * 2 * pnorm(-max(abs(kfas_z))))
# 
# dlm_build <- function(pars) {
#   dlmModARMA(ar = 0.999 * tanh(pars[1L]), sigma2 = exp(pars[2L]), dV = 0)
# }
# dlm_fit <- dlmMLE(y_fit - y_bar, parm = c(0, 0), build = dlm_build)
# dlm_filt <- dlmFilter(y - y_bar, dlm_build(dlm_fit$par))
# dlm_var <- unlist(dlmSvd2var(dlm_filt$U.R, dlm_filt$D.R))
# dlm_z <- (y[last] - y_bar - dlm_filt$f[last]) / sqrt(dlm_var[last])
# p_dlm <- min(1, m * 2 * pnorm(-max(abs(dlm_z))))
# 
# c("proxymix chi-square" = p_chisq, "proxymix subsampling" = p_andrews,
#   "strucchange sup-F" = p_supf, "strucchange OLS-CUSUM" = p_cusum,
#   "changepoint cpt.mean" = p_cpt, "KFAS forecast interval" = p_kfas,
#   "dlm forecast interval" = p_dlm, "t-test on forecast residuals" = p_t)

## ----session-info, collapse = FALSE, class.output = "session-info"------------
sessionInfo()

