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
res <- readRDS("results/quickstart.rds")
major_minor <- function(v) paste(unlist(package_version(v))[1:2],
                                 collapse = ".")
if (major_minor(res$proxymix_version) !=
    major_minor(as.character(packageVersion("proxymix")))) {
  stop("results/quickstart.rds was built under proxymix ",
       res$proxymix_version, ", but this is proxymix ",
       packageVersion("proxymix"), ". Rerun the simulation and ",
       "data-raw/vignette_results/quickstart.R.", call. = FALSE)
}

## Small numbers are written as plain decimals rather than in the
## scientific notation that knitr's inline hook would otherwise use.
fixed <- function(v, digits) {
  format(round(v, digits), nsmall = digits, scientific = FALSE)
}

## ----regime-table, echo = FALSE-----------------------------------------------
regime_tbl <- data.frame(
  have = c(
    "Data, and you want one bell curve",
    "Data, and you want several bell curves",
    "Only the formula, no data"
  ),
  does = c(
    "Fits one bell curve with the same centre and spread as the data",
    paste("Moves and reshapes the bell curves, step by step, until together",
          "they match the data"),
    paste("Draws trial points, weights each one by the formula, then fits",
          "the bell curves to the weighted points")
  ),
  setting = c("`\"moment\"`", "`\"sample\"`", "`\"kld\"`"),
  stringsAsFactors = FALSE
)
knitr::kable(
  regime_tbl,
  col.names = c("What you have", "What the package does", "`regime` setting"),
  caption = paste(
    "How `fit_proxymix()` chooses a fitting method. The `regime` argument",
    "can also name one directly."
  )
)

## ----banana-------------------------------------------------------------------
tgt <- banana_target()
tgt

## ----fit----------------------------------------------------------------------
proposal <- proposal_mvt(n_dim = 2L, mean = c(0, 0),
                         sigma = 4 * diag(2), df = 5)
fit <- fit_proxymix(tgt, N = 3L, regime = "kld",
                    proposal = proposal,
                    is_size = 2000L,
                    max_iter = 60L,
                    seed = 1L)
fit

## ----certificate--------------------------------------------------------------
quality <- gmm_fit_quality(fit)

## ----certificate-table, echo = FALSE------------------------------------------
cert_tbl <- data.frame(
  Check = c("fitting method", "rounds settled before the limit",
            "weights collapsed onto a few draws",
            "effective sample size",
            "effective sample size as a share of all draws",
            "smallest effective sample size of any component",
            "largest share of the weight held by one draw",
            "KL divergence on the fitting draws",
            "half the variance of the log density ratio (a local KL approximation)",
            "KL divergence on fresh draws",
            "fresh draws minus fitting draws"),
  Value = c(
    quality$regime,
    as.character(quality$converged),
    as.character(quality$degenerate),
    format(round(quality$ess, 1L), nsmall = 1L),
    format(round(quality$ess_relative, 3L), nsmall = 3L),
    format(round(quality$min_component_ess, 1L), nsmall = 1L),
    format(signif(quality$max_weight, 3L)),
    format(signif(quality$kld_final, 3L)),
    format(signif(quality$kld_approx, 3L)),
    format(signif(quality$heldout_kld, 3L)),
    format(signif(quality$validation_gap, 3L))
  ),
  stringsAsFactors = FALSE
)
knitr::kable(
  cert_tbl,
  caption = "The fit certificate returned by `gmm_fit_quality()`."
)

## ----overlay-grid-------------------------------------------------------------
grid_x <- seq(-3, 3, length.out = 120L)
grid_g <- expand.grid(x1 = grid_x, x2 = grid_x)
grid_mat <- as.matrix(grid_g)
grid_g$target <- exp(tgt@log_density(grid_mat))
grid_g$proxy <- dgmm(grid_mat, fit)

## ----overlay, eval = has_ggplot2, echo = has_ggplot2, fig.cap = sprintf("The banana target (filled contours) with the three-component proxy overlaid as dashed contours. The dashed contours follow the curve of the banana, which a single bell curve could not do. The KL divergence on fresh draws is %s.", format(signif(fit@diagnostics$validation_kld, 2L))), fig.alt = "Filled contour map of the curved banana density with dashed contours of the three-component Gaussian-mixture proxy following the same curve."----
ggplot2::ggplot(grid_g, ggplot2::aes(x1, x2)) +
  ggplot2::geom_contour_filled(ggplot2::aes(z = target), bins = 10L,
                               alpha = 0.85) +
  ggplot2::geom_contour(ggplot2::aes(z = proxy), colour = "white",
                        linetype = "dashed", linewidth = 0.45, bins = 8L) +
  ggplot2::scale_fill_viridis_d(option = "mako", guide = "none") +
  ggplot2::coord_equal() +
  ggplot2::labs(
    title = "Target (filled) and fitted proxy (dashed)",
    x = expression(x[1]), y = expression(x[2])
  ) +
  ggplot2::theme_minimal(base_size = 11)

## ----overlay-skip, eval = !has_ggplot2, echo = FALSE, results = "asis"--------
# cat("ggplot2 is not installed on this build, so the target-and-proxy",
#     "overlay figure is skipped.\n")

## ----operations---------------------------------------------------------------
gmm_marginalise(fit, keep = 1L)
gmm_conditionalise(fit, given = c(NA, 0.5))

## ----sample-------------------------------------------------------------------
draws <- rgmm(500L, fit)
dim(draws)

## ----compare-facts, include = FALSE-------------------------------------------
sim_value <- function(method, what) {
  res$sim_tab[[what]][res$sim_tab$method == method]
}
kl_value <- function(method) fixed(sim_value(method, "kl_mean"), 4)
tail_error <- function(method) fixed(sim_value(method, "tail_rmse"), 4)
secs <- function(method) fixed(sim_value(method, "secs"), 2)
tail_others <- vapply(c("NUTS", "NUTS + mclust", "DEzs"), sim_value,
                      numeric(1L), what = "tail_rmse")
tail_range <- paste(fixed(min(tail_others), 4), "to",
                    fixed(max(tail_others), 4))
# the three-component fit from above, scored on the simulation's grid
quad_grid <- as.matrix(expand.grid(x1 = seq(-5 + res$h / 2, 5, by = res$h),
                                   x2 = seq(-5 + res$h / 2, 12, by = res$h)))
quad_log_f <- tgt@log_density(quad_grid)
kl_fit_quad <- sum(exp(quad_log_f) *
                     (quad_log_f - dgmm(quad_grid, fit, log = TRUE))) * res$h^2

## ----compare-table, echo = FALSE----------------------------------------------
methods <- c("proxymix", "NUTS + mclust", "Laplace", "NUTS", "DEzs")
cmp_tbl <- data.frame(
  method = c("proxymix", "NUTS draws, then mclust", "Laplace approximation",
             "NUTS draws", "DEzs draws"),
  kl = vapply(methods, function(s1) {
    if (is.na(sim_value(s1, "kl_mean"))) return("--")
    paste0(kl_value(s1), " (", fixed(sim_value(s1, "kl_sd"), 4), ")")
  }, character(1L)),
  tail = vapply(methods, function(s1) {
    fixed(sim_value(s1, "tail_rmse"), 5)
  }, character(1L)),
  secs = vapply(methods, function(s1) {
    if (sim_value(s1, "secs") < 0.01) return("< 0.01")
    secs(s1)
  }, character(1L)),
  stringsAsFactors = FALSE
)
knitr::kable(
  cmp_tbl, row.names = FALSE,
  align = c("l", "r", "r", "r"),
  col.names = c("Method", "KL divergence, mean (sd)",
                "Error in $P(x_1 > 2)$", "Seconds per run"),
  caption = paste0(
    "Results over ", res$n_rep, " runs of each method on the banana ",
    "target. The KL divergence is averaged over the runs, with its ",
    "standard deviation in brackets. The error in $P(x_1 > 2)$ is the root ",
    "mean squared error over the runs. Seconds per run is the median; the ",
    "proxymix time includes choosing the number of components, and the ",
    "NUTS times leave out the one-off compilation of the Stan program."
  )
)

## ----compare-code, eval = FALSE-----------------------------------------------
# library(proxymix)
# library(cmdstanr)
# library(mclust)
# library(BayesianTools)
# 
# tgt <- banana_target()
# 
# # midpoint-rule quadrature; the target's mass outside the box is below 1e-5
# h <- 0.05
# grid <- as.matrix(expand.grid(x1 = seq(-5 + h / 2, 5, by = h),
#                               x2 = seq(-5 + h / 2, 12, by = h)))
# log_f <- tgt@log_density(grid)
# f_grid <- exp(log_f)
# tail_ref <- sum(f_grid[grid[, 1L] > 2]) * h^2
# kl_of <- function(log_g) sum(f_grid * (log_f - log_g)) * h^2
# 
# # mass above x1 = 2 under a Gaussian mixture, from its x1 marginal
# tail_of <- function(w, mean1, sd1) {
#   sum(w * pnorm(2, mean1, sd1, lower.tail = FALSE))
# }
# 
# # the same target as a Stan program
# stan_file <- file.path(tempdir(), "banana.stan")
# writeLines(c(
#   "parameters {",
#   "  vector[2] x;",
#   "}",
#   "model {",
#   "  target += -0.5 * (square(x[1])",
#   "                    + square(x[2] - 0.5 * (square(x[1]) - 1)));",
#   "}"), stan_file)
# banana_stan <- cmdstan_model(stan_file)
# 
# box <- createBayesianSetup(likelihood = function(x) tgt@log_density(x),
#                            lower = c(-6, -6), upper = c(6, 14))
# 
# set.seed(1L)
# fit <- select_N(tgt, seed = 1L)$best_fit
# fit_1 <- gmm_marginalise(fit, keep = 1L)
# 
# la <- optim(c(0, 0), function(x) -tgt@log_density(x),
#             method = "BFGS", hessian = TRUE)
# la_mean <- la$par
# la_cov <- solve(la$hessian)
# 
# nuts <- banana_stan$sample(seed = 1L, refresh = 0L, show_messages = FALSE)
# draws <- nuts$draws("x", format = "matrix")
# 
# mc <- Mclust(draws, G = 1:6, modelNames = "VVV", verbose = FALSE)
# mc_par <- mc$parameters
# 
# de <- runMCMC(box, sampler = "DEzs", settings = list(message = FALSE))
# de_draws <- getSample(de, start = 1000L)
# 
# # KL divergence from the target, for the methods that return a density
# c(proxymix = kl_of(dgmm(grid, fit, log = TRUE)),
#   Laplace = kl_of(dmvnorm(grid, la_mean, la_cov, log = TRUE)),
#   "NUTS + mclust" = kl_of(dens(grid, mc$modelName, parameters = mc_par,
#                                logarithm = TRUE)))
# 
# # error in the estimated probability that x1 > 2
# c(proxymix = tail_of(gmm_weights(fit_1),
#                      vapply(gmm_means(fit_1), `[[`, numeric(1L), 1L),
#                      sqrt(vapply(gmm_covariances(fit_1), `[[`,
#                                  numeric(1L), 1L))),
#   Laplace = pnorm(2, la_mean[1L], sqrt(la_cov[1L, 1L]), lower.tail = FALSE),
#   NUTS = mean(draws[, 1L] > 2),
#   "NUTS + mclust" = tail_of(mc_par$pro, mc_par$mean[1L, ],
#                             sqrt(mc_par$variance$sigma[1L, 1L, ])),
#   DEzs = mean(de_draws[, 1L] > 2)) - tail_ref

## ----session-info, collapse = FALSE, class.output = "session-info"------------
sessionInfo()

