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
res <- readRDS("results/posterior_proxy.rds")
major_minor <- function(v) paste(unlist(package_version(v))[1:2],
                                 collapse = ".")
if (major_minor(res$proxymix_version) !=
    major_minor(as.character(packageVersion("proxymix")))) {
  stop("results/posterior_proxy.rds was built under proxymix ",
       res$proxymix_version, ", but this is proxymix ",
       packageVersion("proxymix"), ". Rerun the simulation and ",
       "data-raw/vignette_results/posterior_proxy.R.", call. = FALSE)
}

## Small numbers are written as plain decimals rather than in the
## scientific notation that knitr's inline hook would otherwise use.
fixed <- function(v, digits) {
  format(round(v, digits), nsmall = digits, scientific = FALSE)
}
## A small probability as "1 in N", with N rounded to two figures.
one_in <- function(prob) {
  format(signif(1 / prob, 2L), big.mark = ",", scientific = FALSE)
}

## ----posterior----------------------------------------------------------------
set.seed(20260705)

y <- mtcars$am
w <- mtcars$wt

log_post <- function(theta) {
  if (is.null(dim(theta))) theta <- matrix(theta, ncol = 2L)
  eta <- outer(rep(1, length(y)), theta[, 1L]) + outer(w, theta[, 2L])
  colSums(y * eta - log1p(exp(eta)))
}

tgt <- gmm_target(
  n_dim = 2L,
  log_density = log_post,
  normalised = FALSE,
  name = "logistic(am ~ wt)"
)

## ----fit----------------------------------------------------------------------
mle <- stats::glm(am ~ wt, data = mtcars, family = stats::binomial())
q0 <- proposal_mvt(2L, mean = stats::coef(mle),
                   sigma = 9 * stats::vcov(mle), df = 5)

fit_notes <- character(0)
fit <- withCallingHandlers(
  fit_proxymix(tgt, N = 2L, regime = "kld", proposal = q0,
               is_size = 3000L, max_iter = 60L, seed = 1L,
               adapt = "pmc"),
  proxymix_nonmonotone = function(cond) {
    fit_notes <<- c(fit_notes, conditionMessage(cond))
    invokeRestart("muffleWarning")
  }
)
cat(sub("\n.*", "", fit_notes), sep = "\n")

## ----certificate--------------------------------------------------------------
cert <- gmm_fit_quality(fit)

## ----certificate-table, echo = FALSE------------------------------------------
cert_tbl <- data.frame(
  Check = c("fitting method", "rounds settled before the limit",
            "weights collapsed onto a few draws",
            "effective sample size",
            "effective sample size as a share of all draws",
            "smallest effective sample size of any component",
            "largest share of the weight held by one draw",
            "share of draws where the posterior could be evaluated"),
  Value = c(
    cert$regime,
    as.character(cert$converged),
    as.character(cert$degenerate),
    format(round(cert$ess, 1L), nsmall = 1L),
    format(round(cert$ess_relative, 3L), nsmall = 3L),
    format(round(cert$min_component_ess, 1L), nsmall = 1L),
    format(signif(cert$max_weight, 3L), scientific = FALSE),
    format(cert$support_fraction)
  ),
  stringsAsFactors = FALSE
)
knitr::kable(
  cert_tbl,
  caption = "The fit certificate returned by `gmm_fit_quality()`."
)

## ----evidence-----------------------------------------------------------------
ev <- gmm_evidence(fit, n = 4000L, seed = 2L)

## Laplace approximation: log f(theta_hat) + (d/2) log(2 pi)
## - (1/2) log det(-Hessian).
H <- -solve(stats::vcov(mle))
log_z_laplace <- log_post(matrix(stats::coef(mle), nrow = 1L)) +
  log(2 * pi) -
  0.5 * as.numeric(determinant(-H, logarithm = TRUE)$modulus)

c(proxymix = round(ev$log_z, 3), se = signif(ev$se_log_z, 2),
  laplace = round(log_z_laplace, 3))

## ----kl-fresh-----------------------------------------------------------------
kl_fresh <- fit@diagnostics$validation_kld - ev$log_z
# the held-out estimate and log Z come from separate sets of draws
kl_fresh_se <- sqrt(fit@diagnostics$validation_mc_se^2 + ev$se_log_z^2)
c(kl = signif(kl_fresh, 2), se = signif(kl_fresh_se, 1))

## ----reads--------------------------------------------------------------------
slope <- gmm_marginalise(fit, keep = 2L)
p_negative <- pgmm(0, slope)
interval <- qgmm(c(0.05, 0.95), slope)
c(p_slope_negative = round(p_negative, 5),
  lower = round(interval[1L], 2), upper = round(interval[2L], 2))

## ----quadrature---------------------------------------------------------------
a_grid <- seq(-2, 70, length.out = 600L)     # intercept
b_grid <- seq(-21, 0, length.out = 400L)     # slope
quad <- as.matrix(expand.grid(a = a_grid, b = b_grid))

log_dens <- log_post(quad)
dens <- matrix(exp(log_dens - max(log_dens)), nrow = length(a_grid))
da <- a_grid[2L] - a_grid[1L]
db <- b_grid[2L] - b_grid[1L]
log_z_grid <- log(sum(dens) * da * db) + max(log_dens)

marg_quad <- colSums(dens) * da
marg_quad <- marg_quad / (sum(marg_quad) * db)
marg_proxy <- dgmm(matrix(b_grid, ncol = 1L), slope)

cdf_quad <- cumsum(marg_quad) * db
interval_quad <- stats::approx(cdf_quad, b_grid + db / 2,
                               xout = c(0.05, 0.95), ties = mean)$y
curve_gap <- max(abs(marg_quad - marg_proxy))
c(log_z = round(log_z_grid, 3), lower = round(interval_quad[1L], 2),
  upper = round(interval_quad[2L], 2),
  gap_pct_of_peak = round(100 * curve_gap / max(marg_quad), 1))

## ----grid-edge, include = FALSE-----------------------------------------------
## The log-density on the edges of the grid, relative to its peak.
edge_vec <- c(log_dens[quad[, "a"] %in% range(a_grid)],
              log_dens[quad[, "b"] %in% range(b_grid)])
edge_drop <- max(edge_vec) - max(log_dens)

## ----marginal-figure, eval = has_ggplot2, echo = has_ggplot2, fig.height = 4.2, fig.cap = "The posterior of the slope from the proxy (orange, dashed) and from direct integration of the posterior on a grid (blue, solid). The grid curve does not use the proxy. Both curves are drawn where the grid density exceeds a thousandth of its peak.", fig.alt = "Two close density curves for the slope, one from the fitted Gaussian-mixture proxy and one from grid integration of the exact posterior."----
shown <- marg_quad > 1e-3 * max(marg_quad)
marg_df <- rbind(
  data.frame(slope = b_grid[shown], density = marg_quad[shown],
             source = "Grid integration of the posterior"),
  data.frame(slope = b_grid[shown], density = marg_proxy[shown],
             source = "Mixture proxy")
)

ggplot2::ggplot(marg_df,
                ggplot2::aes(slope, density, colour = source,
                             linetype = source)) +
  ggplot2::geom_line(linewidth = 0.8) +
  ggplot2::scale_colour_manual(
    name = NULL,
    values = c("Grid integration of the posterior" = "#0072B2",
               "Mixture proxy" = "#D55E00")
  ) +
  ggplot2::scale_linetype_manual(
    name = NULL,
    values = c("Grid integration of the posterior" = "solid",
               "Mixture proxy" = "dashed")
  ) +
  ggplot2::labs(
    title = "Posterior of the slope, two ways",
    x = expression(paste("slope ", beta, " (log-odds per 1000 lb)")),
    y = "posterior density"
  ) +
  ggplot2::theme_minimal(base_size = 11) +
  ggplot2::theme(plot.title = ggplot2::element_text(face = "bold"),
                 legend.position = "top")

## ----marginal-figure-skip, eval = !has_ggplot2, echo = FALSE, results = "asis"----
# cat("ggplot2 is not installed on this build, so the figure comparing the",
#     "two slope curves is skipped.\n")

## ----quadrature-tail----------------------------------------------------------
a_tail <- seq(-10 + 0.01, 10, by = 0.02)     # cell centres, intercept
b_tail <- seq(0.005, 3, by = 0.01)           # cell centres, slope > 0
log_tail <- vapply(b_tail, function(b) log_post(cbind(a_tail, b)),
                   numeric(length(a_tail)))
p_positive_quad <- sum(exp(log_tail - log_z_grid)) * 0.02 * 0.01
p_positive <- pgmm(0, slope, lower.tail = FALSE)
c(grid = signif(p_positive_quad, 2), proxy = signif(p_positive, 2))

## ----ensemble, warning = TRUE-------------------------------------------------
ens <- gmm_fit_ensemble(fit, B = 80L, seed = 3L)
ci_mean <- proxy_functional_ci(ens, gmm_mean, level = 0.9)
ci_tail <- proxy_functional_ci(
  ens, function(g) pgmm(0, gmm_marginalise(g, keep = 2L)), level = 0.9
)

## ----ensemble-table, echo = FALSE---------------------------------------------
## Three decimals for the two means, five for the probability.
ens_fmt <- function(v) {
  vapply(seq_along(v), function(i1) {
    formatC(v[i1], format = "f", digits = c(3L, 3L, 5L)[i1])
  }, character(1L))
}
knitr::kable(
  data.frame(
    quantity = c("posterior mean of the intercept",
                 "posterior mean of the slope",
                 "probability that the slope is negative"),
    estimate = ens_fmt(c(ci_mean$estimate, ci_tail$estimate)),
    conf_low = ens_fmt(c(ci_mean$conf.low, ci_tail$conf.low)),
    conf_high = ens_fmt(c(ci_mean$conf.high, ci_tail$conf.high))
  ),
  align = c("l", "r", "r", "r"),
  row.names = FALSE,
  col.names = c("Quantity", "Estimate", "Lower 5%", "Upper 95%"),
  caption = paste0("Bootstrap intervals over ", ens$B, " refits of the ",
                   "proxy.")
)

## ----compare-facts, include = FALSE-------------------------------------------
sim_value <- function(p, method, what) {
  res$sim_tab[[what]][res$sim_tab$p == p & res$sim_tab$method == method]
}
lead_z <- function(what) res$lead_tab$z[res$lead_tab$what == what]
lead_first <- function(what) {
  unique(res$lead_tab$first[res$lead_tab$what == what])
}
secs <- function(p, method) {
  format(signif(sim_value(p, method, "secs"), 2L), scientific = FALSE)
}
p1 <- res$p_vec[1L]
p2 <- res$p_vec[2L]
mcse_cov <- sqrt(0.95 * 0.05 / (res$n_rep * res$p_vec))
px_total <- sapply(res$p_vec, function(p) {
  sim_value(p, "proxymix", "secs") +
    res$secs_ensemble[[as.character(p)]]
})
nuts_secs <- sapply(res$p_vec, function(p) {
  sim_value(p, "NUTS, 4000 draws", "secs")
})
lap_ratio <- sapply(res$p_vec, function(p) {
  s1 <- res$sim_tab$p == p & res$sim_tab$method != "Laplace"
  min(res$sim_tab$secs[s1]) / sim_value(p, "Laplace", "secs")
})
stopifnot(lead_first("endpoints") == "proxymix",
          lead_first("err_log_z") == "proxymix",
          lead_first("err_mean") == "INLA",
          lead_first("err_p_pos") == "INLA")
## The sentences below on which differences are clear hold only if:
stopifnot(lead_z("err_p_pos")[2L] > 2, lead_z("err_p_pos")[1L] < 2,
          all(abs(res$z_px_nuts_mean) < 2),
          0.95 - res$coverage[[as.character(p1)]] > 2 * mcse_cov[1L])

## ----compare-table, echo = FALSE----------------------------------------------
cmp_tbl <- res$sim_tab
cmp_tbl$method[cmp_tbl$method == "NUTS, 4000 draws"] <- "Stan NUTS"
cmp_tbl$method[cmp_tbl$method == "DEzs"] <- "BayesianTools DEzs"
cmp_tbl$method[cmp_tbl$method == "Laplace"] <- "Laplace (LearnBayes)"
cmp_tbl$secs <- formatC(signif(cmp_tbl$secs, 2L), format = "fg",
                        digits = 2L)
old_na <- options(knitr.kable.NA = "--")
knitr::kable(
  cmp_tbl, digits = c(0L, 0L, 3L, 3L, 4L, 4L, 0L), row.names = FALSE,
  align = c("r", "l", "r", "r", "r", "r", "r"),
  col.names = c("$p$", "Method", "Mean", "Interval ends",
                "$P(\\beta_j > 0)$", "$\\log Z$", "Seconds"),
  caption = paste0(
    "Mean error against a 40,000-draw NUTS reference over ", res$n_rep,
    " simulated datasets for each number of coefficients $p$. Mean and ",
    "interval ends: error in reference posterior standard deviations, ",
    "averaged over the coefficients. $P(\\beta_j > 0)$: error in the ",
    "probability that a coefficient is positive. $\\log Z$: error in the ",
    "log normalising constant. Seconds: time per dataset on one core. The ",
    "NUTS time includes bridge sampling. DEzs gives no normalising constant."
  )
)
options(old_na)

## ----compare-code, eval = FALSE-----------------------------------------------
# library(proxymix)
# library(cmdstanr)
# library(bridgesampling)
# library(LearnBayes)
# library(INLA)
# library(BayesianTools)
# 
# # one simulated dataset: 200 observations, an intercept and one covariate
# n <- 200L
# p <- 2L
# prior_sd <- 10
# beta_pop <- c(0.5, 1.5, -1, 0.75, 0)   # the first p entries are used
# r <- 1L
# set.seed(r)
# X <- cbind(1, matrix(rnorm(n * (p - 1L)), n))
# y <- rbinom(n, 1L, plogis(X %*% beta_pop[seq_len(p)]))
# data <- list(n = n, p = p, X = X, y = y, prior_sd = prior_sd)
# 
# log_lik_one <- function(theta, data) {
#   eta <- data$X %*% theta
#   sum(data$y * eta - log1p(exp(eta)))
# }
# log_post_one <- function(theta, data) {
#   log_lik_one(theta, data) + sum(dnorm(theta, 0, data$prior_sd, log = TRUE))
# }
# log_post_rows <- function(theta, data) {
#   if (is.null(dim(theta))) theta <- matrix(theta, ncol = data$p)
#   apply(theta, 1L, log_post_one, data = data)
# }
# mle <- glm(y ~ X - 1, family = binomial())
# 
# # proxymix: fit, log normalising constant, bootstrap intervals
# tgt <- gmm_target(n_dim = p,
#                   log_density = function(theta) log_post_rows(theta, data),
#                   normalised = FALSE)
# q0 <- proposal_mvt(p, mean = coef(mle), sigma = 9 * vcov(mle), df = 5)
# fit <- fit_proxymix(tgt, N = 2L, regime = "kld", proposal = q0,
#                     is_size = 3000L, max_iter = 60L, seed = r, adapt = "pmc")
# ev <- gmm_evidence(fit, n = 4000L, seed = r)
# ens <- gmm_fit_ensemble(fit, B = 80L, seed = r)
# ci <- proxy_functional_ci(ens, gmm_mean, level = 0.95)
# 
# # Stan NUTS: four chains of 1000 warm-up and 1000 kept draws, then bridge
# # sampling on the draws for the log normalising constant
# stan_file <- file.path(tempdir(), "logistic.stan")
# writeLines(c(
#   "data {",
#   "  int<lower=1> n;",
#   "  int<lower=1> p;",
#   "  matrix[n, p] X;",
#   "  array[n] int<lower=0, upper=1> y;",
#   "  real<lower=0> prior_sd;",
#   "}",
#   "parameters {",
#   "  vector[p] beta;",
#   "}",
#   "model {",
#   "  beta ~ normal(0, prior_sd);",
#   "  y ~ bernoulli_logit(X * beta);",
#   "}"
# ), stan_file)
# model <- cmdstan_model(stan_file)
# fit_nuts <- model$sample(data, chains = 4L, parallel_chains = 1L,
#                          iter_warmup = 1000L, iter_sampling = 1000L,
#                          seed = r + 100000L, refresh = 0L,
#                          show_messages = FALSE, show_exceptions = FALSE)
# draws <- fit_nuts$draws("beta", format = "matrix")
# draws <- matrix(draws, ncol = data$p, dimnames = list(NULL, colnames(draws)))
# bounds <- setNames(rep(-Inf, data$p), colnames(draws))
# bridge <- bridge_sampler(draws, log_posterior = log_post_one, data = data,
#                          lb = bounds, ub = -bounds, silent = TRUE)
# 
# # Laplace approximation
# lap <- laplace(log_post_one, coef(mle), data)
# 
# # INLA
# df <- data.frame(y = y, X[, -1L, drop = FALSE])
# names(df) <- c("y", paste0("x", seq_len(p - 1L)))
# fi <- inla(reformulate(names(df)[-1L], "y"), family = "binomial",
#            Ntrials = 1, data = df,
#            control.fixed = list(prec = 1 / prior_sd^2,
#                                 prec.intercept = 1 / prior_sd^2),
#            num.threads = "1:1")
# 
# # BayesianTools, differential-evolution sampler DEzs
# setup <- createBayesianSetup(
#   likelihood = function(theta) log_lik_one(theta, data),
#   prior = createPrior(
#     density = function(theta) sum(dnorm(theta, 0, prior_sd, log = TRUE)),
#     sampler = function(n = 1L) matrix(rnorm(n * p, 0, prior_sd), n)
#   )
# )
# set.seed(r)
# bt <- runMCMC(setup, sampler = "DEzs",
#               settings = list(iterations = 15000L, message = FALSE))
# draws_bt <- getSample(bt, start = 1000L)

## ----session-info, collapse = FALSE, class.output = "session-info"------------
sessionInfo()

