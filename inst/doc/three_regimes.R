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
res <- readRDS("results/three_regimes.rds")
major_minor <- function(v) paste(unlist(package_version(v))[1:2],
                                 collapse = ".")
if (major_minor(res$proxymix_version) !=
    major_minor(as.character(packageVersion("proxymix")))) {
  stop("results/three_regimes.rds was built under proxymix ",
       res$proxymix_version, ", but this is proxymix ",
       packageVersion("proxymix"), ". Rerun the simulation and ",
       "data-raw/vignette_results/three_regimes.R.", call. = FALSE)
}

## Small numbers are written as plain decimals rather than in the
## scientific notation that knitr's inline hook would otherwise use.
fixed <- function(v, digits) {
  format(round(v, digits), nsmall = digits, scientific = FALSE)
}
## A very small number as a power of ten in LaTeX, to two significant digits.
power_ten <- function(v) {
  if (v == 0) return("$0$")
  e <- floor(log10(abs(v)))
  m <- signif(v / 10^e, 2L)
  if (m == 1) sprintf("$10^{%d}$", e) else sprintf("$%s \\times 10^{%d}$", format(m), e)
}

## ----target-------------------------------------------------------------------
tgt <- mixture_target(with_samples = TRUE, n = 1500L, seed = 1L)
tgt

## ----moment-------------------------------------------------------------------
m_fit <- fit_proxymix(tgt, N = 1L, regime = "moment")
m_fit

## ----moment-recover-----------------------------------------------------------
m_fit_bare <- fit_proxymix(tgt, N = 1L, regime = "moment", ridge_eps = 0)
moment_gap <- c(
  mean = max(abs(m_fit@means[[1L]] - colMeans(tgt@samples))),
  covariance = max(abs(m_fit@covariances[[1L]] - cov(tgt@samples))),
  covariance_no_ridge = max(abs(m_fit_bare@covariances[[1L]] -
                                  cov(tgt@samples)))
)
signif(moment_gap, 3L)

## ----em-----------------------------------------------------------------------
s_fit <- fit_proxymix(tgt, N = 3L, regime = "sample",
                      max_iter = 200L, n_starts = 4L, seed = 1L)
s_fit

## ----kld----------------------------------------------------------------------
k_fit <- fit_proxymix(tgt, N = 3L, regime = "kld",
                      proposal = proposal_mvt(n_dim = 2L,
                                              mean = c(0, 0),
                                              sigma = 6 * diag(2),
                                              df = 5),
                      is_size = 3000L,
                      max_iter = 60L,
                      seed = 1L)
k_fit

## ----overlay-grid-------------------------------------------------------------
grid_x <- seq(-4.5, 4.5, length.out = 100L)
grid_base <- expand.grid(x1 = grid_x, x2 = grid_x)
grid_mat <- as.matrix(grid_base)
target_d <- exp(tgt@log_density(grid_mat))
panel_of <- function(fit, label) {
  data.frame(
    x1 = grid_base$x1, x2 = grid_base$x2,
    target = target_d,
    proxy = dgmm(grid_mat, fit),
    regime = label,
    stringsAsFactors = FALSE
  )
}
overlay_df <- rbind(
  panel_of(m_fit, "(i) moment, N = 1"),
  panel_of(s_fit, "(ii) sample EM, N = 3"),
  panel_of(k_fit, "(iii) KLD-EM, N = 3")
)
overlay_df$regime <- factor(overlay_df$regime,
                            levels = unique(overlay_df$regime))

## ----overlay, eval = has_ggplot2, echo = has_ggplot2, fig.height = 3.4, fig.cap = "The three-peak target (filled contours, the same in all three panels) with each method's proxy overlaid as dashed contours. Method (i) spreads one normal distribution across all three peaks. Methods (ii) and (iii) each place one component on each peak and look alike, although (ii) used only the sample and (iii) only the formula.", fig.alt = "Three side-by-side contour panels of the same three-peak target, overlaid with the single-normal moment fit, the sample-EM fit and the KLD-EM fit."----
ggplot2::ggplot(overlay_df, ggplot2::aes(x1, x2)) +
  ggplot2::geom_contour_filled(ggplot2::aes(z = target), bins = 10L,
                               alpha = 0.85) +
  ggplot2::geom_contour(ggplot2::aes(z = proxy), colour = "white",
                        linetype = "dashed", linewidth = 0.4, bins = 5L) +
  ggplot2::scale_fill_viridis_d(option = "mako", guide = "none") +
  ggplot2::facet_wrap(~ regime) +
  ggplot2::coord_equal() +
  ggplot2::labs(x = expression(x[1]), y = expression(x[2])) +
  ggplot2::theme_minimal(base_size = 11)

## ----overlay-skip, eval = !has_ggplot2, echo = FALSE, results = "asis"--------
# cat("ggplot2 is not installed on this build, so the three-panel overlay",
#     "figure is skipped.\n")

## ----traces-------------------------------------------------------------------
trace_df <- rbind(
  data.frame(
    iteration = seq_along(s_fit@diagnostics$loglik_trace),
    value = s_fit@diagnostics$loglik_trace,
    panel = "(ii) sample EM: log-likelihood (up)",
    stringsAsFactors = FALSE
  ),
  data.frame(
    iteration = seq_along(kld_trace(k_fit)),
    value = kld_trace(k_fit),
    panel = "(iii) KLD-EM: KL on the fitting draws (down)",
    stringsAsFactors = FALSE
  )
)

## ----traces-plot, eval = has_ggplot2, echo = has_ggplot2, fig.height = 3.2, fig.cap = "The value each method improves, round by round. Method (ii) raises the log-likelihood of the sample. Method (iii) lowers an estimate of the KL divergence, scored on the same trial points the fit was tuned to. That makes the estimate read low, and it can fall below zero, although a true KL divergence cannot.", fig.alt = "Two panels of iteration traces: an increasing log-likelihood curve for sample EM and a decreasing Kullback-Leibler curve for KLD-EM."----
ggplot2::ggplot(trace_df, ggplot2::aes(iteration, value)) +
  ggplot2::geom_line(colour = "#0072B2", linewidth = 0.8) +
  ggplot2::geom_point(colour = "#0072B2", size = 1.1) +
  ggplot2::facet_wrap(~ panel, scales = "free") +
  ggplot2::scale_x_continuous(breaks = function(lim) unique(floor(pretty(lim)))) +
  ggplot2::labs(x = "round", y = "value") +
  ggplot2::theme_minimal(base_size = 11)

## ----traces-skip, eval = !has_ggplot2, echo = FALSE, results = "asis"---------
# cat("ggplot2 is not installed on this build, so the two-panel iteration",
#     "trace figure is skipped.\n")

## ----n1-fits------------------------------------------------------------------
banana <- banana_target(with_samples = TRUE, n = 2000L, seed = 1L)
m_b <- fit_proxymix(banana, N = 1L, regime = "moment")
k_b <- fit_proxymix(banana, N = 1L, regime = "kld",
                    proposal = proposal_mvt(n_dim = 2L,
                                            sigma = 4 * diag(2),
                                            df = 5),
                    is_size = 3000L, max_iter = 50L, seed = 1L)
tr_of <- function(f) sum(diag(f@covariances[[1L]]))
sample_trace <- sum(diag(cov(banana@samples)))
exact_trace <- 1 + (1 + 2 * 0.5^2)

## ----n1-table, echo = FALSE---------------------------------------------------
n1_traces <- c(tr_of(m_b), tr_of(k_b), sample_trace, exact_trace)
n1_tbl <- data.frame(
  Source = c("method (i), moment match", "method (iii), KLD-EM",
             "the attached sample", "exact value"),
  Uses = c("the sample", "the formula", "the sample",
           "the construction of the target"),
  `Trace of covariance` = round(n1_traces, 3L),
  `Error` = round(n1_traces - exact_trace, 3L),
  check.names = FALSE,
  stringsAsFactors = FALSE
)
knitr::kable(
  n1_tbl,
  caption = paste(
    "The single normal proxy for the banana target, fitted two ways. The",
    "error is the trace of the covariance minus its exact value, 5/2, so a",
    "negative error means the spread is understated."
  )
)

## ----n1-spread----------------------------------------------------------------
sample_traces <- vapply(seq_len(500L), function(s) {
  sum(diag(cov(banana_target(with_samples = TRUE, n = 2000L,
                             seed = s)@samples)))
}, numeric(1L))
kld_traces <- vapply(seq_len(100L), function(s) {
  tr_of(fit_proxymix(banana, N = 1L, regime = "kld",
                     proposal = proposal_mvt(n_dim = 2L,
                                             sigma = 4 * diag(2),
                                             df = 5),
                     is_size = 3000L, max_iter = 50L, seed = s))
}, numeric(1L))
spread <- data.frame(
  seeds = c(length(sample_traces), length(kld_traces)),
  mean = c(mean(sample_traces), mean(kld_traces)),
  sd = c(sd(sample_traces), sd(kld_traces)),
  row.names = c("trace of a 2000-point sample", "method (iii) trace")
)
round(spread, 3L)

## ----n1-grid, include = FALSE-------------------------------------------------
# a grid sum of the same trace, for the Limitations section
quad_x <- seq(-8, 8, length.out = 400L)
quad_g <- as.matrix(expand.grid(x1 = quad_x, x2 = quad_x))
quad_cell <- (quad_x[2L] - quad_x[1L])^2
quad_f <- exp(banana@log_density(quad_g))
quad_mass <- sum(quad_f) * quad_cell
quad_trace <- sum(quad_f * (quad_g[, 1L]^2 + quad_g[, 2L]^2)) *
  quad_cell / quad_mass

## ----compare-facts, include = FALSE-------------------------------------------
sim_value <- function(n, method, what) {
  s1 <- res$sim_tab$n == n & res$sim_tab$method == method
  res$sim_tab[[what]][s1]
}
misfits <- function(n, method) {
  s1 <- res$misfit_tab$n == n & res$misfit_tab$method == method
  res$misfit_tab$misfit[s1]
}
n_small <- res$n_sizes[1L]
n_large <- res$n_sizes[2L]
ii <- "proxymix, regime (ii)"
ll <- res$faithful_loglik

## ----compare-table, echo = FALSE----------------------------------------------
pkgs <- c(ii, "mclust", "mixtools", "flexmix")
sim_col <- function(n, what, digits) {
  vapply(pkgs, function(s1) {
    if (s1 == "mixtools") return("not run")
    if (what == "misfit") return(as.character(misfits(n, s1)))
    fixed(sim_value(n, s1, what), digits)
  }, character(1L))
}
cmp_tbl <- data.frame(
  method = c("proxymix, method (ii)", "mclust", "mixtools", "flexmix"),
  loglik = fixed(ll[c("proxymix", "mclust", "mixtools", "flexmix")], 3),
  kl_small = sim_col(n_small, "kl", 4),
  kl_large = sim_col(n_large, "kl", 4),
  mis_small = sim_col(n_small, "misfit", 0),
  mis_large = sim_col(n_large, "misfit", 0),
  stringsAsFactors = FALSE
)
knitr::kable(
  cmp_tbl, row.names = FALSE,
  align = c("l", "r", "r", "r", "r", "r"),
  col.names = c("Package", "Old Faithful: held-out log-likelihood",
                paste0("Mean KL, ", n_small, " points"),
                paste0("Mean KL, ", format(n_large, big.mark = ","),
                       " points"),
                paste0("Misfits, ", n_small, " points"),
                paste0("Misfits, ", format(n_large, big.mark = ","),
                       " points")),
  caption = paste0(
    "Two-component fits to half of the Old Faithful eruptions, scored on ",
    "the other half (higher is better), and three-component fits to ",
    res$n_rep, " simulated datasets at each sample size, scored by the KL ",
    "divergence from the true mixture (lower is better). A misfit is a fit ",
    "with a divergence above 0.1, counted out of ", res$n_rep, "."
  )
)

## ----compare-code, eval = FALSE-----------------------------------------------
# library(proxymix)
# library(mclust)
# library(mixtools)
# library(flexmix)
# 
# # split the Old Faithful eruptions into a training half and a held-out half
# faithful_mat <- as.matrix(datasets::faithful)
# set.seed(20260925)
# i_train <- sort(sample.int(nrow(faithful_mat), nrow(faithful_mat) / 2L))
# train <- faithful_mat[i_train, ]
# test <- faithful_mat[-i_train, ]
# train_df <- data.frame(eruptions = train[, 1L], waiting = train[, 2L])
# 
# # convert each package's fit to a proxymix mixture, so dgmm() scores all four
# as_gmm_mclust <- function(fit) {
#   gmm(
#     weights = fit$parameters$pro,
#     means = lapply(seq_len(fit$G), function(k) fit$parameters$mean[, k]),
#     covariances = lapply(seq_len(fit$G), function(k) {
#       fit$parameters$variance$sigma[, , k]
#     })
#   )
# }
# as_gmm_mixtools <- function(fit) {
#   gmm(weights = fit$lambda, means = fit$mu, covariances = fit$sigma)
# }
# as_gmm_flexmix <- function(fit) {
#   comps <- lapply(fit@components, function(cc) cc[[1L]]@parameters)
#   gmm(
#     weights = prior(fit),
#     means = lapply(comps, function(p) unname(p$center)),
#     covariances = lapply(comps, function(p) unname(p$cov))
#   )
# }
# 
# # two-component fits to the training half
# set.seed(1L)
# fits <- list(
#   proxymix = fit_proxymix(gmm_target_from_samples(train), N = 2L,
#                           regime = "sample"),
#   mclust = as_gmm_mclust(Mclust(train, G = 2L, verbose = FALSE)),
#   mixtools = as_gmm_mixtools(mvnormalmixEM(train, k = 2L, verb = FALSE)),
#   flexmix = as_gmm_flexmix(stepFlexmix(
#     cbind(eruptions, waiting) ~ 1, data = train_df, k = 2L, nrep = 5L,
#     model = FLXMCmvnorm(diagonal = FALSE), verbose = FALSE
#   ))
# )
# 
# # mean log-likelihood of the held-out eruptions under each fit
# vapply(fits, function(g) mean(dgmm(test, g, log = TRUE)), numeric(1L))

## ----session-info, collapse = FALSE, class.output = "session-info"------------
sessionInfo()

