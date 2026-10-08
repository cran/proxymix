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
res <- readRDS("results/from_kde.rds")
major_minor <- function(v) paste(unlist(package_version(v))[1:2],
                                 collapse = ".")
if (major_minor(res$proxymix_version) !=
    major_minor(as.character(packageVersion("proxymix")))) {
  stop("results/from_kde.rds was built under proxymix ",
       res$proxymix_version, ", but this is proxymix ",
       packageVersion("proxymix"), ". Rerun the simulation and ",
       "data-raw/vignette_results/from_kde.R.", call. = FALSE)
}

## Small numbers are written as plain decimals rather than in the
## scientific notation that knitr's inline hook would otherwise use.
fixed <- function(v, digits) {
  format(round(v, digits), nsmall = digits, scientific = FALSE)
}

## ----recover-bimodal----------------------------------------------------------
set.seed(20260601)

n_each <- 150L
true_means <- cbind(c(-2, 0), c(2, 0))
true_cov <- diag(2)
x <- rbind(
  mvnfast::rmvn(n_each, mu = true_means[, 1L], sigma = true_cov),
  mvnfast::rmvn(n_each, mu = true_means[, 2L], sigma = true_cov)
)

fit <- from_kde(
  x, N = 2L,
  bandwidth = "silverman",
  is_size = 2000L, max_iter = 60L, seed = 1L,
  validation_size = 2000L
)
fit

## ----recovered-means----------------------------------------------------------
mu_hat <- vapply(fit@means, function(mu) mu, numeric(2L))
comp_order <- order(mu_hat[1L, ])
mu_hat <- mu_hat[, comp_order, drop = FALSE]
weight_hat <- fit@weights[comp_order]
mean_error <- max(abs(mu_hat - true_means))

## ----recovered-means-table, echo = FALSE--------------------------------------
knitr::kable(
  data.frame(
    component = c(1L, 2L),
    fitted_x1 = round(mu_hat[1L, ], 3),
    fitted_x2 = round(mu_hat[2L, ], 3),
    true_x1 = true_means[1L, ],
    true_x2 = true_means[2L, ],
    weight = round(weight_hat, 3)
  ),
  row.names = FALSE,
  col.names = c("Component", "Fitted $x_1$", "Fitted $x_2$", "True $x_1$",
                "True $x_2$", "Weight"),
  caption = paste0(
    "The means and weights of the two fitted components beside the means ",
    "of the two groups the ", nrow(x), " points were drawn from."
  )
)

## ----fit-quality--------------------------------------------------------------
es <- ess_summary(fit)
print(data.frame(is_size = es$is_size, ess = round(es$ess, 1),
                 ess_relative = round(es$ess_relative, 3),
                 max_weight = signif(es$max_weight, 3)),
      row.names = FALSE)
print(data.frame(validation_size = es$validation_size,
                 validation_kld = signif(es$validation_kld, 3),
                 validation_se = signif(fit@diagnostics$validation_mc_se, 3)),
      row.names = FALSE)

## ----compression-cost---------------------------------------------------------
g1 <- seq(-6, 6, length.out = 160L)
g2 <- seq(-5, 5, length.out = 140L)
grid <- expand.grid(x1 = g1, x2 = g2)
gm <- as.matrix(grid)
cell <- (g1[2L] - g1[1L]) * (g2[2L] - g2[1L])

grid$kde <- exp(fit@target@log_density(gm))
grid$proxy <- dgmm(gm, fit)

total_variation <- 0.5 * sum(abs(grid$kde - grid$proxy)) * cell
hell <- hellinger_mc(fit, n_mc = 10000L, seed = 1L)
c(total_variation = signif(total_variation, 3),
  hellinger_sq = signif(hell$h2, 3), hellinger_se = signif(hell$se, 3))

## ----visualise, eval = has_ggplot2, echo = has_ggplot2, fig.height = 4.5, fig.cap = "Contours of the log-density of the kernel estimate (blue, solid) and of the two-component proxy (orange, dashed), on shared levels, over the 300 data points (grey). Around the two group centres the contours nearly coincide. On the outer, low-density levels the kernel estimate bends around single outlying points, and the proxy stays smooth.", fig.alt = "Contour plot of the kernel-density log-density and the Gaussian-mixture proxy log-density on a planar grid, with sample points overlaid. The inner contours around the two group centres nearly coincide; the outer contours of the kernel estimate are wavy and those of the proxy are smooth ellipses."----
plot_grid <- expand.grid(
  x1 = seq(-5, 5, length.out = 80L),
  x2 = seq(-4, 4, length.out = 60L)
)
pm <- as.matrix(plot_grid)
plot_grid$kde <- fit@target@log_density(pm)
plot_grid$proxy <- log(dgmm(pm, fit))

## Shared contour levels so the two log-densities are directly comparable.
brks <- pretty(range(c(plot_grid$kde, plot_grid$proxy), finite = TRUE), 9L)
sample_df <- data.frame(x1 = x[, 1L], x2 = x[, 2L])

ggplot2::ggplot() +
  ggplot2::geom_point(data = sample_df, ggplot2::aes(x1, x2),
                      colour = "grey60", alpha = 0.25, size = 0.5) +
  ggplot2::geom_contour(
    data = plot_grid,
    ggplot2::aes(x1, x2, z = kde, colour = "Kernel estimate",
                 linetype = "Kernel estimate"),
    breaks = brks, linewidth = 0.45
  ) +
  ggplot2::geom_contour(
    data = plot_grid,
    ggplot2::aes(x1, x2, z = proxy, colour = "Mixture proxy",
                 linetype = "Mixture proxy"),
    breaks = brks, linewidth = 0.55
  ) +
  ggplot2::scale_colour_manual(
    name = NULL,
    values = c("Kernel estimate" = "#0072B2", "Mixture proxy" = "#D55E00")
  ) +
  ggplot2::scale_linetype_manual(
    name = NULL,
    values = c("Kernel estimate" = "solid", "Mixture proxy" = "dashed")
  ) +
  ggplot2::coord_equal(expand = FALSE) +
  ggplot2::labs(
    title = "Kernel estimate and two-component proxy",
    subtitle = "log-density contours on shared levels",
    x = expression(x[1]), y = expression(x[2])
  ) +
  ggplot2::theme_minimal(base_size = 11) +
  ggplot2::theme(plot.title = ggplot2::element_text(face = "bold"),
                 legend.position = "top",
                 panel.grid.minor = ggplot2::element_blank())

## ----visualise-skip, eval = !has_ggplot2, echo = FALSE, results = "asis"------
# cat("ggplot2 is not installed on this build, so the contour-comparison",
#     "figure is skipped.\n")

## ----compose------------------------------------------------------------------
slice <- gmm_conditionalise(fit, given = c(NA, 0))
draws <- rgmm(200L, slice)
c(components = gmm_n_components(slice), dimension = gmm_dim(slice),
  draws = nrow(draws))

## ----bandwidth-sweep----------------------------------------------------------
bandwidth_grid <- c(0.2, 0.5, 1.0)
fits <- lapply(bandwidth_grid, function(h) {
  from_kde(x, N = 2L, bandwidth = h,
           is_size = 1500L, max_iter = 40L, seed = 1L)
})
sweep_is <- fits[[1L]]@diagnostics$is_size
sweep_vs <- vapply(fits, function(f) f@diagnostics$validation_size,
                   numeric(1L))
sweep_val <- vapply(fits, function(f) f@diagnostics$validation_kld,
                    numeric(1L))
sweep_se <- vapply(fits, function(f) f@diagnostics$validation_mc_se,
                   numeric(1L))

## ----bandwidth-check, include = FALSE-----------------------------------------
## spread of the left-hand component: the stored order of components can
## differ between fits
trace_left <- function(f) {
  j <- which.min(vapply(f@means, function(mu) mu[1L], numeric(1L)))
  sum(diag(f@covariances[[j]]))
}
## the prose below states that the KL divergence falls as the bandwidth
## widens and that at the widest it is more than two standard errors above zero
if (is.unsorted(rev(sweep_val)) || sweep_val[3L] <= 2 * sweep_se[3L]) {
  stop("The text on the bandwidth sweep no longer matches the fits.",
       call. = FALSE)
}

## ----bandwidth-table, echo = FALSE--------------------------------------------
knitr::kable(
  data.frame(
    bandwidth = bandwidth_grid,
    ess = round(vapply(fits, function(f) f@diagnostics$ess, numeric(1L)), 1),
    max_weight = signif(
      vapply(fits, function(f) f@diagnostics$max_weight, numeric(1L)), 3
    ),
    trace_sigma = round(vapply(fits, trace_left, numeric(1L)), 3),
    fresh = sweep_vs,
    val_kl = vapply(sweep_val, function(v) {
      formatC(signif(v, 2), digits = 2, format = "fg", flag = "#")
    }, character(1L)),
    val_se = vapply(sweep_se, function(v) {
      format(signif(v, 2), scientific = FALSE)
    }, character(1L))
  ),
  row.names = FALSE, align = "r",
  col.names = c("Bandwidth", "Effective sample size", "Largest weight",
                "Spread of left component", "Fresh draws", "KL, fresh draws",
                "Standard error"),
  caption = paste0(
    "The same sample compressed with three bandwidths, at ", sweep_is,
    " weighted draws. The spread of ",
    "the left-hand component is the sum of its two variances."
  )
)

## ----compare-facts, include = FALSE-------------------------------------------
sim_value <- function(method, what) {
  res$sim_tab[[what]][res$sim_tab$method == method]
}
paired_value <- function(method, what) {
  res$paired[[what]][res$paired$method == method]
}
ise_k <- function(K) {
  fixed(1000 * sim_value(paste0("proxymix, K = ", K), "ise"), 2)
}
query_ratio_ks <- res$query_ms[["ks, kde"]] /
  res$query_ms[["proxymix, K = 3"]]
query_ratio_bk <- res$query_ms[["proxymix, K = 3"]] /
  res$query_ms[["KernSmooth, bkde2D"]]
## the prose below reads these orderings from the stored results
kernel_ise <- 1000 * c(sim_value("ks, kde", "ise"),
                       sim_value("np, npudens", "ise"),
                       sim_value("KernSmooth, bkde2D", "ise"))
proxy_ise <- 1000 * vapply(c(3L, 5L, 8L), function(K) {
  sim_value(paste0("proxymix, K = ", K), "ise")
}, numeric(1L))
mixture_ise <- 1000 * c(sim_value("mclust, G by BIC", "ise"),
                        sim_value("mixtools, k = 3", "ise"))
kernel_secs <- c(sim_value("ks, kde", "secs"),
                 sim_value("np, npudens", "secs"),
                 sim_value("KernSmooth, bkde2D", "secs"))
if (!(max(proxy_ise[1:2]) < min(kernel_ise) &&
      max(mixture_ise) < min(proxy_ise) &&
      sim_value("proxymix, K = 3", "secs") > max(kernel_secs) &&
      query_ratio_ks > 1 && query_ratio_bk > 1)) {
  stop("The text on the comparison no longer matches ",
       "results/from_kde.rds.", call. = FALSE)
}

## ----compare-table, echo = FALSE----------------------------------------------
method_order <- c("ks, kde", "KernSmooth, bkde2D", "np, npudens",
                  "proxymix, K = 3", "proxymix, K = 5", "proxymix, K = 8",
                  "mclust, G by BIC", "mixtools, k = 3")
method_label <- c("ks", "KernSmooth", "np", "proxymix, 3 components",
                  "proxymix, 5 components", "proxymix, 8 components",
                  "mclust, components by BIC", "mixtools, 3 components")
cmp_tbl <- data.frame(
  method = method_label,
  ise = 1000 * vapply(method_order, sim_value, numeric(1L), what = "ise"),
  ise_se = 1000 * vapply(method_order, sim_value, numeric(1L),
                         what = "ise_se"),
  diff = c(NA, vapply(method_order[-1L], paired_value, numeric(1L),
                      what = "mean")),
  diff_se = c(NA, vapply(method_order[-1L], paired_value, numeric(1L),
                         what = "se")),
  secs = vapply(method_order, sim_value, numeric(1L), what = "secs"),
  stringsAsFactors = FALSE
)
cmp_tbl$ise <- sprintf("%.3f (%.3f)", cmp_tbl$ise, cmp_tbl$ise_se)
cmp_tbl$diff <- ifelse(is.na(cmp_tbl$diff), "",
                       sprintf("%.4f (%.4f)", cmp_tbl$diff, cmp_tbl$diff_se))
cmp_tbl$secs <- sprintf("%.3f", cmp_tbl$secs)
knitr::kable(
  cmp_tbl[, c("method", "ise", "diff", "secs")],
  row.names = FALSE, align = c("l", "r", "r", "r"),
  col.names = c("Method", "Error", "Difference from ks",
                "Seconds per fit"),
  caption = paste0(
    "Integrated squared error against the true density, in thousandths, ",
    "averaged over ", res$n_rep, " simulated datasets of ", res$n,
    " points, with its standard error in brackets. The difference from ks ",
    "is taken dataset by dataset. A negative value means a smaller error ",
    "than ks. Seconds per fit is the mean time to fit one dataset."
  )
)

## ----fit-time-check, include = FALSE------------------------------------------
## the sentence below orders the fit times
px3 <- sim_value("proxymix, K = 3", "secs")
stopifnot(px3 > sim_value("mclust, G by BIC", "secs"),
          px3 > sim_value("ks, kde", "secs"),
          sim_value("proxymix, K = 8", "secs") < sim_value("mixtools, k = 3", "secs"))

## ----compare-code, eval = FALSE-----------------------------------------------
# library(proxymix)
# library(mclust)
# library(mixtools)
# library(np)
# options(np.messages = FALSE)
# 
# # one dataset of 500 points from a mixture of three normal distributions
# truth <- list(
#   weights = c(0.35, 0.35, 0.30),
#   means = list(c(-2, 0), c(2, 1), c(0, 4)),
#   covs = list(diag(2), matrix(c(1, 0.5, 0.5, 1), 2L), 0.6 * diag(2))
# )
# set.seed(1L)
# n <- 500L
# k <- sample.int(3L, n, replace = TRUE, prob = truth$weights)
# x <- matrix(NA_real_, n, 2L)
# for (j in seq_len(3L)) {
#   s <- k == j
#   x[s, ] <- mvnfast::rmvn(sum(s), truth$means[[j]], truth$covs[[j]])
# }
# 
# # integrated squared error against the true density, on a grid
# g1 <- seq(-6, 6, length.out = 121L)
# g2 <- seq(-4, 8, length.out = 121L)
# grid <- as.matrix(expand.grid(x1 = g1, x2 = g2))
# cell <- (g1[2L] - g1[1L]) * (g2[2L] - g2[1L])
# f_true <- Reduce(`+`, lapply(seq_len(3L), function(j) {
#   truth$weights[j] * mvnfast::dmvn(grid, truth$means[[j]], truth$covs[[j]])
# }))
# ise <- function(f_hat) sum((f_hat - f_true)^2) * cell
# as_gmm <- function(w, mu, sigma) {
#   gmm(weights = w, means = mu, covariances = sigma)
# }
# 
# # ks with its diagonal plug-in bandwidth, and proxymix compressing the same estimate
# h <- sqrt(diag(ks::Hpi.diag(x)))
# kd <- ks::kde(x, H = diag(h^2), eval.points = grid)
# f <- from_kde(x, N = 3L, bandwidth = h, is_size = 20000L, seed = 1L)
# 
# mc <- Mclust(x, verbose = FALSE)
# g_mc <- as_gmm(
#   mc$parameters$pro,
#   lapply(seq_len(mc$G), function(j) mc$parameters$mean[, j]),
#   lapply(seq_len(mc$G), function(j) mc$parameters$variance$sigma[, , j])
# )
# mt <- mvnormalmixEM(x, k = 3L)
# nd <- npudens(npudensbw(x))
# f_np <- predict(nd, newdata = data.frame(grid))
# bk <- KernSmooth::bkde2D(x, bandwidth = h, gridsize = c(121L, 121L),
#                          range.x = list(range(g1), range(g2)))
# 
# 1000 * c(ks = ise(kd$estimate),
#          proxymix = ise(dgmm(grid, f)),
#          mclust = ise(dgmm(grid, g_mc)),
#          mixtools = ise(dgmm(grid, as_gmm(mt$lambda, mt$mu, mt$sigma))),
#          np = ise(f_np),
#          KernSmooth = ise(as.vector(bk$fhat)))

## ----session-info, collapse = FALSE, class.output = "session-info"------------
sessionInfo()

