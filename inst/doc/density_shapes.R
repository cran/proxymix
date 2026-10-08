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

## ----fits---------------------------------------------------------------------
fit_b <- fit_proxymix(banana_target(), N = 4L, regime = "kld",
                      proposal = proposal_mvt(n_dim = 2L,
                                              sigma = 4 * diag(2),
                                              df = 5),
                      is_size = 3000L, max_iter = 50L, seed = 1L)
fit_d <- fit_proxymix(donut_target(), N = 6L, regime = "kld",
                      proposal = proposal_mvt(n_dim = 2L,
                                              sigma = 9 * diag(2),
                                              df = 5),
                      is_size = 3500L, max_iter = 60L, seed = 1L)
fit_m <- fit_proxymix(mixture_target(), N = 3L, regime = "kld",
                      proposal = proposal_mvt(n_dim = 2L,
                                              sigma = 6 * diag(2),
                                              df = 5),
                      is_size = 3000L, max_iter = 50L, seed = 1L)

## ----overlay-grid-------------------------------------------------------------
shape_panel <- function(target, fit, label, lim) {
  gx <- seq(-lim, lim, length.out = 110L)
  base <- expand.grid(x1 = gx, x2 = gx)
  gm <- as.matrix(base)
  data.frame(
    x1 = base$x1, x2 = base$x2,
    target = exp(target@log_density(gm)),
    proxy = dgmm(gm, fit),
    shape = label,
    stringsAsFactors = FALSE
  )
}
overlay_df <- rbind(
  shape_panel(banana_target(), fit_b, "banana, N = 4", 4.5),
  shape_panel(donut_target(), fit_d, "donut, N = 6", 4.5),
  shape_panel(mixture_target(), fit_m, "three peaks, N = 3", 4.5)
)
overlay_df$shape <- factor(overlay_df$shape,
                           levels = unique(overlay_df$shape))

## ----overlay, eval = has_ggplot2, echo = has_ggplot2, fig.height = 3.4, fig.cap = "Each target (filled contours) with its fitted proxy overlaid as dashed contours. The proxies for the banana and the three peaks lie close to their targets. The donut's proxy is a ring of six ellipses, and its dashed contours break into separate lobes where neighbouring ellipses meet.", fig.alt = "Three contour panels showing the banana, donut and three-peak densities, each overlaid with dashed contours of its Gaussian-mixture proxy."----
ggplot2::ggplot(overlay_df, ggplot2::aes(x1, x2)) +
  ggplot2::geom_contour_filled(ggplot2::aes(z = target), bins = 10L,
                               alpha = 0.85) +
  ggplot2::geom_contour(ggplot2::aes(z = proxy), colour = "white",
                        linetype = "dashed", linewidth = 0.4, bins = 8L) +
  ggplot2::scale_fill_viridis_d(option = "mako", guide = "none") +
  ggplot2::facet_wrap(~ shape) +
  ggplot2::coord_equal() +
  ggplot2::labs(x = expression(x[1]), y = expression(x[2])) +
  ggplot2::theme_minimal(base_size = 11)

## ----overlay-skip, eval = !has_ggplot2, echo = FALSE, results = "asis"--------
# cat("ggplot2 is not installed on this build, so the three-shape overlay",
#     "figure is skipped.\n")

## ----kl-grid------------------------------------------------------------------
grid_points <- function(lim, step = 0.04) {
  gx <- seq(-lim, lim, by = step)
  list(x = as.matrix(expand.grid(x1 = gx, x2 = gx)), cell = step^2)
}
grid_mass <- function(target, lim) {
  g <- grid_points(lim)
  sum(exp(target@log_density(g$x))) * g$cell
}
kl_grid <- function(target, fit, lim) {
  g <- grid_points(lim)
  log_f <- target@log_density(g$x)
  log_g <- dgmm(g$x, fit, log = TRUE)
  dens_f <- exp(log_f)
  ok <- is.finite(log_f) & is.finite(log_g) & dens_f > 1e-300
  sum(dens_f[ok] * (log_f[ok] - log_g[ok])) * g$cell
}
grid_lim <- c(banana = 12, donut = 7, mixture = 8)
format(round(c(
  banana = grid_mass(banana_target(), grid_lim[["banana"]]),
  donut = grid_mass(donut_target(), grid_lim[["donut"]]),
  mixture = grid_mass(mixture_target(), grid_lim[["mixture"]])
), 6L), nsmall = 6L)
kl_quad <- c(
  banana = kl_grid(banana_target(), fit_b, grid_lim[["banana"]]),
  donut = kl_grid(donut_target(), fit_d, grid_lim[["donut"]]),
  mixture = kl_grid(mixture_target(), fit_m, grid_lim[["mixture"]])
)

## ----summary-table, echo = FALSE----------------------------------------------
fits <- list(fit_b, fit_d, fit_m)
hell <- lapply(fits, hellinger_mc, n_mc = 2000L, seed = 1L)
val_kl <- vapply(fits, function(f) f@diagnostics$validation_kld, numeric(1L))
val_se <- vapply(fits, function(f) f@diagnostics$validation_mc_se,
                 numeric(1L))
hell_h2 <- vapply(hell, function(h) h$h2, numeric(1L))
hell_se <- vapply(hell, function(h) h$se, numeric(1L))
## distance of each held-out estimate from its grid value, and of each
## Hellinger estimate from zero, in standard errors
val_z_max <- max(abs(val_kl - kl_quad) / val_se)
hell_z_min <- min(hell_h2 / hell_se)
hi <- which.max(hell_h2)
lo <- which.min(hell_h2)
hell_gap <- hell_h2[hi] - hell_h2[lo]
hell_gap_se <- sqrt(hell_se[hi]^2 + hell_se[lo]^2)
## held-out difference between two shapes, in standard errors of the difference
pair_z <- function(i, j) abs(val_kl[i] - val_kl[j]) / sqrt(val_se[i]^2 + val_se[j]^2)
z_donut <- min(pair_z(2L, 1L), pair_z(2L, 3L))
z_banana_peaks <- pair_z(1L, 3L)
## the text below states each of these
stopifnot(z_donut > 2, z_banana_peaks < 2, hell_gap < hell_gap_se)
summary_tbl <- data.frame(
  Shape = c("banana", "donut", "three peaks"),
  Components = vapply(fits, gmm_n_components, integer(1L)),
  `Trial points` = vapply(fits, function(f) f@diagnostics$is_size,
                          numeric(1L)),
  `Effective sample size` = round(vapply(fits, function(f) f@diagnostics$ess,
                                         numeric(1L)), 0L),
  `KL, held-out` = round(val_kl, 4L),
  `KL SE` = round(val_se, 4L),
  `KL, grid` = round(unname(kl_quad), 4L),
  `Squared Hellinger` = round(hell_h2, 4L),
  `Hellinger SE` = round(hell_se, 4L),
  check.names = FALSE,
  stringsAsFactors = FALSE
)
knitr::kable(
  summary_tbl,
  caption = paste(
    "Fit quality on the three shapes that cover the plane. The held-out",
    "and grid columns measure the same KL divergence, the first from fresh",
    "trial points and the second by summing over a grid. SE is the standard",
    "error."
  )
)

## ----donut-sweep--------------------------------------------------------------
donut_n <- c(3L, 4L, 6L, 10L)
donut_draws <- c(3500L, 20000L)
donut_kl <- sapply(donut_draws, function(m) {
  vapply(donut_n, function(k) {
    f <- fit_proxymix(donut_target(), N = k, regime = "kld",
                      proposal = proposal_mvt(n_dim = 2L,
                                              sigma = 9 * diag(2),
                                              df = 5),
                      is_size = m, max_iter = 60L, seed = 1L)
    kl_grid(donut_target(), f, grid_lim[["donut"]])
  }, numeric(1L))
})

## ----donut-sweep-table, echo = FALSE------------------------------------------
knitr::kable(
  data.frame(
    Components = donut_n,
    `KL, 3,500 trial points` = round(donut_kl[, 1L], 4L),
    `KL, 20,000 trial points` = round(donut_kl[, 2L], 4L),
    check.names = FALSE
  ),
  caption = paste(
    "KL divergence of the donut proxy, computed on the grid, for four",
    "numbers of components and two numbers of trial points."
  )
)

## ----epan-target--------------------------------------------------------------
epan <- epanechnikov_target(n_dim = 1L)
epan

## ----epan-fit-----------------------------------------------------------------
fit_e <- fit_proxymix(epan, N = 3L, regime = "kld",
                      is_size = 4000L, max_iter = 300L, seed = 1L)
c(proposal = fit_e@diagnostics$proposal_name,
  support_fraction = fit_e@diagnostics$support_fraction,
  iterations = length(kld_trace(fit_e)),
  converged = gmm_fit_quality(fit_e)$converged)

## ----epan-grid----------------------------------------------------------------
epan_x <- seq(-1.4, 1.4, length.out = 400L)
epan_df <- rbind(
  data.frame(x = epan_x, density = exp(epan@log_density(
    matrix(epan_x, ncol = 1L))), series = "target", stringsAsFactors = FALSE),
  data.frame(x = epan_x, density = dgmm(matrix(epan_x, ncol = 1L), fit_e),
             series = "mixture proxy", stringsAsFactors = FALSE)
)
## probability the proxy places outside [-1, 1]
leak_x <- seq(-6, 6, by = 0.001)
leak <- sum(dgmm(matrix(leak_x, ncol = 1L), fit_e)[abs(leak_x) > 1]) * 0.001

## ----epan-plot, eval = has_ggplot2, echo = has_ggplot2, fig.height = 3.4, fig.cap = sprintf("The Epanechnikov target and its three-component proxy. The proxy follows the parabola across the interval with small ripples, and places %s per cent of its probability outside the interval, beyond the grey lines at -1 and 1.", format(signif(100 * leak, 2L))), fig.alt = "Line plot of the upside-down parabola of the Epanechnikov density and its dashed Gaussian-mixture proxy, which extends slightly beyond the interval from -1 to 1."----
ggplot2::ggplot(epan_df, ggplot2::aes(x, density, colour = series,
                                      linetype = series)) +
  ggplot2::geom_line(linewidth = 0.8) +
  ggplot2::geom_vline(xintercept = c(-1, 1), colour = "grey60",
                      linewidth = 0.3) +
  ggplot2::scale_colour_manual(
    name = NULL, values = c("target" = "#0072B2",
                            "mixture proxy" = "#D55E00")
  ) +
  ggplot2::scale_linetype_manual(
    name = NULL, values = c("target" = "solid", "mixture proxy" = "dashed")
  ) +
  ggplot2::labs(x = expression(x), y = "density") +
  ggplot2::theme_minimal(base_size = 11) +
  ggplot2::theme(legend.position = "top")

## ----epan-plot-skip, eval = !has_ggplot2, echo = FALSE, results = "asis"------
# cat("ggplot2 is not installed on this build, so the bounded-support",
#     "figure is skipped.\n")

## ----refusal------------------------------------------------------------------
bad <- tryCatch(
  fit_kld_em(
    mixture_target(), N = 3L,
    proposal = proposal_mvn(n_dim = 2L, mean = c(25, 25), cov = diag(2)),
    is_size = 2000L, max_iter = 10L, seed = 1L,
    min_ess = 50, on_low_ess = "abort"
  ),
  error = function(e) e
)
class(bad)[1L]
cat(conditionMessage(bad))

## ----session-info, collapse = FALSE, class.output = "session-info"------------
sessionInfo()

