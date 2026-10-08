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
res <- readRDS("results/operator_calculus.rds")
major_minor <- function(v) paste(unlist(package_version(v))[1:2],
                                 collapse = ".")
if (major_minor(res$proxymix_version) !=
    major_minor(as.character(packageVersion("proxymix")))) {
  stop("results/operator_calculus.rds was built under proxymix ",
       res$proxymix_version, ", but this is proxymix ",
       packageVersion("proxymix"), ". Rerun the simulation and ",
       "data-raw/vignette_results/operator_calculus.R.", call. = FALSE)
}

## Small numbers are written as plain decimals rather than in the
## scientific notation that knitr's inline hook would otherwise use.
fixed <- function(v, digits) {
  format(round(v, digits), nsmall = digits, scientific = FALSE)
}
sci <- function(v) formatC(v, format = "e", digits = 1L)
gap_text <- function(v) if (v == 0) "exactly" else paste("to", sci(v))
n_word <- function(k) {
  c("one", "two", "three", "four", "five", "six", "seven", "eight", "nine",
    "ten")[k]
}

## ----seed---------------------------------------------------------------------
set.seed(20260514)

## ----prior--------------------------------------------------------------------
g_prior <- gmm(
  weights = c(0.6, 0.4),
  means = list(c(-1, 0), c(1.5, 0.5)),
  covariances = list(diag(c(0.6, 0.8)), diag(c(0.7, 0.5)))
)
g_prior

## ----sensor-------------------------------------------------------------------
A_sensor <- matrix(
  c(1, 0,
    0, 1,
    1, 1),
  nrow = 3L, byrow = TRUE
)
b_sensor <- c(0, 0, 0)
R_sensor <- 0.05 * diag(3)

g_pushed <- gmm_affine(
  g_prior, A_sensor, b_sensor, noise_cov = R_sensor
)

## ----sensor-check-------------------------------------------------------------
mu_hand <- lapply(g_prior@means, function(mu) {
  as.numeric(A_sensor %*% mu + b_sensor)
})
cov_hand <- lapply(g_prior@covariances, function(s) {
  A_sensor %*% s %*% t(A_sensor) + R_sensor
})
affine_gap <- max(
  abs(unlist(g_pushed@means) - unlist(mu_hand)),
  abs(unlist(g_pushed@covariances) - unlist(cov_hand)),
  abs(g_pushed@weights - g_prior@weights)
)

## ----observe------------------------------------------------------------------
A_obs <- matrix(c(1, 0), nrow = 1L)
g_post <- gmm_observe(
  g_prior, A = A_obs, y = 0.8, noise_cov = matrix(0.25, 1L, 1L)
)
g_post

## ----fig-prior-posterior, eval = has_ggplot2, echo = has_ggplot2, fig.height = 4, fig.cap = sprintf("The mixture before (prior) and after (posterior) the first variable is measured as 0.8. The right-hand component, centred at $x_1 = %s$, is closer to the measurement, and its weight rises from %s to %s. Within each component, the variance of $x_1$ shrinks and the mean of $x_1$ moves toward the measurement. The covariance matrices are diagonal, so the mean and variance of $x_2$ in each component do not change.", format(g_prior@means[[2L]][1L]), format(g_prior@weights[2L]), format(round(g_post@weights[2L], 3L))), fig.alt = "Two side-by-side density maps on the same axes. The right-hand posterior panel is more concentrated than the left-hand prior panel, and more of its mass sits in the right-hand component."----
grid <- expand.grid(
  x = seq(-4, 5, length.out = 80L),
  y = seq(-3, 3, length.out = 60L)
)
gm <- as.matrix(grid)
long <- rbind(
  data.frame(x = grid$x, y = grid$y, d = dgmm(gm, g_prior), part = "Prior"),
  data.frame(
    x = grid$x, y = grid$y, d = dgmm(gm, g_post), part = "Posterior"
  )
)
long$part <- factor(long$part, levels = c("Prior", "Posterior"))
ggplot2::ggplot(long, ggplot2::aes(x, y)) +
  ggplot2::geom_raster(ggplot2::aes(fill = d), interpolate = TRUE) +
  ggplot2::geom_contour(
    ggplot2::aes(z = d), colour = "white", linewidth = 0.2,
    alpha = 0.6, bins = 8L
  ) +
  ggplot2::facet_wrap(~ part) +
  ggplot2::coord_equal(expand = FALSE) +
  ggplot2::scale_fill_viridis_c(name = "density") +
  ggplot2::labs(
    x = expression(x[1]), y = expression(x[2]),
    title = "Before and after measuring the first variable"
  ) +
  ggplot2::theme_minimal(base_size = 11) +
  ggplot2::theme(
    strip.text = ggplot2::element_text(face = "bold"),
    panel.grid = ggplot2::element_blank()
  )

## ----fig-prior-posterior-skip, eval = !has_ggplot2, echo = FALSE, results = "asis"----
# cat("ggplot2 is not installed on this build, so the prior-and-posterior",
#     "figure is skipped.\n")

## ----kalman-parity------------------------------------------------------------
g_single <- gmm(
  weights = 1, means = list(c(0, 0)),
  covariances = list(diag(c(1, 2)))
)
g_one_obs <- gmm_observe(
  g_single, A = matrix(c(1, 0), nrow = 1L), y = 0.5,
  noise_cov = matrix(0.5, 1L, 1L)
)

s_prior <- diag(c(1, 2))
h_obs <- matrix(c(1, 0), nrow = 1L)
r_obs <- matrix(0.5, 1L, 1L)
s_innov <- h_obs %*% s_prior %*% t(h_obs) + r_obs
gain <- s_prior %*% t(h_obs) %*% solve(s_innov)
mu_kalman <- as.numeric(gain * 0.5)
cov_kalman <- s_prior - gain %*% h_obs %*% s_prior

kalman_gap <- max(
  abs(mu_kalman - g_one_obs@means[[1L]]),
  abs(cov_kalman - g_one_obs@covariances[[1L]])
)
ridge_default <- 1e-6

## ----aggregate----------------------------------------------------------------
g_fine <- gmm(
  weights = c(0.3, 0.4, 0.3),
  means = list(c(0, 0, 0), c(2, 1, -1), c(-1, -1, 2)),
  covariances = list(diag(3), diag(3), diag(3))
)
A_agg <- matrix(
  c(1, 1, 0,
    0, 0, 1),
  nrow = 2L, byrow = TRUE
)
g_coarse <- gmm_aggregate(g_fine, A_agg)
weights_kept <- max(abs(g_coarse@weights - g_fine@weights))

## ----aggregate-kable, echo = FALSE--------------------------------------------
knitr::kable(
  data.frame(
    component = seq_len(gmm_n_components(g_coarse)),
    weight = round(g_coarse@weights, 3L),
    mean_1 = round(vapply(g_coarse@means, function(m) m[1L], numeric(1L)), 3L),
    mean_2 = round(vapply(g_coarse@means, function(m) m[2L], numeric(1L)), 3L)
  ),
  col.names = c("Component", "Weight", "Mean of x1 + x2", "Mean of x3"),
  caption = paste(
    "The mixture after adding the first two variables together. It has the",
    "same number of components and the same weights as before. The means",
    "and covariance matrices are passed through the summing matrix."
  )
)

## ----conditioning-------------------------------------------------------------
g_cond_index <- gmm_missing(g_prior, observed = 2L, values = 0.5)
g_cond_given <- gmm_conditionalise(g_prior, given = c(NA, 0.5))
cond_gap <- max(
  abs(unlist(g_cond_index@means) - unlist(g_cond_given@means)),
  abs(unlist(g_cond_index@covariances) -
        unlist(g_cond_given@covariances)),
  abs(g_cond_index@weights - g_cond_given@weights)
)

## ----compose------------------------------------------------------------------
g_a <- gmm_observe(
  g_prior, A = matrix(c(1, 0), nrow = 1L), y = 0.5,
  noise_cov = matrix(0.25, 1L, 1L)
)
g_ab <- gmm_observe(
  g_a, A = matrix(c(0, 1), nrow = 1L), y = 0.2,
  noise_cov = matrix(0.25, 1L, 1L)
)
g_stack <- gmm_observe(
  g_prior, A = diag(2), y = c(0.5, 0.2), noise_cov = 0.25 * diag(2)
)
compose_gap <- max(
  abs(g_ab@weights - g_stack@weights),
  abs(unlist(g_ab@means) - unlist(g_stack@means)),
  abs(unlist(g_ab@covariances) - unlist(g_stack@covariances))
)

## ----track--------------------------------------------------------------------
dt <- 1
A_dyn <- matrix(c(1, dt, 0, 1), 2L, 2L, byrow = TRUE)
C_obs <- matrix(c(1, 0), 1L, 2L)
Q_proc <- 0.01 * diag(2)
R_meas <- matrix(0.5, 1L, 1L)

n_steps <- 30L
truth <- matrix(0, n_steps, 2L)
truth[1L, ] <- c(0, 1)
for (k1 in 2:n_steps) {
  truth[k1, ] <- as.numeric(A_dyn %*% truth[k1 - 1L, ]) +
    mvnfast::rmvn(1L, c(0, 0), Q_proc)
} # ends k1, over the simulated state track
y_track <- truth[, 1L] + rnorm(n_steps, 0, sqrt(R_meas[1L, 1L]))

## ----filter-loop--------------------------------------------------------------
g_state <- gmm(
  weights = 1, means = list(c(0, 0)), covariances = list(diag(2))
)
filtered <- numeric(n_steps)
for (k1 in seq_len(n_steps)) {
  if (k1 > 1L) {
    g_state <- gmm_affine(
      g_state, A = A_dyn, b = c(0, 0), noise_cov = Q_proc
    )                                                       # predict
  }
  g_state <- gmm_observe(
    g_state, A = C_obs, y = y_track[k1], noise_cov = R_meas
  )                                                         # update
  filtered[k1] <- g_state@means[[1L]][1L]
} # ends k1, over the filtering recursion

## ----filter-parity------------------------------------------------------------
mu_kf <- c(0, 0)
p_kf <- diag(2)
kf_track <- numeric(n_steps)
for (k1 in seq_len(n_steps)) {
  if (k1 > 1L) {
    mu_kf <- as.numeric(A_dyn %*% mu_kf)
    p_kf <- A_dyn %*% p_kf %*% t(A_dyn) + Q_proc            # predict
  }
  gain_kf <- p_kf %*% t(C_obs) %*%
    solve(C_obs %*% p_kf %*% t(C_obs) + R_meas)             # gain
  mu_kf <- mu_kf + as.numeric(gain_kf %*% (y_track[k1] - C_obs %*% mu_kf))
  p_kf <- (diag(2) - gain_kf %*% C_obs) %*% p_kf
  kf_track[k1] <- mu_kf[1L]
} # ends k1, over the hand-coded Kalman recursion
loop_gap <- max(abs(filtered - kf_track))

## ----fig-kalman, eval = has_ggplot2, echo = has_ggplot2, fig.height = 3.6, fig.cap = "An object moving along a line, tracked by the predict and update loop with one component, which is the Kalman filter. The points are the noisy position readings, and the two lines are the true position and the filtered estimate.", fig.alt = "A time series with scattered grey noisy position readings, a line for the true position, and a filtered-estimate line closely following it."----
track_df <- data.frame(
  t = seq_len(n_steps), truth = truth[, 1L], y = y_track,
  filtered = filtered
)
ggplot2::ggplot(track_df, ggplot2::aes(t)) +
  ggplot2::geom_point(
    ggplot2::aes(y = y, colour = "noisy reading"), size = 1.3, alpha = 0.7
  ) +
  ggplot2::geom_line(
    ggplot2::aes(y = truth, colour = "true position"), linewidth = 0.8
  ) +
  ggplot2::geom_line(
    ggplot2::aes(y = filtered, colour = "filtered, one component"),
    linewidth = 0.9
  ) +
  ggplot2::scale_colour_manual(
    name = NULL,
    values = c(
      "noisy reading" = "grey60",
      "true position" = "#0072B2",
      "filtered, one component" = "#D55E00"
    )
  ) +
  ggplot2::labs(
    x = "time step", y = "position",
    title = "Predict and update over time: the Kalman filter"
  ) +
  ggplot2::theme_minimal(base_size = 11) +
  ggplot2::theme(legend.position = "top")

## ----fig-kalman-skip, eval = !has_ggplot2, echo = FALSE, results = "asis"-----
# cat("ggplot2 is not installed on this build, so the Kalman-track figure",
#     "is skipped.\n")

## ----reduce-------------------------------------------------------------------
g_six <- gmm(
  weights = rep(1 / 6, 6L),
  means = list(
    c(-5, 0), c(-5, 0.15), c(5, 0), c(5.1, -0.1), c(0, 6), c(0.1, 6.1)
  ),
  covariances = rep(list(0.5 * diag(2)), 6L)
)
g_three <- gmm_reduce(g_six, k_max = 3L)

mix_mean <- function(g) Reduce(`+`, Map(`*`, g@weights, g@means))
reduce_shift <- max(abs(mix_mean(g_six) - mix_mean(g_three)))
reduce_divergence <- gmm_divergence(g_six, g_three, type = "cs")

## ----reduce-kable, echo = FALSE-----------------------------------------------
knitr::kable(
  data.frame(
    quantity = c(
      "components before", "components after",
      "change in the mean of the mixture",
      "Cauchy-Schwarz divergence from the original"
    ),
    value = c(
      formatC(gmm_n_components(g_six), format = "d"),
      formatC(gmm_n_components(g_three), format = "d"),
      formatC(reduce_shift, format = "e", digits = 1L),
      formatC(reduce_divergence, format = "e", digits = 1L)
    )
  ),
  col.names = c("Quantity", "Value"),
  caption = paste(
    "Reducing six components to three by merging pairs, and how much the",
    "mixture changes."
  )
)

## ----verb-kalman--------------------------------------------------------------
prior_state <- gmm(
  weights = 1, means = list(c(0, 0)), covariances = list(diag(2))
)
out_verb <- gmm_filter(
  prior_state,
  dynamics = list(A = A_dyn, Q = Q_proc),
  measurement = list(C = C_obs, R = R_meas),
  y = y_track, ridge_eps = 0
)

mu_v <- c(0, 0)
p_v <- diag(2)
kf_verb <- numeric(n_steps)
for (k1 in seq_len(n_steps)) {
  mu_v <- as.numeric(A_dyn %*% mu_v)
  p_v <- A_dyn %*% p_v %*% t(A_dyn) + Q_proc                # predict
  gain_v <- p_v %*% t(C_obs) %*%
    solve(C_obs %*% p_v %*% t(C_obs) + R_meas)              # gain
  mu_v <- mu_v + as.numeric(gain_v %*% (y_track[k1] - C_obs %*% mu_v))
  p_v <- (diag(2) - gain_v %*% C_obs) %*% p_v               # update
  kf_verb[k1] <- mu_v[1L]
} # ends k1, over the predict-then-update reference recursion
verb_gap <- max(abs(out_verb$mean[, 1L] - kf_verb))

## ----verb-gsf-----------------------------------------------------------------
q_heavy <- gmm(
  weights = c(0.9, 0.1),
  means = list(c(0, 0), c(0, 0)),
  covariances = list(0.01 * diag(2), 0.5 * diag(2))
)
out_gsf <- gmm_filter(
  prior_state,
  dynamics = list(A = A_dyn, Q = q_heavy),
  measurement = list(C = C_obs, R = R_meas),
  y = y_track, k_max = 6L
)
gsf_max_k <- max(out_gsf$summary$n_components)
gsf_rmse <- sqrt(mean((out_gsf$mean[, 1L] - truth[, 1L])^2))
kalman_rmse <- sqrt(mean((out_verb$mean[, 1L] - truth[, 1L])^2))

## ----gsf-counts, include = FALSE----------------------------------------------
gsf_k <- out_gsf$summary$n_components
gsf_first_cap <- which(gsf_k == gsf_max_k)[1L]
gsf_early <- paste(paste0(gsf_k[seq_len(gsf_first_cap - 1L)],
                          " after step ", seq_len(gsf_first_cap - 1L)),
                   collapse = ", ")

## ----fig-gsf, eval = has_ggplot2, echo = has_ggplot2, fig.height = 3.6, fig.cap = sprintf("The same track filtered by the Gaussian-sum filter with a two-component process noise, capped at six components. The number of components is %s, and %s from step %s to step %s.", gsf_early, gsf_max_k, gsf_first_cap, n_steps), fig.alt = "A time series with grey noisy position readings, a line for the true position, and a Gaussian-sum-filter estimate line following it closely."----
gsf_df <- data.frame(
  t = seq_len(n_steps), truth = truth[, 1L], y = y_track,
  filtered = out_gsf$mean[, 1L]
)
ggplot2::ggplot(gsf_df, ggplot2::aes(t)) +
  ggplot2::geom_point(
    ggplot2::aes(y = y, colour = "noisy reading"), size = 1.3, alpha = 0.7
  ) +
  ggplot2::geom_line(
    ggplot2::aes(y = truth, colour = "true position"), linewidth = 0.8
  ) +
  ggplot2::geom_line(
    ggplot2::aes(y = filtered, colour = "Gaussian-sum filter"),
    linewidth = 0.9
  ) +
  ggplot2::scale_colour_manual(
    name = NULL,
    values = c(
      "noisy reading" = "grey60",
      "true position" = "#0072B2",
      "Gaussian-sum filter" = "#D55E00"
    )
  ) +
  ggplot2::labs(
    x = "time step", y = "position",
    title = "A Gaussian-sum filter capped at six components"
  ) +
  ggplot2::theme_minimal(base_size = 11) +
  ggplot2::theme(legend.position = "top")

## ----fig-gsf-skip, eval = !has_ggplot2, echo = FALSE, results = "asis"--------
# cat("ggplot2 is not installed on this build, so the Gaussian-sum-filter",
#     "figure is skipped.\n")

## ----parity-kable, echo = FALSE-----------------------------------------------
knitr::kable(
  data.frame(
    check = c(
      "linear map against the hand formula",
      "one-component update against a hand Kalman update",
      "conditioning by position against conditioning by value",
      "two measurements in turn against both at once",
      "predict and update loop against a hand Kalman filter",
      "`gmm_filter()`, `ridge_eps = 0`, against a hand Kalman filter",
      "weights after adding variables together"
    ),
    gap = formatC(
      c(
        affine_gap, kalman_gap, cond_gap, compose_gap, loop_gap,
        verb_gap, weights_kept
      ),
      format = "e", digits = 1L
    )
  ),
  col.names = c("Check", "Largest absolute difference"),
  caption = paste(
    "Each exact operation against a calculation written out by hand.",
    "`gmm_affine()` and `gmm_observe()` add a ridge of $10^{-6}$ to each",
    "covariance matrix unless `ridge_eps = 0`."
  )
)

## ----compare-facts, include = FALSE-------------------------------------------
sim_value <- function(method, what) {
  res$sim_tab[[what]][res$sim_tab$method == method]
}
pm <- "proxymix, two-component noise"
pf <- "particle filter, two-component noise"
secs_pm <- sim_value(pm, "secs")
ratio_to <- function(method) secs_pm / sim_value(method, "secs")
z_t <- res$paired$student_t[["diff"]] / res$paired$student_t[["se"]]

## ----compare-table, echo = FALSE----------------------------------------------
method_order <- c(pm, pf, "nimble, two-component noise",
                  "nimble, Student-t noise", "dlm, Gaussian noise",
                  "KFAS, Gaussian noise")
cmp_tbl <- res$sim_tab[match(method_order, res$sim_tab$method), ]
knitr::kable(
  data.frame(
    filter = c("proxymix", "particle filter, by hand",
               "particle filter, nimbleSMC", "Student-t filter, nimbleSMC",
               "dlm", "KFAS"),
    noise = c("two-component", "two-component", "two-component",
              "Student-t", "normal", "normal"),
    rmse = fixed(cmp_tbl$rmse, 4),
    rmse_se = fixed(cmp_tbl$rmse_se, 4),
    log_lik = fixed(cmp_tbl$log_lik, 1),
    secs = fixed(cmp_tbl$secs, 4),
    stringsAsFactors = FALSE
  ),
  align = c("l", "l", "r", "r", "r", "r"), row.names = FALSE,
  col.names = c("Filter", "Noise model", "RMSE", "Standard error",
                "Log-likelihood", "Seconds"),
  caption = paste0(
    "Filtering ", res$n_rep, " simulated series of ", res$n_t, " steps ",
    "with outliers in the readings. RMSE is the root mean squared error ",
    "of the filtered level, averaged over the series, with its standard ",
    "error. Log-likelihood and seconds are averages per series."
  )
)

## ----compare-code, eval = FALSE-----------------------------------------------
# library(proxymix)
# library(dlm)
# library(KFAS)
# 
# n_t <- 200L           # steps per series
# n_particles <- 10000L
# q_state <- 0.25       # variance of the level increments
# r_sd <- c(1, 5)       # noise standard deviations, core and outlier
# r_w <- c(0.9, 0.1)    # their weights
# r_var <- sum(r_w * r_sd^2)
# m0 <- 0
# c0 <- 10
# 
# # one series: a random-walk level and readings with occasional outliers
# set.seed(1L)
# x <- m0 + sqrt(c0) * rnorm(1L) + cumsum(rnorm(n_t, 0, sqrt(q_state)))
# outlier <- runif(n_t) < r_w[2L]
# y <- x + rnorm(n_t, 0, ifelse(outlier, r_sd[2L], r_sd[1L]))
# 
# # proxymix: the two-component noise, capped at four components per step
# prior_sim <- gmm(weights = 1, means = list(m0), covariances = list(matrix(c0)))
# r_mix_sim <- gmm(weights = r_w, means = list(0, 0),
#                  covariances = list(matrix(r_sd[1L]^2), matrix(r_sd[2L]^2)))
# f_pm <- gmm_filter(prior_sim,
#                    dynamics = list(A = matrix(1), Q = matrix(q_state)),
#                    measurement = list(C = matrix(1), R = r_mix_sim),
#                    y = y, k_max = 4L)
# 
# # dlm and KFAS: normal noise with the same variance
# mod_dlm_sim <- dlmModPoly(1L, dV = r_var, dW = q_state, m0 = m0, C0 = c0)
# f_dlm <- dlmFilter(y, mod_dlm_sim)
# mod <- SSModel(y ~ SSMtrend(1L, Q = q_state, a1 = m0,
#                             P1 = c0 + q_state), H = r_var)
# f_kfas <- KFS(mod, filtering = "state", smoothing = "none")
# 
# # a bootstrap particle filter under the two-component noise
# bootstrap_filter <- function(y, n_particles, m0, c0, q, r_sd, r_w) {
#   n <- length(y)
#   particles <- rnorm(n_particles, m0, sqrt(c0))
#   filtered <- numeric(n)
#   log_lik <- 0
#   for (t in seq_len(n)) {
#     particles <- particles + rnorm(n_particles, 0, sqrt(q))
#     lik <- r_w[1L] * dnorm(y[t], particles, r_sd[1L]) +
#       r_w[2L] * dnorm(y[t], particles, r_sd[2L])
#     log_lik <- log_lik + log(mean(lik))
#     filtered[t] <- sum(lik * particles) / sum(lik)
#     particles <- particles[sample.int(n_particles, n_particles,
#                                       replace = TRUE, prob = lik)]
#   }
#   list(mean = filtered, log_lik = log_lik)
# }
# f_pf <- bootstrap_filter(y, n_particles, m0, c0, q_state, r_sd, r_w)
# 
# # root mean squared error of each filtered level against the true level
# rmse <- function(m) sqrt(mean((m - x)^2))
# c(proxymix = rmse(f_pm$mean[, 1L]),
#   dlm = rmse(as.numeric(f_dlm$m[-1L])),
#   KFAS = rmse(as.numeric(f_kfas$att)),
#   particle = rmse(f_pf$mean))

## ----session-info, collapse = FALSE, class.output = "session-info"------------
sessionInfo()

