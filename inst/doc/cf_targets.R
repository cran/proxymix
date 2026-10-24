## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width = 7,
  fig.height = 4,
  dpi = 150,
  out.width = "100%"
)
## One BLAS thread, so that the check of this vignette does not take more CPU
## time than elapsed time under a multithreaded BLAS.
if (requireNamespace("RhpcBLASctl", quietly = TRUE)) {
  RhpcBLASctl::blas_set_num_threads(1L)
}

## ----library------------------------------------------------------------------
library(proxymix)

## ----helpers, include = FALSE-------------------------------------------------
## Small numbers are written as plain decimals rather than in the
## scientific notation that knitr's inline hook would otherwise use.
fixed <- function(v, digits) {
  format(round(v, digits), nsmall = digits, scientific = FALSE)
}

## ----sum-target---------------------------------------------------------------
cf_sum <- function(t) {
  s <- t[, 1]
  uniform <- ifelse(s == 0, 1, sin(s) / s)
  (1 - 1i * s)^(-3) * uniform
}
tgt <- cf_target(cf_sum, n_dim = 1L, name = "Gamma(3) plus uniform")
fit_sum <- fit_proxymix(tgt, N = 4L, n_freq = 2048L, seed = 1L)
fit_sum

## ----sum-check----------------------------------------------------------------
cdf_sum <- function(x) {
  vapply(x, function(z) {
    stats::integrate(function(u) stats::pgamma(z - u, 3) / 2, -1, 1)$value
  }, numeric(1))
}
x <- seq(-1, 15, length.out = 801L)
exact_gap <- max(abs(cdf_sum(x) - pgmm(x, fit_sum)))
q_sum <- gmm_fit_quality(fit_sum)
c(lower = q_sum$cf_distance, `fine grid` = exact_gap, upper = q_sum$cf_upper)

## ----sum-wide-----------------------------------------------------------------
wide <- gmm_cf_distance(fit_sum, cf_sum, probs = c(1e-4, 1 - 1e-4),
                        n_grid = 801L)
c(lower = wide$distance, upper = wide$upper)

## ----sum-plot, fig.alt = "Difference between the exact distribution function of the Gamma plus uniform sum and that of the four-component mixture, plotted against x. The difference reaches the dashed line at minus the lower bound near x = 0 and stays well inside the dotted lines at plus and minus the upper bound."----
plot(x, cdf_sum(x) - pgmm(x, fit_sum), type = "l", lwd = 1.5,
     ylim = c(-1, 1) * q_sum$cf_upper, xlab = "x",
     ylab = "exact minus mixture distribution function")
abline(h = c(-1, 1) * q_sum$cf_distance, lty = 2)
abline(h = c(-1, 1) * q_sum$cf_upper, lty = 3)
abline(h = 0, col = "grey60")
legend("topright", c("lower bound", "upper bound"), lty = c(2, 3),
       bty = "n")

## ----stable-fit---------------------------------------------------------------
cf_stable <- function(t) exp(-abs(t[, 1])^1.5)
fit_stable <- from_cf(cf_stable, n_dim = 1L, N = 6L, n_freq = 2048L,
                      seed = 1L)

## ----stable-elliptical--------------------------------------------------------
fit_ell <- from_cf_elliptical(cf_stable, n_dim = 1L, seed = 1L)
fit_ell

## ----stable-exact-------------------------------------------------------------
## symmetric stable law with characteristic function exp(-|t|^alpha),
## for 1 < alpha <= 2
pstable_sym <- function(x, alpha) {
  r <- alpha / (alpha - 1)
  v <- function(th) {
    (cos(th) / sin(alpha * th))^r * cos((alpha - 1) * th) / cos(th)
  }
  vapply(x, function(z) {
    if (z == 0) return(0.5)
    upper_tail <- stats::integrate(function(th) exp(-abs(z)^r * v(th)),
                                   0, pi / 2, rel.tol = 1e-10)$value / pi
    if (z > 0) 1 - upper_tail else upper_tail
  }, numeric(1))
}
probs <- c(0.5, 0.75, 0.95, 0.99, 0.999)
exact_q <- vapply(probs, function(p) {
  if (p == 0.5) return(0)
  stats::uniroot(function(x) pstable_sym(x, 1.5) - p, c(0, 200),
                 tol = 1e-10)$root
}, numeric(1))
quantiles <- rbind(
  exact = exact_q,
  `from_cf, 6 components` = qgmm(probs, fit_stable),
  from_cf_elliptical = qgmm(probs, fit_ell)
)
colnames(quantiles) <- probs
round(quantiles, 3)

## ----stable-bounds------------------------------------------------------------
x <- seq(-60, 60, length.out = 4001L)
exact_cdf <- pstable_sym(x, 1.5)
bounds <- sapply(list(`from_cf` = fit_stable, from_cf_elliptical = fit_ell),
                 function(fit) {
  q <- gmm_fit_quality(fit)
  c(lower = q$cf_distance, `fine grid` = max(abs(exact_cdf - pgmm(x, fit))),
    upper = q$cf_upper)
})
signif(bounds, 3)

## ----stable-plot, fig.height = 3.6, fig.alt = "Left: the density of the scale-mixture fit drawn over a histogram of exact draws from the stable law of index 1.5; the curve follows the histogram. Right: upper tail probability of the stable law and of the two mixture fits on log scales; the scale mixture stays close to the stable tail up to x = 100, while the six-component fit drops below it beyond about x = 25."----
set.seed(2)
u <- stats::runif(1e5, -pi / 2, pi / 2)
e <- stats::rexp(1e5)
draws <- sin(1.5 * u) / cos(u)^(1 / 1.5) *
  (cos(-0.5 * u) / e)^(-0.5 / 1.5)
op <- par(mfrow = c(1, 2), mar = c(4, 4, 1.5, 1))
## density over all draws, shown on (-8, 8)
h <- hist(draws[abs(draws) < 8], breaks = seq(-8, 8, by = 0.2), plot = FALSE)
h$density <- h$counts / (length(draws) * 0.2)
plot(h, freq = FALSE, col = "grey80", border = NA, main = "", xlab = "x",
     ylab = "density", ylim = c(0, 0.38))
xs <- seq(-8, 8, length.out = 401L)
lines(xs, dgmm(matrix(xs), fit_ell), lwd = 2)
legend("topleft", c("exact draws", "from_cf_elliptical"),
       fill = c("grey80", NA), border = NA, lwd = c(NA, 2), bty = "n",
       cex = 0.8)
xt <- exp(seq(log(0.5), log(100), length.out = 200L))
## P(X > x) = P(X < -x) by symmetry, which avoids 1 - F(x) rounding to 0
plot(xt, pstable_sym(-xt, 1.5), log = "xy", type = "l", lwd = 3,
     col = "grey60", ylim = c(1e-8, 0.5), xlab = "x", ylab = "P(X > x)")
lines(xt, pgmm(-xt, fit_ell), lwd = 1.5)
lines(xt, pgmm(-xt, fit_stable), lwd = 1.5, lty = 2)
legend("bottomleft", c("stable law", "from_cf_elliptical", "from_cf"),
       lwd = c(3, 1.5, 1.5), lty = c(1, 1, 2),
       col = c("grey60", "black", "black"), bty = "n")
par(op)

## ----cauchy2------------------------------------------------------------------
S <- matrix(c(2, 0.5, 0.5, 1), 2L)
cf_cauchy2 <- function(t) exp(1i * t[, 1] - sqrt(rowSums((t %*% S) * t)))
fit_c2 <- from_cf_elliptical(cf_cauchy2, n_dim = 2L, seed = 1L)
fit_c2
probs <- c(0.01, 0.05, 0.5, 0.95, 0.99)
round(rbind(
  `first, exact` = 1 + sqrt(S[1, 1]) * stats::qcauchy(probs),
  `first, mixture` = qgmm(probs, gmm_marginalise(fit_c2, 1L)),
  `second, exact` = sqrt(S[2, 2]) * stats::qcauchy(probs),
  `second, mixture` = qgmm(probs, gmm_marginalise(fit_c2, 2L))
), 3)

## ----sources------------------------------------------------------------------
A <- matrix(c(0.7, 0.3, 0, -0.7), 2L)
cf_mix <- function(t) {
  u <- t %*% A
  (1 - 1i * u[, 1])^(-3) * (1 - 1i * u[, 2])^(-3) *
    exp(-0.125 * rowSums(t^2))
}
fit_src <- from_cf_sources(cf_mix, n_dim = 2L, N = 3L, n_freq = 1024L,
                           seed = 1L)
fit_src

## ----sources-project----------------------------------------------------------
a <- c(1, 1) / sqrt(2)
along <- fit_src$project(a)
bracket <- gmm_cf_distance(along, function(s) cf_mix(s %*% t(a)))
c(lower = bracket$distance, upper = bracket$upper)

## ----sources-check------------------------------------------------------------
set.seed(3)
n <- 5e5
s_draws <- cbind(stats::rgamma(n, 3), stats::rgamma(n, 3))
x_draws <- s_draws %*% t(A) + matrix(stats::rnorm(2 * n, sd = 0.5), n)
y <- sort(as.vector(x_draws %*% a))
p_mix <- pgmm(y, along)
max(abs(p_mix - seq_len(n) / n), abs(p_mix - (seq_len(n) - 1) / n))

## ----filter-predict-----------------------------------------------------------
prior <- gmm(weights = c(0.5, 0.5), means = list(-1, 2),
             covariances = list(matrix(0.3), matrix(0.5)))
cf_pred <- function(t) gmm_cf(prior, 0.9 * t[, 1]) * exp(-0.2 * abs(t[, 1]))
pred <- suppressWarnings(
  from_cf(cf_pred, n_dim = 1L, N = 6L, n_freq = 2048L, seed = 1L),
  classes = "proxymix_cf_heavy_tail"
)

## ----filter-update------------------------------------------------------------
post <- gmm_observe(pred, A = matrix(1), y = 1.2, noise_cov = matrix(0.25))
post_mean <- sum(post@weights * unlist(post@means))
post_sd <- sqrt(sum(post@weights * (unlist(post@covariances) +
                                      unlist(post@means)^2)) - post_mean^2)

## ----filter-check-------------------------------------------------------------
dens_pred <- function(x) {
  vapply(x, function(z) {
    stats::integrate(function(e) {
      dgmm(matrix((z - e) / 0.9), prior) / 0.9 * stats::dcauchy(e, 0, 0.2)
    }, -Inf, Inf, subdivisions = 500L)$value
  }, numeric(1))
}
grid <- seq(-6, 8, length.out = 1401L)
h <- grid[2] - grid[1]
w <- dens_pred(grid) * stats::dnorm(1.2, grid, 0.5)
w <- w / sum(w * h)
exact_mean <- sum(grid * w * h)
exact_sd <- sqrt(sum((grid - exact_mean)^2 * w * h))
round(rbind(exact = c(mean = exact_mean, sd = exact_sd),
            mixture = c(post_mean, post_sd)), 4)

## ----bracket-table------------------------------------------------------------
fit_sum@diagnostics$bracket

