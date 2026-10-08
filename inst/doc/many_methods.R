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

## ----shared-helpers-----------------------------------------------------------
## Mean of the first variable (y) when the second (x) is fixed at each value
## of xv: the component means of y, weighted by how likely each component is
## at that x.
cond_mean <- function(fit, xv) {
  vapply(xv, function(xx) {
    g <- gmm_conditionalise(fit, given = c(NA, xx))
    sum(g@weights * vapply(g@means, function(m) m[1L], numeric(1L)))
  }, numeric(1L))
}

## ----sci-notation, include = FALSE--------------------------------------------
## Very small differences are typeset as powers of ten in LaTeX.
sci <- function(v, digits = 2L) {
  v <- signif(v, digits)
  e <- floor(log10(abs(v)))
  paste0("$", signif(v / 10^e, digits), " \\times 10^{", e, "}$")
}
sci_plain <- function(v) formatC(v, format = "e", digits = 1L)

## ----reg-fit------------------------------------------------------------------
set.seed(20260617)
n <- 400L
x <- runif(n, -3, 3)
y <- 0.3 * x + 1.2 * pmax(x, 0) + rnorm(n, sd = 0.4)  # bends at x = 0
dat <- data.frame(y = y, x = x)

joint <- gmm_target_from_samples(cbind(y, x))
fit1 <- fit_proxymix(joint, N = 1L, regime = "moment", ridge_eps = 0)
fit3 <- fit_proxymix(joint, N = 3L, regime = "sample", max_iter = 150L)

## With one component, the slope of E[y | x] should equal the lm slope.
slope_mix <- gmm_conditionalise(fit1, given = c(NA, 1))@means[[1L]] -
  gmm_conditionalise(fit1, given = c(NA, 0))@means[[1L]]
slope_lm <- unname(coef(lm(y ~ x, dat))["x"])
diff_reg <- abs(slope_mix - slope_lm)

## ----fig-reg, eval = has_ggplot2, echo = has_ggplot2, fig.cap = "The least-squares line (one component) and the conditional mean of a three-component mixture, on data whose true relationship bends at zero. The mixture follows the bend, and the straight line does not.", fig.alt = "Scatter of y against x with a straight least-squares line and a curved mixture conditional mean that follows a bend in the data at x equal to zero."----
grid_reg <- data.frame(x = seq(-3, 3, length.out = 200L))
grid_reg$ols <- as.numeric(predict(lm(y ~ x, dat), newdata = grid_reg))
grid_reg$mix <- cond_mean(fit3, grid_reg$x)
ggplot2::ggplot() +
  ggplot2::geom_point(data = dat, ggplot2::aes(x, y), colour = "grey60",
                      alpha = 0.4, size = 0.7) +
  ggplot2::geom_line(data = grid_reg,
                     ggplot2::aes(x, ols, colour = "lm (K = 1)"),
                     linewidth = 0.9) +
  ggplot2::geom_line(data = grid_reg,
                     ggplot2::aes(x, mix, colour = "mixture (K = 3)"),
                     linewidth = 0.9) +
  ggplot2::scale_colour_manual(
    name = NULL,
    values = c("lm (K = 1)" = "#0072B2", "mixture (K = 3)" = "#D55E00")
  ) +
  ggplot2::labs(
    x = "x", y = "y",
    title = "Regression: a straight line and a conditioned mixture"
  ) +
  ggplot2::theme_minimal(base_size = 11) +
  ggplot2::theme(legend.position = "top")

## ----fig-reg-skip, eval = !has_ggplot2, echo = FALSE, results = "asis"--------
# cat("ggplot2 is not installed, so this figure is skipped.\n")

## ----nw-fit-------------------------------------------------------------------
h <- 0.4                                     # bandwidth
nw <- function(xq) {
  vapply(xq, function(q) {
    w <- dnorm(q, x, h)                      # Nadaraya-Watson weights
    sum(w * y) / sum(w)
  }, numeric(1L))
}

## One normal component per data point, then condition on x.
kde <- gmm(weights = rep(1 / n, n),
           means = lapply(seq_len(n), function(i) c(y[i], x[i])),
           covariances = rep(list(diag(c(h^2, h^2))), n))
xq <- seq(-2.5, 2.5, length.out = 21L)
diff_nw <- max(abs(nw(xq) - cond_mean(kde, xq)))

## ----fig-nw, eval = has_ggplot2, echo = has_ggplot2, fig.cap = "The two ends of one scale: the least-squares line (one component) and the Nadaraya-Watson smoother (one component per data point). The same conditioning step produces both.", fig.alt = "Scatter of y against x with the straight least-squares line and the curved Nadaraya-Watson kernel-regression line."----
grid_nw <- data.frame(x = seq(-3, 3, length.out = 200L))
grid_nw$ols <- as.numeric(predict(lm(y ~ x, dat), newdata = grid_nw))
grid_nw$nw <- nw(grid_nw$x)
ggplot2::ggplot() +
  ggplot2::geom_point(data = dat, ggplot2::aes(x, y), colour = "grey60",
                      alpha = 0.4, size = 0.7) +
  ggplot2::geom_line(data = grid_nw,
                     ggplot2::aes(x, ols, colour = "least squares (K = 1)"),
                     linewidth = 0.9) +
  ggplot2::geom_line(data = grid_nw,
                     ggplot2::aes(x, nw, colour = "kernel (K = n)"),
                     linewidth = 0.9) +
  ggplot2::scale_colour_manual(
    name = NULL,
    values = c("least squares (K = 1)" = "#0072B2",
               "kernel (K = n)" = "#D55E00")
  ) +
  ggplot2::labs(
    x = "x", y = "y",
    title = "From a straight line to a kernel smoother"
  ) +
  ggplot2::theme_minimal(base_size = 11) +
  ggplot2::theme(legend.position = "top")

## ----fig-nw-skip, eval = !has_ggplot2, echo = FALSE, results = "asis"---------
# cat("ggplot2 is not installed, so this figure is skipped.\n")

## ----clust-fit----------------------------------------------------------------
set.seed(20260617)
x_clust <- rbind(
  mvnfast::rmvn(150L, c(-2, -1), 0.5 * diag(2)),
  mvnfast::rmvn(150L, c(2, 0), matrix(c(0.6, 0.3, 0.3, 0.4), 2L)),
  mvnfast::rmvn(150L, c(0, 2.5), 0.3 * diag(2))
)
colnames(x_clust) <- c("V1", "V2")
target_clust <- gmm_target_from_samples(x_clust)
fitc <- fit_proxymix(target_clust, N = 3L, regime = "sample",
                     max_iter = 150L)

## Responsibility of each component for each row.
responsibilities <- function(fit, xx) {
  comp <- vapply(seq_len(gmm_n_components(fit)), function(k) {
    fit@weights[k] *
      mvnfast::dmvn(xx, mu = fit@means[[k]], sigma = fit@covariances[[k]])
  }, numeric(nrow(xx)))
  comp / rowSums(comp)
}
resp <- responsibilities(fitc, x_clust)
mean_confidence <- mean(apply(resp, 1L, max))

## ----clust-table, echo = FALSE------------------------------------------------
knitr::kable(
  head(round(resp, 3L), 4L),
  col.names = paste("component", seq_len(3L)),
  caption = paste0(
    "Responsibilities of the three components for the first four rows. ",
    "Each row sums to one."
  )
)

## ----pca----------------------------------------------------------------------
fit_pca <- fit_proxymix(target_clust, N = 1L, regime = "moment",
                        ridge_eps = 0)
ev <- eigen(fit_pca@covariances[[1L]])$vectors
pr <- prcomp(x_clust)$rotation
## Each direction may point either way, so compare absolute values.
diff_pca <- max(abs(abs(ev) - abs(unname(pr))))

## ----fig-pca, eval = has_ggplot2, echo = has_ggplot2, fig.cap = "Clusters (colour) from a three-component fit and principal directions (arrows) from a one-component fit to the same data.", fig.alt = "Three coloured point clusters with two principal-axis arrows drawn from the overall centre of the data."----
vals <- eigen(fit_pca@covariances[[1L]])$values
mu_pca <- unname(fit_pca@means[[1L]])
axes <- data.frame(
  x = mu_pca[1L], y = mu_pca[2L],
  xend = mu_pca[1L] + 2 * sqrt(vals) * ev[1L, ],
  yend = mu_pca[2L] + 2 * sqrt(vals) * ev[2L, ]
)
pts <- data.frame(x_clust, cluster = factor(max.col(resp)))
ggplot2::ggplot() +
  ggplot2::geom_point(data = pts,
                      ggplot2::aes(V1, V2, colour = cluster),
                      alpha = 0.6, size = 0.9) +
  ggplot2::geom_segment(
    data = axes, ggplot2::aes(x = x, y = y, xend = xend, yend = yend),
    arrow = grid::arrow(length = grid::unit(0.2, "cm")), linewidth = 0.8
  ) +
  ggplot2::scale_colour_viridis_d(name = "cluster", end = 0.85) +
  ggplot2::coord_equal() +
  ggplot2::labs(x = expression(x[1]), y = expression(x[2]),
                title = "Clusters and principal directions") +
  ggplot2::theme_minimal(base_size = 11)

## ----fig-pca-skip, eval = !has_ggplot2, echo = FALSE, results = "asis"--------
# cat("ggplot2 is not installed, so this figure is skipped.\n")

## ----ridge--------------------------------------------------------------------
lambda <- c(0, 0.5, 2, 8)
slope_pm <- vapply(lambda, function(lam) {
  f <- fit_proxymix(joint, N = 1L, regime = "moment", ridge_eps = lam)
  gmm_conditionalise(f, given = c(NA, 1))@means[[1L]] -
    gmm_conditionalise(f, given = c(NA, 0))@means[[1L]]
}, numeric(1L))
slope_formula <- cov(x, y) / (var(x) + lambda)
diff_ridge <- max(abs(slope_pm - slope_formula))

## ----ridge-table, echo = FALSE------------------------------------------------
knitr::kable(
  data.frame(lambda = lambda, proxymix = slope_pm,
             ridge_formula = slope_formula),
  digits = 4L,
  col.names = c("Penalty (lambda)", "proxymix slope",
                "Ridge formula"),
  caption = paste0(
    "The conditional slope after adding lambda to the variances, and the ",
    "ridge estimate cov(x, y) / (var(x) + lambda)."
  )
)

## ----uplift-fit---------------------------------------------------------------
set.seed(20260902)
n_up <- 600L
x_up <- rnorm(n_up)
t_up <- rbinom(n_up, 1L, 0.5)
tau_true <- function(v) 0.5 + v
y_up <- 1 + tau_true(x_up) * t_up + rnorm(n_up, sd = 0.5)
dat_up <- data.frame(y = y_up, t = t_up, x = x_up)

model <- fit_uplift(dat_up, "y", "t", "x", N = 2L, regime = "sample",
                    max_iter = 80L, seed = 1L)
model

## ----uplift-cate--------------------------------------------------------------
grid_up <- data.frame(x = seq(-2, 2, length.out = 41L))
cate <- proxy_cate(model, grid_up)
err_cate <- max(abs(cate$tau - tau_true(grid_up$x)))

## T-learner: one least-squares line per arm, fitted as one model.
fit_arms <- lm(y ~ t * x, data = dat_up)
b_arms <- coef(fit_arms)[c("t", "t:x")]
diff_arms <- max(abs(cate$tau - (b_arms[[1L]] + b_arms[[2L]] * grid_up$x)))
a_arms <- cbind(0, 1, 0, grid_up$x)          # picks out t + x * t:x
se_arms <- sqrt(rowSums((a_arms %*% vcov(fit_arms)) * a_arms))
se_ratio <- range(cate$se / se_arms)
## Distance of the fitted intercept and slope from 0.5 and 1, in
## standard errors.
z_arms <- (b_arms - c(0.5, 1)) / sqrt(diag(vcov(fit_arms)))[c("t", "t:x")]
## Treatment value at the centre of each mixture component.
arm_of_component <- vapply(model@fit@means, function(m) {
  m[model@roles$treatment]
}, numeric(1L))

## ----fig-cate, eval = has_ggplot2, echo = has_ggplot2, fig.cap = "The estimated treatment effect with its 95 per cent interval, and the true effect 0.5 + x. The estimate is a straight line, the difference between the two arms' least-squares lines.", fig.alt = "Estimated treatment effect against the covariate, shown as a straight line with a shaded interval band, and a dashed line for the true effect; the two lines nearly coincide."----
cate_df <- as.data.frame(cate)
cate_df$x <- grid_up$x
cate_df$truth <- tau_true(grid_up$x)
ggplot2::ggplot(cate_df, ggplot2::aes(x)) +
  ggplot2::geom_ribbon(ggplot2::aes(ymin = ci_lo, ymax = ci_hi),
                       fill = "#56B4E9", alpha = 0.35) +
  ggplot2::geom_line(ggplot2::aes(y = tau, colour = "proxymix estimate"),
                     linewidth = 0.9) +
  ggplot2::geom_line(ggplot2::aes(y = truth, colour = "true effect"),
                     linewidth = 0.9, linetype = "dashed") +
  ggplot2::geom_hline(yintercept = 0, colour = "grey60",
                      linewidth = 0.3) +
  ggplot2::scale_colour_manual(
    name = NULL,
    values = c("proxymix estimate" = "#0072B2",
               "true effect" = "#D55E00")
  ) +
  ggplot2::labs(
    x = "covariate x", y = "treatment effect on y",
    title = "Treatment effect from one mixture fit"
  ) +
  ggplot2::theme_minimal(base_size = 11) +
  ggplot2::theme(legend.position = "top")

## ----fig-cate-skip, eval = !has_ggplot2, echo = FALSE, results = "asis"-------
# cat("ggplot2 is not installed, so this figure is skipped.\n")

## ----uplift-decide------------------------------------------------------------
decision <- proxy_decide(model, grid_up, value = 1, cost = 0.5)
switch_x <- grid_up$x[min(which(decision$action == 1L))]

## ----uplift-decide-table, echo = FALSE----------------------------------------
sel <- c(1L, 11L, 21L, 31L, 41L)
knitr::kable(
  data.frame(
    x = grid_up$x[sel],
    tau = cate$tau[sel],
    truth = tau_true(grid_up$x[sel]),
    action = decision$action[sel],
    expected_value = decision$expected_value[sel]
  ),
  digits = 3L,
  col.names = c("Covariate x", "Estimated effect", "True effect",
                "Recommended arm", "Net value"),
  caption = paste0(
    "Recommendations at five values of x, for a value of 1 per unit of ",
    "outcome and a treatment cost of 0.5. Treating pays when the effect ",
    "exceeds 0.5. The net value of treating is the estimated effect times ",
    "the value, minus the cost."
  )
)

## ----uplift-refusal-----------------------------------------------------------
refusal <- tryCatch(
  proxy_cate(model, grid_up, t1 = 100, t0 = 0),
  error = function(e) e
)
msg_lines <- strsplit(conditionMessage(refusal), "\n")[[1L]]
writeLines(strwrap(msg_lines, width = 70L, exdent = 2L))

## ----uplift-report------------------------------------------------------------
proxy_identification_report(model, grid_up)

## ----equality-table, echo = FALSE---------------------------------------------
knitr::kable(
  data.frame(
    claim = c(
      "One-component conditional slope equals the lm slope",
      "Per-point kernel estimate, conditioned, equals Nadaraya-Watson",
      "One-component eigenvectors equal the prcomp directions",
      "Conditional slope with ridge_eps equals the ridge formula",
      "Two-component treatment effect equals the per-arm lm contrast"
    ),
    difference = sci_plain(c(diff_reg, diff_nw, diff_pca, diff_ridge,
                             diff_arms)),
    stringsAsFactors = FALSE
  ),
  col.names = c("Claim", "Largest absolute difference"),
  caption = paste0(
    "Each mixture result compared with the usual tool, on the data ",
    "fitted above."
  )
)

## ----summary-table, echo = FALSE----------------------------------------------
knitr::kable(
  data.frame(
    method = c("Regression", "Kernel regression", "Clustering",
               "Principal components", "Ridge", "Treatment effects"),
    usual = c("lm, glm", "ksmooth, np", "kmeans, mclust", "prcomp",
              "glmnet, lm.ridge", "T-learner, grf, DoubleML"),
    route = c("mixture over (y, x), then gmm_conditionalise()",
              "one component per point, then gmm_conditionalise()",
              "fit_proxymix(regime = \"sample\")",
              "eigen() of the one-component covariance",
              "ridge_eps",
              "fit_uplift(), then proxy_cate()"),
    gain = c("curved means; full conditional distribution",
             "full conditional distribution; works from a formula alone",
             "elliptical clusters with soft assignments",
             "directions within each cluster",
             "shrinkage from the same fit",
             "one fit for all queries; a stated list of assumptions"),
    give_up = c("standard errors; normal components",
                "cost grows with n unless compressed; bandwidth choice",
                "number of clusters; speed on very large data",
                "loadings, scree plots and biplots",
                "lasso and variable selection",
                "per-unit accuracy when an arm needs several components"),
    stringsAsFactors = FALSE
  ),
  col.names = c("Analysis", "Usual tool", "proxymix route", "Gain",
                "Cost"),
  caption = "Six analyses from fitted mixtures, and what each gains and costs."
)

## ----session-info, collapse = FALSE, class.output = "session-info"------------
sessionInfo()

