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
res <- readRDS("results/entropy.rds")
major_minor <- function(v) paste(unlist(package_version(v))[1:2],
                                 collapse = ".")
if (major_minor(res$proxymix_version) !=
    major_minor(as.character(packageVersion("proxymix")))) {
  stop("results/entropy.rds was built under proxymix ",
       res$proxymix_version, ", but this is proxymix ",
       packageVersion("proxymix"), ". Rerun the simulation and ",
       "data-raw/vignette_results/entropy.R.", call. = FALSE)
}

## Small numbers are written as plain decimals rather than in the
## scientific notation that knitr's inline hook would otherwise use.
fixed <- function(v, digits) {
  format(round(v, digits), nsmall = digits, scientific = FALSE)
}
## a small number as a \times 10^{b}, in LaTeX math
sci <- function(v, digits) {
  e <- floor(log10(abs(v)))
  paste0("$", formatC(v / 10^e, format = "f", digits = digits),
         " \\times 10^{", e, "}$")
}

## ----seed---------------------------------------------------------------------
set.seed(20260618)

## ----entropy------------------------------------------------------------------
g <- gmm(
  weights = c(0.5, 0.5),
  means = list(c(-2, 0), c(2, 0)),
  covariances = list(diag(2), diag(2))
)
h2_mixture <- gmm_entropy(g)

sigma_one <- matrix(c(1, 0.3, 0.3, 1), 2L, 2L)
one <- gmm(weights = 1, means = list(c(0, 0)), covariances = list(sigma_one))
h2_closed <- gmm_entropy(one)
h2_analytic <- 0.5 * (2 * log(4 * pi) +
  as.numeric(determinant(sigma_one, logarithm = TRUE)$modulus))
h2_gap <- abs(h2_closed - h2_analytic)

sh <- gmm_entropy(g, order = "shannon", n_mc = 5000L, seed = 1L)
sh_slack <- sh$upper_bound - sh$mc
sh_slack_se <- sh_slack / sh$mc_se

## ----entropy-kable, echo = FALSE----------------------------------------------
knitr::kable(
  data.frame(
    quantity = c(
      "Rényi-2, two-component mixture",
      "Rényi-2, single normal, from the package",
      "Rényi-2, single normal, from the formula",
      "Shannon, two-component mixture, simulation estimate",
      "Shannon, standard error of the estimate",
      "Shannon, exact upper bound"
    ),
    value = formatC(c(h2_mixture, h2_closed, h2_analytic, sh$mc, sh$mc_se,
                      sh$upper_bound), format = "f", digits = 4L)
  ),
  col.names = c("Quantity", "Value (nats)"), align = c("l", "r"),
  caption = paste0(
    "Entropy of the two-component mixture and of the single normal ",
    "distribution. The Shannon estimate uses ", sh$n_mc, " draws."
  )
)

## ----divergence---------------------------------------------------------------
q <- gmm(
  weights = 1, means = list(c(0, 0)), covariances = list(diag(2) * 2)
)
d_cs <- gmm_divergence(g, q)
d_self <- gmm_divergence(g, g)
d_kl <- gmm_divergence(g, q, type = "kl", n_mc = 2000L)

## ----divergence-kable, echo = FALSE-------------------------------------------
knitr::kable(
  data.frame(
    quantity = c(
      "Cauchy-Schwarz, g against q", "Cauchy-Schwarz, g against itself",
      "Kullback-Leibler, simulation estimate",
      "Kullback-Leibler, standard error of the estimate",
      "Kullback-Leibler, Hershey-Olsen approximation"
    ),
    value = formatC(c(d_cs, d_self, d_kl$mc, d_kl$mc_se, d_kl$variational),
                    format = "f", digits = 4L)
  ),
  col.names = c("Quantity", "Value (nats)"), align = c("l", "r"),
  caption = "Two divergences between the same pair of mixtures."
)

## ----mutual-information-------------------------------------------------------
sigma_joint <- matrix(c(1, 0.7, 0.7, 1), 2L, 2L)
joint <- gmm(
  weights = 1, means = list(c(0, 0)), covariances = list(sigma_joint)
)
mi <- gmm_mutual_information(joint, 1L, 2L)

independent <- gmm(
  weights = 1, means = list(c(0, 0)), covariances = list(diag(2))
)
mi_independent <- gmm_mutual_information(independent, 1L, 2L)

given_grid <- rbind(c(NA, 0), c(NA, 1), c(NA, 2))
h_cond <- gmm_conditional_entropy(joint, given = given_grid)

## ----mutual-information-kable, echo = FALSE-----------------------------------
knitr::kable(
  data.frame(
    quantity = c(
      "mutual information, correlation 0.7",
      "mutual information, independent variables",
      "conditional entropy of $x_1$ at $x_2 = 0$",
      "conditional entropy of $x_1$ at $x_2 = 1$",
      "conditional entropy of $x_1$ at $x_2 = 2$"
    ),
    value = formatC(c(mi, mi_independent, h_cond), format = "f",
                    digits = 4L)
  ),
  col.names = c("Quantity", "Value (nats)"), align = c("l", "r"),
  caption = paste(
    "Cauchy-Schwarz mutual information between two variables, and the",
    "Rényi-2 entropy of the first given the second."
  )
)

## ----anneal-fit---------------------------------------------------------------
x_three <- rbind(
  matrix(rnorm(200L), ncol = 2L) +
    matrix(rep(c(-7, -7), each = 100L), ncol = 2L),
  matrix(rnorm(200L), ncol = 2L) +
    matrix(rep(c(7, -7), each = 100L), ncol = 2L),
  matrix(rnorm(200L), ncol = 2L) +
    matrix(rep(c(0, 8), each = 100L), ncol = 2L)
)
tgt_three <- gmm_target_from_samples(x_three)
fit_annealed <- fit_em_samples(tgt_three, N = 3L, anneal = TRUE, seed = 1L)
fit_annealed@diagnostics$annealed

## ----anneal-path--------------------------------------------------------------
path <- gmm_anneal_path(x_three, k_max = 6L, n_steps = 60L, seed = 1L)
k_found <- path$k_selected
t_empirical <- path$first_critical_temperature
t_analytic <- path$t_critical_analytic

## ----anneal-kable, echo = FALSE-----------------------------------------------
knitr::kable(
  data.frame(
    quantity = c(
      "components found", "first critical temperature, recorded",
      "first critical temperature, exact", "cooling steps"
    ),
    value = c(
      formatC(k_found, format = "d"),
      formatC(c(t_empirical, t_analytic), format = "f", digits = 2L),
      formatC(nrow(path$path), format = "d")
    )
  ),
  col.names = c("Quantity", "Value"), align = c("l", "r"),
  caption = "What the cooling found on three well-separated clusters."
)

## ----fig-anneal, eval = has_ggplot2, echo = has_ggplot2, fig.height = 4.4, fig.cap = "Cooling on three well-separated clusters. Upper panel: the number of distinct component centres at each temperature. Lower panel: the free energy. The dashed line marks the temperature at which the first split was recorded, the dotted line the exact critical temperature.", fig.alt = "Two stacked panels against a logarithmic temperature axis: the upper panel is a staircase of the number of distinct component centres, the lower panel a free-energy curve that is flat at the hot end and then falls, both with two nearly coincident vertical lines marking the first critical temperature."----
anneal_df <- rbind(
  data.frame(
    temperature = path$path$temperature,
    value = path$path$n_effective,
    panel = "distinct centres"
  ),
  data.frame(
    temperature = path$path$temperature,
    value = path$path$free_energy,
    panel = "free energy"
  )
)
ggplot2::ggplot(anneal_df, ggplot2::aes(temperature, value)) +
  ggplot2::geom_vline(
    xintercept = t_analytic, linetype = "dotted",
    colour = "#0072B2", linewidth = 0.7
  ) +
  ggplot2::geom_vline(
    xintercept = t_empirical, linetype = "dashed",
    colour = "#D55E00", linewidth = 0.7
  ) +
  # a count recorded at one temperature holds until the next, cooler step
  ggplot2::geom_step(direction = "vh", linewidth = 0.8,
                     colour = "#000000") +
  ggplot2::facet_wrap(~ panel, ncol = 1L, scales = "free_y") +
  ggplot2::scale_x_log10() +
  ggplot2::labs(
    x = "temperature (log scale, cooling from right to left)",
    y = NULL,
    title = "Where the components split as the fit cools"
  ) +
  ggplot2::theme_minimal(base_size = 11)

## ----fig-anneal-skip, eval = !has_ggplot2, echo = FALSE, results = "asis"-----
# cat("ggplot2 is not installed on this build, so the cooling figure",
#     "is skipped.\n")

## ----maxent-------------------------------------------------------------------
me_gauss <- maxent_target(moments = list(mean = c(0, 0), cov = diag(2)))
me_unif <- maxent_target(support = list(lower = c(0, 0), upper = c(1, 1)))
unif_density <- exp(me_unif@log_density(matrix(c(0.5, 0.5), nrow = 1L)))
gauss_density <- exp(me_gauss@log_density(matrix(c(0, 0), nrow = 1L)))

## ----maxent-kable, echo = FALSE-----------------------------------------------
knitr::kable(
  data.frame(
    constraint = c("mean and covariance matrix", "the unit square as range"),
    family = c(me_gauss@metadata$family, me_unif@metadata$family),
    check = formatC(c(gauss_density, unif_density), format = "f",
                    digits = 3L)
  ),
  col.names = c("Constraint", "Family returned", "Density at the centre"),
  align = c("l", "l", "r"),
  caption = paste(
    "The maximum-entropy density under each constraint, with mean zero and",
    "identity covariance matrix for the first. The unit square has area",
    "one, so the uniform density on it is one."
  )
)

## ----icl----------------------------------------------------------------------
x_two <- rbind(
  matrix(rnorm(200L, -4), ncol = 2L),
  matrix(rnorm(200L, 4), ncol = 2L)
)
fit_two <- fit_em_samples(gmm_target_from_samples(x_two), N = 2L, seed = 1L)
crit <- bic_aic(fit_two)

## ----icl-kable, echo = FALSE--------------------------------------------------
knitr::kable(
  data.frame(
    criterion = c("BIC", "AIC", "ICL", "$E_N$", "free parameters"),
    value = c(
      formatC(c(crit$bic, crit$aic, crit$icl), format = "f", digits = 2L),
      sci(crit$classification_entropy, 1L),
      formatC(crit$n_params, format = "d")
    )
  ),
  col.names = c("Criterion", "Value"), align = c("l", "r"),
  caption = paste(
    "Information criteria for a two-component fit to two well-separated",
    "clusters."
  )
)

## ----graph-chain--------------------------------------------------------------
omega <- diag(4L)
for (i1 in seq_len(3L)) {
  omega[i1, i1 + 1L] <- -0.5
  omega[i1 + 1L, i1] <- -0.5
} # ends i1, over the off-diagonal band of the chain precision
g_chain <- gmm(
  weights = 1, means = list(rep(0, 4L)), covariances = list(solve(omega))
)
adj_chain <- gmm_independence_graph(g_chain)
edges_chain <- sum(adj_chain) / 2L

## ----graph-chain-kable, echo = FALSE------------------------------------------
knitr::kable(
  as.data.frame(adj_chain[, ]),
  caption = "Edges of the four-variable chain (1 = edge)."
)

## ----graph-field--------------------------------------------------------------
energy <- function(x_mat) {
  x_mat <- matrix(x_mat, ncol = 3L)
  rowSums((x_mat^2 - 1)^2) -
    0.7 * (x_mat[, 1L] * x_mat[, 2L] + x_mat[, 2L] * x_mat[, 3L])
}
field <- gmm_target(n_dim = 3L, log_density = function(x_mat) -energy(x_mat))
fit_field <- fit_kld_em(
  field,
  N = 8L,
  proposal = proposal_uniform(3L, -3, 3),
  is_size = 6000L,
  anneal = TRUE,
  seed = 1L,
  support_warn = FALSE
)
adj_field <- gmm_independence_graph(fit_field)
edges_field <- sum(adj_field) / 2L

## ----graph-field-facts, include = FALSE---------------------------------------
## the target's partial correlations, by summation over a grid on [-3, 3]^3
g_axis <- seq(-3, 3, length.out = 61L)
g_pts <- as.matrix(expand.grid(g_axis, g_axis, g_axis))
g_w <- exp(-energy(g_pts))
g_w <- g_w / sum(g_w)
g_mean <- colSums(g_pts * g_w)
prec_exact <- solve(crossprod(g_pts * sqrt(g_w)) - tcrossprod(g_mean))
pcor_exact <- -prec_exact / tcrossprod(sqrt(diag(prec_exact)))
pcor_field <- attr(adj_field, "pcor")
q_field <- gmm_fit_quality(fit_field)
## the prose quotes the KL on fresh draws, so the fit must have drawn them
stopifnot(is.finite(q_field$heldout_kld))
field_flagged <- isTRUE(q_field$degenerate) || isFALSE(q_field$converged) ||
  q_field$heldout_kld > 0.3
field_is_chain <- edges_field == 2L && adj_field[1L, 3L] == 0L

## ----graph-field-kable, echo = FALSE------------------------------------------
knitr::kable(
  as.data.frame(adj_field[, ]),
  caption = "Edges of the mixture fitted to the three-variable density (1 = edge)."
)

## ----compare-facts, include = FALSE-------------------------------------------
d2 <- "2-d, three components"
d4 <- "4-d, two components"
e_value <- function(design, quantity, method, what) {
  s1 <- res$entropy_tab$design == design &
    res$entropy_tab$quantity == quantity & res$entropy_tab$method == method
  res$entropy_tab[[what]][s1]
}
k_value <- function(design, selector) {
  res$k_tab$share_true[res$k_tab$design == design &
                         res$k_tab$selector == selector]
}
g_value <- function(method, what) {
  res$graph_tab[[what]][res$graph_tab$method == method]
}
pair_z <- abs(res$pair_tab$diff / res$pair_tab$diff_se)
mi_4d <- res$pair_tab$design == d4 &
  res$pair_tab$quantity == "Shannon mutual information"
k_two <- function(design, selector) {
  res$k_two_tab$share_two[res$k_two_tab$design == design &
                            res$k_two_tab$selector == selector]
}
## pairs with no edge in the six-variable chain of the graph design
n_absent <- choose(6L, 2L) - 5L
t_s <- function(s1) fixed(res$time_secs[[s1]], 3)
## a share of datasets as a percentage
pct <- function(v, digits = 0L) {
  paste(format(round(100 * v, digits), nsmall = digits), "per cent")
}

## ----compare-table, echo = FALSE----------------------------------------------
cmp <- res$pair_tab
cmp$truth <- mapply(function(d, q) {
  res$truth$truth[res$truth$design == d & res$truth$quantity == q]
}, cmp$design, cmp$quantity)
cmp$err_proxymix <- mapply(e_value, cmp$design, cmp$quantity, "proxymix",
                           "abs_error")
cmp$err_fnn <- mapply(e_value, cmp$design, cmp$quantity, cmp$rival,
                      "abs_error")
cmp_tbl <- data.frame(
  design = ifelse(cmp$design == d2, "three components, 2 variables",
                  "two components, 4 variables"),
  quantity = sub("Shannon ", "", cmp$quantity),
  truth = fixed(cmp$truth, 3),
  err_proxymix = fixed(cmp$err_proxymix, 3),
  err_fnn = fixed(cmp$err_fnn, 3),
  diff = paste0(fixed(cmp$diff, 3), " (", fixed(cmp$diff_se, 3), ")"),
  stringsAsFactors = FALSE
)
cmp_tbl <- cmp_tbl[order(cmp_tbl$quantity, decreasing = FALSE), ]
knitr::kable(
  cmp_tbl, row.names = FALSE,
  align = c("l", "l", "r", "r", "r", "r"),
  col.names = c("Design", "Shannon quantity", "True value",
                "Error, proxymix", "Error, FNN",
                "Difference (standard error)"),
  caption = paste0(
    "Mean absolute error in nats over ", res$n_rep, " datasets of ", res$n,
    " rows per design. The difference is the proxymix error minus the FNN ",
    "error, paired by dataset. A negative value favours proxymix."
  )
)

## ----compare-code, eval = FALSE-----------------------------------------------
# library(proxymix)
# library(FNN)
# library(mclust)
# library(pcalg)
# 
# # one dataset of 500 rows from the design with three overlapping components
# w <- c(0.4, 0.35, 0.25)
# mu <- list(c(-2, 0), c(2, 1), c(0, 3))
# sig <- list(matrix(c(1, 0.5, 0.5, 1), 2L),
#             matrix(c(1.5, -0.6, -0.6, 0.8), 2L),
#             diag(c(0.5, 1.2)))
# set.seed(1L)
# z <- sample.int(3L, 500L, replace = TRUE, prob = w)
# x <- matrix(rnorm(500L * 2L), 500L, 2L)
# for (k in 1:3) {
#   x[z == k, ] <- x[z == k, , drop = FALSE] %*% chol(sig[[k]]) +
#     matrix(mu[[k]], sum(z == k), 2L, byrow = TRUE)
# }
# 
# # Shannon entropy and mutual information by nearest neighbours
# entropy(x, k = 10L)[10L]
# mutinfo(x[, 1L, drop = FALSE], x[, 2L, drop = FALSE], k = 10L)
# 
# # the same two quantities from a fitted three-component mixture
# fit <- fit_em_samples(gmm_target_from_samples(x), N = 3L, seed = 1L)
# h_mc <- function(g) {
#   gmm_entropy(g, order = "shannon", n_mc = 5000L, seed = 1L)$mc
# }
# h_mc(fit)
# h_mc(gmm_marginalise(fit, 1L)) + h_mc(gmm_marginalise(fit, 2L)) - h_mc(fit)
# 
# # number of components by mclust's BIC and ICL, one to five offered
# Mclust(x, G = 1:5, verbose = FALSE)$G
# icl_m <- mclustICL(x, G = 1:5, verbose = FALSE)
# which(icl_m == max(icl_m, na.rm = TRUE), arr.ind = TRUE)[1L]
# 
# # one dataset from the graph design: a six-variable chain, two components
# omega <- diag(6L)
# omega[abs(row(omega) - col(omega)) == 1L] <- -0.4
# set.seed(1L)
# z_g <- sample.int(2L, 500L, replace = TRUE, prob = c(0.5, 0.5))
# x_g <- matrix(rnorm(500L * 6L), 500L, 6L) %*% chol(solve(omega))
# x_g[z_g == 2L, 1L] <- x_g[z_g == 2L, 1L] + 4
# 
# # the graph by the PC algorithm, and by proxymix
# pc_fit <- pc(suffStat = list(C = cor(x_g), n = nrow(x_g)),
#              indepTest = gaussCItest, alpha = 0.01, p = ncol(x_g))
# adj <- as(pc_fit@graph, "matrix")
# ((adj + t(adj)) > 0) * 1L
# sel <- select_N(gmm_target_from_samples(x_g), candidates = 1:4, seed = 1L)
# gmm_independence_graph(sel$best_fit)

## ----compare-install, eval = FALSE--------------------------------------------
# install.packages(c("FNN", "mclust", "BiocManager"))
# BiocManager::install(c("graph", "RBGL"))
# install.packages("pcalg")

## ----field-seeds, include = FALSE---------------------------------------------
## refits of the same example over seeds 1 to 10, stored by the builder
seed_count <- function(size, what) sum(res$field_tab[[what]][res$field_tab$is_size == size])
seed_n <- function(size) sum(res$field_tab$is_size == size)
ft <- res$field_tab
## every stored fit drew fresh draws, so the flag used their KL
stopifnot("heldout_kld" %in% names(ft), all(is.finite(ft$heldout_kld)),
          identical(ft$flagged, !ft$converged | ft$degenerate | ft$heldout_kld > 0.3))
n_flag <- sum(ft$flagged)
flag_chain <- sum(ft$flagged & ft$chain)
n_miss <- sum(!ft$chain)
miss_flag <- sum(ft$flagged & !ft$chain)
## TRUE when every flag came from a fit that reached its round limit
flag_by_rounds <- n_flag > 0L && all(!ft$converged[ft$flagged]) &&
  !any(ft$degenerate[ft$flagged]) && all(ft$heldout_kld[ft$flagged] <= 0.3)
miss_text <- if (n_miss == 0L) {
  "No fit missed the chain."
} else if (n_miss == 1L) {
  if (miss_flag == 1L) "The one fit that missed the chain was flagged." else
    "The one fit that missed the chain was not flagged."
} else {
  paste0(miss_flag, " of the ", n_miss, " fits that missed the chain were flagged.")
}

## ----session-info, collapse = FALSE, class.output = "session-info"------------
sessionInfo()

