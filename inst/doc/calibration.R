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

## ----bimodal------------------------------------------------------------------
set.seed(20260619)

f <- function(v) (v[1]^2 - 4)^2

fit <- from_objective(f, lower = -5, upper = 5, N = 6L,
                      is_size = 3000L, n_steps = 5L, seed = 1L)

modes <- gmm_modes(fit)
sort(round(modes$modes[, 1], 3))

## ----bimodal-ess--------------------------------------------------------------
ess_1d <- ess_summary(fit)
c(ess = round(ess_1d$ess, 1), is_size = ess_1d$is_size,
  ess_relative = round(ess_1d$ess_relative, 3))

## ----himmelblau---------------------------------------------------------------
himmelblau <- function(v) {
  x <- v[1]
  y <- v[2]
  (x * x + y - 11)^2 + (x + y * y - 7)^2
}

fit2 <- from_objective(himmelblau, lower = c(-5, -5), upper = c(5, 5),
                       N = 10L, is_size = 4000L, n_steps = 6L, seed = 1L)

found <- gmm_modes(fit2)
found$n

## ----match--------------------------------------------------------------------
truth <- rbind(c(3, 2), c(-2.805118, 3.131312),
               c(-3.779310, -3.283186), c(3.584428, -1.848126))

pair_dist <- as.matrix(dist(rbind(found$modes, truth)))
n_found <- nrow(found$modes)
pair_dist <- pair_dist[seq_len(n_found), n_found + seq_len(nrow(truth))]
nearest <- apply(pair_dist, 1L, which.min)
gap <- apply(pair_dist, 1L, min)
f_at_mode <- apply(found$modes, 1L, himmelblau)

all_distinct <- length(unique(nearest)) == nrow(truth)
worst_gap <- max(gap)
c(modes_found = found$n, one_per_minimum = all_distinct,
  worst_gap = round(worst_gap, 3))

## ----match-table, echo = FALSE------------------------------------------------
knitr::kable(
  data.frame(
    recovered_x1 = round(found$modes[, 1], 3),
    recovered_x2 = round(found$modes[, 2], 3),
    true_x1 = round(truth[nearest, 1], 3),
    true_x2 = round(truth[nearest, 2], 3),
    distance = round(gap, 3),
    objective = round(f_at_mode, 3),
    height = round(found$density, 3)
  ),
  col.names = c("Mode $x_1$", "Mode $x_2$", "True $x_1$", "True $x_2$",
                "Distance", "Objective at mode", "Height of map"),
  caption = paste(
    "Each mode of the map beside the true minimum nearest to it, with the",
    "distance between them, the objective at the mode (zero at a true",
    "minimum) and the height of the map at the mode."
  )
)

## ----himmelblau-map, eval = has_ggplot2, echo = has_ggplot2, fig.height = 5.4, fig.cap = "The four modes of the map (orange circles) sit on the four true minima of Himmelblau's function (black crosses). The background shows the objective on a log scale, light where it is low.", fig.alt = "Surface of Himmelblau's function on a log scale as a shaded raster with contour lines, the four modes of the map as orange circles and the four true minima as black crosses, each pair close together."----
xs <- seq(-5, 5, length.out = 140L)
ys <- seq(-5, 5, length.out = 140L)
surface <- expand.grid(x1 = xs, x2 = ys)
surface$log_f <- log1p(
  apply(as.matrix(surface[, c("x1", "x2")]), 1L, himmelblau)
)

mode_df <- data.frame(x1 = found$modes[, 1], x2 = found$modes[, 2])
truth_df <- data.frame(x1 = truth[, 1], x2 = truth[, 2])

ggplot2::ggplot(surface, ggplot2::aes(x1, x2)) +
  ggplot2::geom_raster(ggplot2::aes(fill = log_f), interpolate = TRUE) +
  ggplot2::geom_contour(ggplot2::aes(z = log_f), colour = "white",
                        linewidth = 0.2, alpha = 0.6, bins = 8L) +
  ggplot2::geom_point(data = truth_df, shape = 4, size = 4, stroke = 1.4,
                      colour = "#000000") +
  ggplot2::geom_point(data = mode_df, shape = 21, size = 3, stroke = 1,
                      fill = "#D55E00", colour = "#000000") +
  ggplot2::scale_fill_viridis_c(name = "log(1 + f)", option = "mako",
                                direction = -1) +
  ggplot2::coord_equal(expand = FALSE) +
  ggplot2::labs(
    title = "One fit, four basins of Himmelblau's function",
    subtitle = "crosses: true minima; filled circles: modes of the map",
    x = expression(x[1]), y = expression(x[2])
  ) +
  ggplot2::theme_minimal(base_size = 11) +
  ggplot2::theme(plot.title = ggplot2::element_text(face = "bold"),
                 panel.grid = ggplot2::element_blank())

## ----himmelblau-map-skip, eval = !has_ggplot2, echo = FALSE, results = "asis"----
# cat("ggplot2 is not installed on this build, so the figure of the",
#     "objective surface is skipped.\n")

## ----session-info, collapse = FALSE, class.output = "session-info"------------
sessionInfo()

