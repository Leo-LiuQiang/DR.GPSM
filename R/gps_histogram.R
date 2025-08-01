#' Plot overlaid histograms of Generalized Propensity Scores
#'
#' @param data    Data frame returned by \code{gps_pre_process()} (must contain \code{gps_} columns)
#' @param bins    Number of histogram bins (default = 30)
#' @param palette Character vector of fill colors (recycled if shorter than the number of GPS columns)
#' @param alpha   Fill transparency (0–1, default = 0.4) for overlapping areas
#' @param eps     Small value to truncate probabilities (avoids Inf after logit; default = 1e-6)
#'
#' @return A ggplot object showing overlaid histograms and density curves of GPS distributions
#' @export
#'
#' @importFrom ggplot2 ggsave
#' @importFrom rlang .data
#' @importFrom ggplot2 aes after_stat
gps_histogram <- function(data,
                          bins    = 30,
                          palette = NULL,
                          alpha   = 0.4,
                          eps     = 1e-6) {

  gps_cols <- grep("^gps_", names(data), value = TRUE)
  if (length(gps_cols) < 2)
    stop("Need two or more gps_* columns to plot.")

  df_long <- tidyr::pivot_longer(
    data[, gps_cols],
    cols      = tidyselect::everything(),
    names_to  = "level",
    values_to = "prob"
  ) |>
    dplyr::mutate(
      level      = sub("^gps_", "", .data$level),
      prob_clip  = pmin(pmax(.data$prob, eps), 1 - eps),
      logit_prob = stats::qlogis(.data$prob_clip)
    )

  n_levels <- length(unique(df_long$level))
  if (is.null(palette))
    palette <- scales::hue_pal()(n_levels)
  if (length(palette) < n_levels)
    palette <- rep(palette, length.out = n_levels)

  p <- ggplot2::ggplot(
    df_long,
    ggplot2::aes(x = .data$logit_prob,
                 fill  = .data$level,
                 colour = .data$level)) +

    ggplot2::geom_histogram(
      bins      = bins,
      alpha     = alpha,
      position  = "identity",
      aes(y = after_stat(.data$density))) +

    ggplot2::geom_density(
      fill        = NA,
      linewidth   = 0.9,
      linetype    = "dashed",
      adjust      = 1.2,
      show.legend = FALSE) +

    ggplot2::scale_fill_manual(values = palette, name = "Treatment") +
    ggplot2::scale_colour_manual(values = palette, guide = "none") +
    ggplot2::labs(
      x = "logit(gps)",
      y = "Density",
      title = "Distribution of Propensity Score by Treatments") +
    ggplot2::theme_bw()

  return(p)
}
