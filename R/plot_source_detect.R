utils::globalVariables(".data")

#' Visualize Source Detection Results
#'
#' S3 plot method for objects of class \code{"source_detection"}.
#' Produces a bar chart of the \eqn{T}-index for each source dataset, with
#' bars extending from a solid baseline at zero and a horizontal reference
#' line at the transferability threshold. Sources whose \eqn{T}-index falls
#' on or below the threshold are considered transferable.
#'
#' @param x An object of class \code{"source_detection"}, returned by
#'   \code{\link{source_detection}}.
#' @param fill_transferable Fill color for transferable-source bars.
#'   Default is \code{"white"}.
#' @param fill_nontransferable Fill color for non-transferable-source bars.
#'   Default is \code{"white"}.
#' @param color Bar border color. Default is \code{"black"}.
#' @param linewidth Bar border line width. Default is \code{0.5}.
#' @param size.text Numeric; size of value labels on bars. Default is \code{3.5}.
#' @param show.values Logical; if \code{TRUE} (default), \eqn{T}-index values
#'   are displayed on each bar.
#' @param digits Integer; number of decimal places for displayed values.
#'   Default is \code{3}.
#' @param base_size Numeric; base font size passed to
#'   \code{\link[ggplot2]{theme_classic}}. Default is \code{12}.
#' @param threshold.linetype Linetype for the threshold reference line.
#'   Default is \code{"dashed"}.
#' @param threshold.color Color for the threshold reference line.
#'   Default is \code{"gray30"}.
#' @param ... Additional arguments (currently unused).
#'
#' @return A \code{\link[ggplot2]{ggplot}} object.
#'
#' @seealso \code{\link{source_detection}}
#'
#' @examples
#' \donttest{
#' set.seed(123)
#' target_data <- data.frame(
#'   Y  = rnorm(200), D  = rnorm(200),
#'   M1 = rnorm(200), M2 = rnorm(200), X1 = rnorm(200)
#' )
#' source1 <- data.frame(
#'   Y  = rnorm(300), D  = rnorm(300),
#'   M1 = rnorm(300), M2 = rnorm(300), X1 = rnorm(300)
#' )
#' source2 <- data.frame(
#'   Y  = rnorm(250), D  = rnorm(250),
#'   M1 = rnorm(250), M2 = rnorm(250), X1 = rnorm(250)
#' )
#' result <- source_detection(
#'   target_data = target_data,
#'   source_data = list(source1, source2),
#'   Y = "Y", D = "D", M = c("M1", "M2"), X = "X1",
#'   kfold = 5, verbose = FALSE
#' )
#' plot(result)
#' }
#'
#' @method plot source_detection
#' @export
plot.source_detection <- function(x,
                                  fill_transferable = "white",
                                  fill_nontransferable = "white",
                                  color = "black",
                                  linewidth = 0.5,
                                  size.text = 3,
                                  show.values = TRUE,
                                  digits = 3,
                                  base_size = 12,
                                  threshold.linetype = "dashed",
                                  threshold.color = "gray30",
                                  ...) {

  if (!inherits(x, "source_detection")) {
    stop("'x' must be an object of class 'source_detection'.")
  }

  .plot_source_detect(
    x, fill_transferable, fill_nontransferable, color, linewidth,
    size.text, show.values, digits, base_size,
    threshold.linetype, threshold.color
  )
}


# ---- internal: source detection bar chart ----
.plot_source_detect <- function(x,
                                fill_transferable,
                                fill_nontransferable,
                                color,
                                linewidth,
                                size.text,
                                show.values,
                                digits,
                                base_size,
                                threshold.linetype,
                                threshold.color) {

  n_source <- length(x$source.loss)

  if (n_source == 0L) {
    p <- ggplot2::ggplot() +
      ggplot2::annotate(
        "text", x = 0, y = 0,
        label = "No source data provided",
        size = size.text + 1
      ) +
      ggplot2::theme_void()
    return(p)
  }

  transferable <- seq_len(n_source) %in% x$transfer.source.id

  df <- data.frame(
    source       = factor(paste0("S", seq_len(n_source)),
                          levels = paste0("S", seq_len(n_source))),
    T_index      = x$T_index,
    transferable = transferable,
    stringsAsFactors = FALSE
  )

  threshold <- x$threshold

  # Y-axis range straddling zero so the baseline axis sits at y = 0
  yr       <- range(df$T_index, na.rm = TRUE)
  yr_span  <- if (diff(yr) > 0) diff(yr) else max(abs(yr), 0.1)
  y_lo     <- min(0, yr[1]) - yr_span * 0.05
  y_hi     <- max(yr[2] + yr_span * 0.25, threshold + yr_span * 0.05)

  # When both fill colors are identical, use a constant fill (no legend)
  two_fills <- !identical(fill_transferable, fill_nontransferable)

  p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$source, y = .data$T_index))

  # Baseline at y = 0 (solid) — bars extend up (positive) and down (negative)
  p <- p + ggplot2::geom_hline(
    yintercept = 0,
    linewidth = 0.4,
    color = "black"
  )

  # Bars: white fill, black border by default
  if (two_fills) {
    p <- p + ggplot2::geom_col(
      ggplot2::aes(fill = .data$transferable),
      color = color, linewidth = linewidth, width = 0.65
    ) +
      ggplot2::scale_fill_manual(
        values  = c("TRUE" = fill_transferable,
                    "FALSE" = fill_nontransferable),
        labels  = c("TRUE" = "Transferable",
                    "FALSE" = "Not transferable"),
        name    = NULL
      )
  } else {
    p <- p + ggplot2::geom_col(
      fill = fill_transferable,
      color = color, linewidth = linewidth, width = 0.65
    )
  }

  # Threshold reference line drawn after the bars so it stays visible on top
  p <- p + ggplot2::geom_hline(
    yintercept = threshold,
    linetype = threshold.linetype,
    linewidth = 0.6,
    color = threshold.color
  )

  p <- p +
    ggplot2::scale_y_continuous(limits = c(y_lo, y_hi), expand = c(0, 0)) +
    ggplot2::labs(x = NULL, y = expression(italic(T) * "-index")) +
    ggplot2::theme_classic(base_size = base_size) +
    ggplot2::theme(
      axis.text.x  = ggplot2::element_text(color = "black"),
      axis.text.y  = ggplot2::element_text(color = "black"),
      axis.ticks.length = ggplot2::unit(0.15, "cm"),
      plot.margin  = ggplot2::margin(5, 40, 5, 5)
    )

  if (two_fills) {
    p <- p + ggplot2::theme(
      legend.position    = "top",
      legend.justification = "left",
      legend.text        = ggplot2::element_text(size = base_size * 0.8),
      legend.key.size    = ggplot2::unit(0.6, "cm")
    )
  }

  if (show.values) {
    offset <- yr_span * 0.04
    p <- p + ggplot2::geom_text(
      ggplot2::aes(
        y = .data$T_index + offset,
        label = formatC(.data$T_index, format = "f", digits = digits)
      ),
      size = size.text,
      vjust = -0.3
    )
  }

  # Threshold annotation placed to the right of the rightmost bar
  p <- p + ggplot2::annotate(
    "text",
    x     = n_source + 0.6,
    y     = threshold,
    label = paste0("threshold = ", round(threshold, digits)),
    size  = size.text * 0.85,
    hjust = 0,
    vjust = -0.5,
    color = threshold.color
  )

  # Extend x-axis limit so the threshold label is never clipped
  p <- p + ggplot2::expand_limits(x = n_source + 1.8)

  return(p)
}
