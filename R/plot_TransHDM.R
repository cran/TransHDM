utils::globalVariables(".data")

#' Visualize TransHDM Mediation Analysis Results
#'
#' S3 plot method for objects of class \code{"TransHDM"}. Produces a bar chart
#' of the overall effects (indirect, direct, total, proportion mediated),
#' individual mediator effects, or p-values.
#'
#' @param x An object of class \code{"TransHDM"}, returned by \code{TransHDM()}.
#' @param type Character string specifying the type of plot:
#'   \code{"overall"} (default) for overall effect estimates,
#'   \code{"mediator"} for individual mediator contributions
#'     (\eqn{\hat\alpha \times \hat\beta}),
#'   \code{"pvalue"} for joint p-values (\eqn{-\log_{10}(\text{ab\_pv})}),
#'   or \code{"alpha_beta"} for a scatter plot of exposure–mediator versus
#'     mediator–outcome effects, coloured by significance.
#' @param fill Bar fill color. Default is \code{"white"}.
#' @param color Bar border color. Default is \code{"black"}.
#' @param linewidth Bar border line width. Default is \code{0.5}.
#' @param size.text Numeric; size of value labels on bars. Default is \code{3.5}.
#' @param show.values Logical; if \code{TRUE} (default), numeric values are
#'   displayed on each bar.
#' @param digits Integer; number of decimal places for displayed values.
#'   Default is \code{3}.
#' @param base_size Numeric; base font size passed to
#'   \code{\link[ggplot2]{theme_classic}}. Default is \code{12}.
#' @param ... Additional arguments (currently unused).
#'
#' @return A \code{\link[ggplot2]{ggplot}} object. The plot is printed
#'   automatically when called interactively.
#'
#' @param label For \code{type = "alpha_beta"} only: \code{"significant"}
#'   (default) labels only significant mediators, \code{"all"} labels every
#'   mediator, \code{"none"} suppresses labels. Requires the
#'   \pkg{ggrepel} package.
#'
#'
#' @examples
#' \donttest{
#' set.seed(123)
#' target_data <- gen_simData_homo(n = 50, p_x = 3, p_m = 20, rho = 0.1)$data
#' source_data <- gen_simData_homo(n = 100, p_x = 3, p_m = 20, rho = 0.1,
#'   source = TRUE, transferable = TRUE)$data
#' result <- TransHDM(
#'   target_data = target_data, source_data = source_data,
#'   Y = "Y", D = "D", M = paste0("M", 1:20), X = paste0("X", 1:3),
#'   transfer = TRUE, ncore = 1, topN = 10
#' )
#' plot(result)
#' plot(result, type = "mediator")
#' plot(result, type = "pvalue")
#' plot(result, type = "alpha_beta")
#' }
#'
#' @method plot TransHDM
#' @export
plot.TransHDM <- function(x,
                          type = c("overall", "mediator", "pvalue", "alpha_beta"),
                          fill = "white",
                          color = "black",
                          linewidth = 0.5,
                          size.text = 3,
                          show.values = TRUE,
                          digits = 3,
                          base_size = 12,
                          label = c("significant", "all", "none"),
                          ...) {

  type  <- match.arg(type)
  label <- match.arg(label)

  if (type == "overall") {
    .plot_effects(
      x, fill, color, linewidth, size.text,
      show.values, digits, base_size
    )
  } else if (type == "mediator") {
    .plot_mediation(
      x, fill, color, linewidth, size.text,
      show.values, digits, base_size
    )
  } else if (type == "pvalue") {
    .plot_pvalue(
      x, fill, color, linewidth, size.text,
      show.values, digits, base_size
    )
  } else {
    .plot_alpha_beta(x, base_size, label)
  }
}


# ---- internal: overall effects bar chart ----
.plot_effects <- function(x, fill, color, linewidth, size.text,
                          show.values, digits, base_size) {

  effects <- x$effects

  if (!is.data.frame(effects) || nrow(effects) != 4L) {
    stop("'effects' component is missing or malformed in the TransHDM object.")
  }

  eff_vec <- effects$estimate
  names(eff_vec) <- effects$effect

  # Separate PE annotation from bar chart (only indirect / direct / total)
  pe_val <- eff_vec[["pe"]]
  pe_str <- formatC(pe_val, format = "f", digits = digits)

  ord <- c("indirect", "direct", "total")
  eff_vec_bar <- eff_vec[ord]
  eff_lab <- c("Indirect", "Direct", "Total")

  df <- data.frame(
    effect   = factor(eff_lab, levels = eff_lab),
    estimate = eff_vec_bar,
    stringsAsFactors = FALSE
  )

  # Y-axis range: force zero baseline, add headroom for value labels
  yr       <- range(df$estimate, na.rm = TRUE)
  yr_span  <- if (diff(yr) > 0) diff(yr) else max(abs(yr), 0.1)
  y_lo     <- min(0, yr[1])
  y_hi     <- yr[2] + yr_span * 0.25

  p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$effect, y = .data$estimate)) +
    ggplot2::geom_col(
      fill = fill, color = color, linewidth = linewidth, width = 0.65
    ) +
    ggplot2::scale_y_continuous(limits = c(y_lo, y_hi), expand = c(0, 0)) +
    ggplot2::coord_cartesian(clip = "off") +
    ggplot2::labs(x = NULL, y = "Estimate") +
    ggplot2::theme_classic(base_size = base_size) +
    ggplot2::theme(
      axis.text.x  = ggplot2::element_text(color = "black"),
      axis.text.y  = ggplot2::element_text(color = "black"),
      axis.ticks.length = ggplot2::unit(0.15, "cm"),
      plot.margin  = ggplot2::margin(5, 5, 5, 20)
    )

  # PE annotation in the top-left corner
  p <- p + ggplot2::annotate(
    "text",
    x     = 0.55,
    y     = y_hi * 0.95,
    label = paste0("P[E] == ", pe_str),
    parse  = TRUE,
    size  = size.text * 0.9,
    hjust = 0,
    vjust = 0,
    color = "black"
  )

  if (show.values) {
    offset <- yr_span * 0.04
    p <- p + ggplot2::geom_text(
      ggplot2::aes(
        y = .data$estimate + offset,
        label = formatC(.data$estimate, format = "f", digits = digits)
      ),
      size = size.text,
      vjust = -0.3
    )
  }

  return(p)
}


# ---- internal: mediation contributions horizontal bar chart ----
.plot_mediation <- function(x, fill, color, linewidth, size.text,
                            show.values, digits, base_size) {

  contrib <- x$contributions

  if (!is.data.frame(contrib)) {
    stop("'contributions' component is missing or malformed in the TransHDM object.")
  }

  if (nrow(contrib) == 0L) {
    p <- ggplot2::ggplot() +
      ggplot2::annotate(
        "text", x = 0, y = 0,
        label = "No mediators selected",
        size = size.text + 1
      ) +
      ggplot2::theme_void()
    return(p)
  }

  # Sort by absolute mediation effect: largest to smallest
  contrib <- contrib[order(abs(contrib$alpha_beta), decreasing = TRUE), , drop = FALSE]
  contrib$mediator <- factor(contrib$mediator, levels = rev(contrib$mediator))

  # Adaptive value-axis range: include zero when all values share the same sign
  xr       <- range(contrib$alpha_beta, na.rm = TRUE)
  xr_span  <- if (diff(xr) > 0) diff(xr) else max(abs(xr), 0.1)

  # Reserve space for value labels: place each label just outside its bar
  offset   <- xr_span * 0.02
  label_y  <- ifelse(
    contrib$alpha_beta >= 0,
    contrib$alpha_beta + offset,
    contrib$alpha_beta - offset
  )

  # Axis limits from the data and the label positions, forcing a zero baseline
  # so bars start from zero; then add a small margin for the text itself
  margin <- xr_span * 0.10
  x_lo <- if (xr[1] >= 0) {0} else {
    min(xr[1], min(label_y, na.rm = TRUE)) - margin
  }
  x_hi <- if (xr[2] <= 0) {0} else {
    max(xr[2], max(label_y, na.rm = TRUE)) + margin
  }

  p <- ggplot2::ggplot(
    contrib,
    ggplot2::aes(x = .data$mediator, y = .data$alpha_beta)
  ) +
    ggplot2::geom_hline(
      yintercept = 0, linetype = "dashed",
      linewidth = 0.35, color = "grey70"
    ) +
    ggplot2::geom_col(
      fill = fill, color = color, linewidth = linewidth, width = 0.65
    ) +
    ggplot2::coord_flip(clip = "off") +
    ggplot2::scale_y_continuous(limits = c(x_lo, x_hi), expand = c(0, 0)) +
    ggplot2::labs(
      x = NULL,
      y = expression(hat(alpha) %*% hat(beta))
    ) +
    ggplot2::theme_classic(base_size = base_size) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(color = "black"),
      axis.text.y = ggplot2::element_text(color = "black"),
      axis.ticks.length = ggplot2::unit(0.15, "cm"),
      plot.margin = ggplot2::margin(5.5, 20, 5.5, 5.5)
    )

  if (show.values) {
    p <- p + ggplot2::geom_text(
      ggplot2::aes(
        y = ifelse(
          .data$alpha_beta >= 0,
          .data$alpha_beta + offset,
          .data$alpha_beta - offset
        ),
        label = formatC(.data$alpha_beta, format = "f", digits = digits),
        hjust = ifelse(.data$alpha_beta >= 0, 0, 1)
      ),
      size = size.text
    )
  }

  return(p)
}


# ---- internal: p-value horizontal bar chart ----
.plot_pvalue <- function(x, fill, color, linewidth, size.text,
                          show.values, digits, base_size) {

  contrib <- x$contributions

  if (!is.data.frame(contrib)) {
    stop("'contributions' component is missing or malformed in the TransHDM object.")
  }

  if (nrow(contrib) == 0L) {
    p <- ggplot2::ggplot() +
      ggplot2::annotate(
        "text", x = 0, y = 0,
        label = "No mediators selected",
        size = size.text + 1
      ) +
      ggplot2::theme_void()
    return(p)
  }

  # Compute -log10 of joint p-value
  contrib$neg_log10_p <- -log10(contrib$ab_pv)

  # Sort by significance: most significant (largest -log10 p) first
  contrib <- contrib[order(contrib$neg_log10_p, decreasing = TRUE), , drop = FALSE]
  contrib$mediator <- factor(contrib$mediator, levels = rev(contrib$mediator))

  # Adaptive value-axis range: reserve space for value labels
  xr       <- range(contrib$neg_log10_p, na.rm = TRUE)
  xr_span  <- if (diff(xr) > 0) diff(xr) else max(abs(xr), 0.1)

  offset   <- xr_span * 0.04
  label_y  <- contrib$neg_log10_p + offset

  x_hi     <- max(xr[2], max(label_y, na.rm = TRUE)) + xr_span * 0.06

  p <- ggplot2::ggplot(
    contrib,
    ggplot2::aes(x = .data$mediator, y = .data$neg_log10_p)
  ) +
    ggplot2::geom_col(
      fill = fill, color = color, linewidth = linewidth, width = 0.65
    ) +
    ggplot2::coord_flip(clip = "off") +
    ggplot2::scale_y_continuous(limits = c(0, x_hi), expand = c(0, 0)) +
    ggplot2::labs(
      x = NULL,
      y = expression(-log[10](p[value]))
    ) +
    ggplot2::theme_classic(base_size = base_size) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(color = "black"),
      axis.text.y = ggplot2::element_text(color = "black"),
      axis.ticks.length = ggplot2::unit(0.15, "cm"),
      plot.margin  = ggplot2::margin(5, 5, 5, 20)
    )

  if (show.values) {
    p <- p + ggplot2::geom_text(
      ggplot2::aes(
        y = .data$neg_log10_p + offset,
        label = formatC(.data$neg_log10_p, format = "f", digits = digits)
      ),
      size = size.text,
      hjust = 0
    )
  }

  return(p)
}


# ---- internal: alpha–beta scatter plot ----
.plot_alpha_beta <- function(x, base_size, label) {

  all_mediators <- x$all_mediators

  if (!is.data.frame(all_mediators)) {
    stop(
      "'all_mediators' component is missing in the TransHDM object. ",
      "Re-run TransHDM() or joint_test() to obtain it."
    )
  }

  all_mediators$significant <- factor(
    all_mediators$significant,
    levels = c(FALSE, TRUE),
    labels = c("No", "Yes")
  )

  p <- ggplot2::ggplot(
    all_mediators,
    ggplot2::aes(x = .data$alpha, y = .data$beta)
  ) +
    ggplot2::geom_hline(
      yintercept = 0, linetype = "dashed",
      linewidth = 0.35, color = "grey70"
    ) +
    ggplot2::geom_vline(
      xintercept = 0, linetype = "dashed",
      linewidth = 0.35, color = "grey70"
    ) +
    ggplot2::geom_point(
      ggplot2::aes(color = .data$significant, shape = .data$significant),
      size = 2.5, alpha = 0.8
    ) +
    ggplot2::scale_color_manual(
      name   = "Significant",
      values = c("No" = "#999999", "Yes" = "#D73027")
    ) +
    ggplot2::scale_shape_manual(
      name   = "Significant",
      values = c("No" = 1, "Yes" = 16)
    ) +
    ggplot2::labs(
      x = expression(hat(alpha) ~ "(exposure-mediator effect)"),
      y = expression(hat(beta) ~ "(mediator-outcome effect)")
    ) +
    ggplot2::theme_classic(base_size = base_size) +
    ggplot2::theme(
      axis.text  = ggplot2::element_text(color = "black"),
      axis.ticks.length = ggplot2::unit(0.15, "cm"),
      legend.position = "top"
    )

  # ---- Mediator name labels via ggrepel ----
  if (label != "none") {
    if (!requireNamespace("ggrepel", quietly = TRUE)) {
      warning(
        "Package 'ggrepel' is required to label mediators. ",
        "Install it with: install.packages(\"ggrepel\")"
      )
    } else {
      label_data <- if (label == "significant") {
        all_mediators[all_mediators$significant == "Yes", , drop = FALSE]
      } else {
        all_mediators
      }

      if (nrow(label_data) > 0L) {
        p <- p + ggrepel::geom_text_repel(
          data = label_data,
          ggplot2::aes(label = .data$mediator),
          size = base_size / 3,
          max.overlaps = 20,
          seed = 1,
          show.legend = FALSE
        )
      }
    }
  }

  return(p)
}
