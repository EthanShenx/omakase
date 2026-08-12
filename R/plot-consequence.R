# ---------------------------------------------------------------------------
# Consequence composition figures.
#
# The donut is the reference design: equal-angle radial bars whose radius
# encodes abundance, dotted magnitude guides, the percentage inside the ring,
# external labels on light leader lines, and colour carrying the category. It
# is drawn here in ggplot2 rather than matplotlib so it lives in the same
# graphics system as everything else in the package.
# ---------------------------------------------------------------------------

#' Draw the consequence composition
#'
#' Visualises how a set of start-site switches divides among consequence
#' classes. Four styles are available; the donut is the default and is the one
#' the package's house design was built around.
#'
#' @section Geometry of the donut:
#' The \eqn{n} subtypes each occupy an equal angular sector of
#' \eqn{(360 - n g)/n} degrees, separated by a gap of \eqn{g}. A sector's inner
#' radius is fixed at `inner_r`; its outer radius is
#' \eqn{r_i = r_0 + h_i} with the bar height \eqn{h_i} scaled linearly from the
#' counts,
#' \deqn{h_i = h_{\min} + \frac{n_i - \min n}{\max n - \min n}(h_{\max} - h_{\min}).}
#' So the *radius* encodes abundance while the angle is constant - which is the
#' point of the design, since equal angles keep every label legible regardless
#' of how rare its category is.
#'
#' @param x A data frame from [classify_consequence()] or
#'   [consequence_summary()]. A raw classified table is summarised
#'   automatically.
#' @param style `"donut"`, `"bar"`, `"lollipop"` or `"stacked"`.
#' @param by `"subtype"` or `"category"`.
#' @param filter An optional expression evaluated in `x` to subset rows, for
#'   example `both_full_length == 1`.
#' @param palette Named vector of category colours, or `NULL` for the default.
#' @param title Plot title, or `NULL`.
#' @param base_size Base font size in points.
#' @param base_family Font family for every piece of text in the figure. The
#'   empty string uses the graphics device's default. Text geoms do not inherit
#'   a theme's family, so this is threaded to each of them explicitly.
#' @param inner_r Inner radius of the ring, in the plot's arbitrary units.
#' @param bar_range Length-2 numeric giving the minimum and maximum bar height.
#' @param gap_deg Angular gap between sectors, in degrees.
#' @param start_deg Angle at which the first sector begins; 90 puts it at the
#'   top.
#' @param show_percent Print each sector's share inside the ring.
#' @param show_n Append the raw count to each external label.
#' @param show_guides Draw the dotted magnitude guide circles.
#' @param legend_position Where to put the category legend.
#'
#' @return A `ggplot` object.
#'
#' @examples
#' d <- data.frame(
#'   category = c(rep("5'UTR change", 3), "Promoter swap", rep("N-terminal/CDS", 2)),
#'   subtype = c("longer", "shorter", "equal", "alt first exon",
#'               "N-term extension", "ORF loss"),
#'   n = c(48, 31, 12, 22, 15, 8)
#' )
#' plot_consequence(d)
#' plot_consequence(d, style = "lollipop")
#'
#' @export
plot_consequence <- function(x, style = c("donut", "bar", "lollipop", "stacked"),
                             by = c("subtype", "category"), filter = NULL,
                             palette = NULL, title = NULL, base_size = 10,
                             base_family = "", inner_r = 0.24, bar_range = c(0.11, 0.185),
                             gap_deg = 2, start_deg = 90,
                             show_percent = TRUE, show_n = FALSE,
                             show_guides = TRUE,
                             legend_position = "bottom") {
  style <- match.arg(style)
  by <- match.arg(by)

  d <- as_df(x)
  if (!has_col(d, "n")) {
    d <- consequence_summary(d, by = by, filter = !!rlang::enquo(filter))
  } else {
    if (!has_col(d, "label")) {
      d$label <- if (by == "subtype") {
        lab <- unname(consequence_labels()[as.character(d$subtype)])
        ifelse(is.na(lab), as.character(d$subtype), lab)
      } else {
        as.character(d$category)
      }
    }
    if (!has_col(d, "proportion")) d$proportion <- 100 * d$n / sum(d$n)
  }
  d$category <- as.character(d$category)
  d$label <- as.character(d$label)
  if (show_n) d$label <- sprintf("%s (%d)", d$label, as.integer(d$n))

  cols <- palette %||% OM_CONSEQUENCE_FILL
  missing <- setdiff(unique(d$category), names(cols))
  if (length(missing)) {
    extra <- omakase_palette("omakase", length(missing))
    cols <- c(cols, stats::setNames(extra, missing))
  }

  p <- switch(style,
    donut = consequence_donut(d, cols, base_size, base_family, inner_r,
                              bar_range, gap_deg, start_deg, show_percent,
                              show_guides),
    bar = consequence_bar(d, cols, base_size, base_family, show_percent),
    lollipop = consequence_lollipop(d, cols, base_size, base_family,
                                    show_percent),
    stacked = consequence_stacked(d, cols, base_size, base_family,
                                  show_percent)
  )

  p <- p +
    ggplot2::scale_fill_manual(values = cols, name = "Consequence type",
                               breaks = intersect(consequence_levels()$categories,
                                                  names(cols))) +
    ggplot2::scale_colour_manual(values = cols, guide = "none")

  if (!is.null(title)) p <- p + ggplot2::ggtitle(title)
  p + ggplot2::theme(legend.position = legend_position)
}

#' @noRd
consequence_donut <- function(d, cols, base_size, base_family, inner_r,
                              bar_range, gap_deg, start_deg, show_percent,
                              show_guides) {
  n <- nrow(d)
  seg <- (360 - n * gap_deg) / n
  txt <- pt_to_mm(base_size)

  d$h <- scale_bar_heights(d$n, bar_range)
  d$r_out <- inner_r + d$h

  # Sectors run clockwise from start_deg, which is how the reference figure
  # reads.
  starts <- start_deg - (seq_len(n) - 1) * (seg + gap_deg)
  d$theta_hi <- starts
  d$theta_lo <- starts - seg
  d$mid <- (d$theta_hi + d$theta_lo) / 2
  rad <- d$mid * pi / 180
  d$ux <- cos(rad)
  d$uy <- sin(rad)

  # ggforce measures arcs clockwise from twelve o'clock, whereas the angles
  # above are the usual anticlockwise-from-east convention.
  to_arc <- function(deg) (90 - deg) * pi / 180

  leader_r <- max(d$r_out) + 0.055
  label_r <- leader_r + 0.105
  # Labels run outward horizontally but are only one line tall, so a square
  # box would leave a deep empty band above and below the ring. The aspect
  # ratio stays fixed, so the ring is still a circle.
  lim_x <- label_r + 0.42
  lim_y <- label_r + 0.14

  guides <- if (show_guides) {
    data.frame(r = unique(c(min(d$r_out), stats::median(d$r_out),
                            max(d$r_out))))
  } else {
    data.frame(r = numeric(0))
  }

  leaders <- do.call(rbind, lapply(seq_len(n), function(i) {
    data.frame(
      x = c(d$r_out[i], d$r_out[i] + 0.10, leader_r) * d$ux[i],
      y = c(d$r_out[i], d$r_out[i] + 0.10, leader_r) * d$uy[i],
      g = i
    )
  }))
  ticks <- data.frame(
    x = (inner_r - 0.020) * d$ux, xend = inner_r * d$ux,
    y = (inner_r - 0.020) * d$uy, yend = inner_r * d$uy
  )

  ang <- d$mid %% 360
  d$hjust <- ifelse(ang >= 45 & ang < 135, 0.5,
             ifelse(ang >= 135 & ang < 225, 1,
             ifelse(ang >= 225 & ang < 315, 0.5, 0)))
  d$vjust <- ifelse(ang >= 45 & ang < 135, 0,
             ifelse(ang >= 135 & ang < 225, 0.5,
             ifelse(ang >= 225 & ang < 315, 1, 0.5)))
  d$lx <- label_r * d$ux
  d$ly <- label_r * d$uy
  d$px <- inner_r * 0.72 * d$ux
  d$py <- inner_r * 0.72 * d$uy

  p <- ggplot2::ggplot()

  if (nrow(guides)) {
    p <- p + ggforce::geom_circle(
      data = guides, ggplot2::aes(x0 = 0, y0 = 0, r = .data$r),
      linetype = "dotted", colour = "#CCCCCC", linewidth = 0.25,
      inherit.aes = FALSE
    )
  }

  p +
    ggplot2::geom_path(
      data = leaders,
      ggplot2::aes(x = .data$x, y = .data$y, group = .data$g),
      colour = "#999999", linewidth = 0.2, lineend = "round"
    ) +
    ggforce::geom_arc_bar(
      data = d,
      ggplot2::aes(x0 = 0, y0 = 0, r0 = inner_r, r = .data$r_out,
                   start = to_arc(.data$theta_hi), end = to_arc(.data$theta_lo),
                   fill = .data$category),
      colour = "white", linewidth = 0.3
    ) +
    ggforce::geom_circle(
      data = data.frame(r = inner_r), ggplot2::aes(x0 = 0, y0 = 0, r = .data$r),
      fill = "white", colour = "black", linewidth = 0.4, inherit.aes = FALSE
    ) +
    ggplot2::geom_segment(
      data = ticks,
      ggplot2::aes(x = .data$x, xend = .data$xend, y = .data$y,
                   yend = .data$yend),
      colour = "black", linewidth = 0.35, lineend = "round"
    ) +
    {
      if (show_percent) {
        ggplot2::geom_text(
          data = d,
          ggplot2::aes(x = .data$px, y = .data$py,
                       label = format_percent(.data$proportion)),
          size = txt * 0.7, colour = "black", family = base_family
        )
      }
    } +
    ggplot2::geom_text(
      data = d,
      ggplot2::aes(x = .data$lx, y = .data$ly, label = .data$label,
                   hjust = .data$hjust, vjust = .data$vjust),
      size = txt * 0.8, colour = "black", family = base_family
    ) +
    ggplot2::coord_fixed(xlim = c(-lim_x, lim_x), ylim = c(-lim_y, lim_y)) +
    theme_omakase(base_size, base_family, legend = "bottom") +
    ggplot2::theme(
      legend.title = ggplot2::element_text(size = base_size * 0.8),
      legend.text = ggplot2::element_text(size = base_size * 0.8)
    )
}

#' @noRd
consequence_bar <- function(d, cols, base_size, base_family, show_percent) {
  txt <- pt_to_mm(base_size)
  d$label <- factor(d$label, levels = rev(d$label))
  p <- ggplot2::ggplot(d, ggplot2::aes(x = .data$n, y = .data$label,
                                       fill = .data$category)) +
    ggplot2::geom_col(width = 0.72) +
    ggplot2::labs(x = "Switches", y = NULL)
  if (show_percent) {
    p <- p + ggplot2::geom_text(
      ggplot2::aes(label = format_percent(.data$proportion)),
      hjust = -0.15, size = txt * 0.8, colour = "black", family = base_family
    ) + ggplot2::scale_x_continuous(expand = ggplot2::expansion(c(0, 0.12)))
  }
  p + theme_omakase_axes(base_size, base_family, legend = "bottom")
}

#' @noRd
consequence_lollipop <- function(d, cols, base_size, base_family, show_percent) {
  txt <- pt_to_mm(base_size)
  d$label <- factor(d$label, levels = rev(d$label))
  p <- ggplot2::ggplot(d, ggplot2::aes(y = .data$label)) +
    ggplot2::geom_segment(
      ggplot2::aes(x = 0, xend = .data$n, yend = .data$label,
                   colour = .data$category),
      linewidth = 0.5
    ) +
    ggplot2::geom_point(
      ggplot2::aes(x = .data$n, fill = .data$category),
      shape = 21, size = 3.2, colour = "white", stroke = 0.5
    ) +
    ggplot2::labs(x = "Switches", y = NULL)
  if (show_percent) {
    p <- p + ggplot2::geom_text(
      ggplot2::aes(x = .data$n, label = format_percent(.data$proportion)),
      hjust = -0.4, size = txt * 0.8, colour = "black", family = base_family
    ) + ggplot2::scale_x_continuous(expand = ggplot2::expansion(c(0.02, 0.16)))
  }
  p + theme_omakase_axes(base_size, base_family, legend = "bottom")
}

#' @noRd
consequence_stacked <- function(d, cols, base_size, base_family, show_percent) {
  txt <- pt_to_mm(base_size)
  d$label <- factor(d$label, levels = d$label)
  p <- ggplot2::ggplot(d, ggplot2::aes(x = 1, y = .data$proportion,
                                       fill = .data$category)) +
    ggplot2::geom_col(width = 0.6, colour = "white", linewidth = 0.3) +
    ggplot2::labs(x = NULL, y = "Share of switches (%)")
  if (show_percent) {
    p <- p + ggplot2::geom_text(
      ggplot2::aes(label = .data$label),
      position = ggplot2::position_stack(vjust = 0.5),
      size = txt * 0.7, colour = "black", family = base_family
    )
  }
  p + ggplot2::scale_x_continuous(breaks = NULL) +
    theme_omakase_axes(base_size, base_family, legend = "bottom")
}

# Linear scaling of counts onto the template's radial bar heights. When every
# count is the same, the mid-height is used - a ring of minimum-height bars
# would read as "no data".
#' @noRd
scale_bar_heights <- function(n, range = c(0.11, 0.185)) {
  v <- as.numeric(n)
  lo <- min(v, na.rm = TRUE)
  hi <- max(v, na.rm = TRUE)
  if (!is.finite(lo) || !is.finite(hi) || isTRUE(all.equal(lo, hi))) {
    return(rep(mean(range), length(v)))
  }
  range[1] + (v - lo) / (hi - lo) * (range[2] - range[1])
}

#' Format a percentage
#'
#' Whole numbers when the value is exact, one decimal otherwise, so a donut of
#' round percentages is not littered with trailing zeros.
#'
#' @param v Numeric vector of percentages.
#' @return A character vector.
#' @examples
#' format_percent(c(25, 12.5, 33.333))
#' @export
format_percent <- function(v) {
  ifelse(is.na(v), NA_character_,
         ifelse(abs(v - round(v)) < 1e-8, sprintf("%.0f%%", v),
                sprintf("%.1f%%", v)))
}
