# ---------------------------------------------------------------------------
# Themes.
#
# The house style: 9 pt type, black and unbolded everywhere including the panel
# labels, hairline rules, and no chartjunk between the panels. Colour is
# carried by the tracks themselves and never by the type, which is what keeps a
# stack of six panels readable at column width.
# ---------------------------------------------------------------------------

#' The omakase theme
#'
#' A blank canvas for genome-track panels: no axes, no grid, no legend, tight
#' margins, and every piece of text black, unbolded and at the base size. Track
#' panels sit directly on top of one another, so any horizontal rule or panel
#' border between them reads as a feature of the data.
#'
#' @param base_size Base font size in points. Everything else is expressed
#'   relative to it.
#' @param base_family Font family. The empty string uses the device default.
#' @param margin Plot margin in points, given as a length-4 numeric in the
#'   order top, right, bottom, left.
#' @param legend Legend position; `"none"` by default.
#' @return A [ggplot2::theme()] object.
#' @examples
#' library(ggplot2)
#' ggplot(mtcars, aes(wt, mpg)) + geom_point() + theme_omakase()
#' @export
theme_omakase <- function(base_size = 9, base_family = "",
                          margin = c(1, 4, 1, 4), legend = "none") {
  ggplot2::theme_void(base_size = base_size, base_family = base_family) +
    ggplot2::theme(
      text = ggplot2::element_text(size = base_size, colour = "black",
                                   face = "plain", family = base_family),
      plot.margin = ggplot2::margin(margin[1], margin[2], margin[3], margin[4]),
      plot.title = ggplot2::element_text(size = base_size, colour = "black",
                                         face = "plain", hjust = 0),
      legend.position = legend,
      legend.title = ggplot2::element_text(size = base_size, colour = "black"),
      legend.text = ggplot2::element_text(size = base_size, colour = "black"),
      legend.key = ggplot2::element_blank(),
      legend.background = ggplot2::element_blank()
    )
}

#' A themed variant that keeps the axes
#'
#' Same type and weight as [theme_omakase()], but with a visible x axis and y
#' axis. Use it when a panel is being read quantitatively - a coverage track
#' whose depth matters, or a standalone plot outside a sashimi stack.
#'
#' @inheritParams theme_omakase
#' @param hairline Line width for the axis rules.
#' @param grid Whether to draw a light horizontal grid.
#' @return A [ggplot2::theme()] object.
#' @examples
#' library(ggplot2)
#' ggplot(mtcars, aes(wt, mpg)) + geom_point() + theme_omakase_axes()
#' @export
theme_omakase_axes <- function(base_size = 9, base_family = "",
                               hairline = 0.3, margin = c(1, 4, 1, 4),
                               legend = "none", grid = FALSE) {
  th <- ggplot2::theme_classic(base_size = base_size, base_family = base_family) +
    ggplot2::theme(
      text = ggplot2::element_text(size = base_size, colour = "black",
                                   face = "plain", family = base_family),
      axis.text = ggplot2::element_text(size = base_size, colour = "black"),
      axis.title = ggplot2::element_text(size = base_size, colour = "black"),
      axis.line = ggplot2::element_line(linewidth = hairline, colour = "black"),
      axis.ticks = ggplot2::element_line(linewidth = hairline, colour = "black"),
      plot.title = ggplot2::element_text(size = base_size, colour = "black",
                                         face = "plain", hjust = 0),
      plot.margin = ggplot2::margin(margin[1], margin[2], margin[3], margin[4]),
      legend.position = legend,
      legend.key = ggplot2::element_blank()
    )
  if (grid) {
    th <- th + ggplot2::theme(
      panel.grid.major.y = ggplot2::element_line(linewidth = hairline / 2,
                                                 colour = "#E6E6E6")
    )
  }
  th
}

# ggplot2 measures geom_text/geom_label in millimetres while themes use points,
# so a label written at `size = base_size` is nearly three times too big. This
# is the conversion every text layer in the package goes through.
#' @noRd
pt_to_mm <- function(pt) pt / .pt

# Resolve the `theme` argument, which accepts a theme object, a name, or NULL.
#' @noRd
resolve_theme <- function(theme, base_size, base_family = "", margin = c(1, 4, 1, 4)) {
  if (inherits(theme, "theme")) return(theme)
  if (is.null(theme) || identical(theme, "omakase")) {
    return(theme_omakase(base_size, base_family, margin))
  }
  if (identical(theme, "axes")) {
    return(theme_omakase_axes(base_size, base_family, margin = margin))
  }
  if (is.function(theme)) return(theme(base_size))
  om_abort("{.arg theme} must be a ggplot2 theme, {.val omakase}, {.val axes}, or a function.")
}
