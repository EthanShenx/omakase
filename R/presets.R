# ---------------------------------------------------------------------------
# Presets.
#
# A preset is a named bundle of plotting arguments. Explicit arguments always
# win, so a preset sets the house style for a kind of figure and the caller
# still adjusts anything they like on top of it.
# ---------------------------------------------------------------------------

OM_PRESETS <- list(
  # The package default: Paired, arcs staggered above the coverage, counts on
  # the curve.
  default = list(),

  # Start-site activity figures. An arc is a pointer from a start site to the
  # body of the transcript it starts and its label is that site's activity, so
  # the numbers are continuous rather than read counts and want adaptive
  # precision. Muted Set2 tracks, a desaturated slate/coral pair for the two
  # model rows, and a right-hand gutter for PSI.
  tss = list(
    palette = "set2",
    role_fill = OM_ROLE_FILL_TSS,
    arc_shape = "sine",
    arc_height_rule = "constant",
    arc_height_frac = c(main = 0.80, ATSS = 1.20, alt = 1.20),
    arc_width_rule = "constant",
    arc_width = 0.5,
    arc_label_format = "activity",
    arc_label_position = "on",
    label_background = TRUE,
    label_padding = 0.6,
    show_psi = TRUE,
    psi_pad = 0.30,
    feature_color = "#E08214",
    chrom_style = "ucsc",
    base_size = 9,
    hairline = 0.3,
    panel_height = 1,
    ann_height = 1.5
  ),

  # Junction-count figures from alignments: integer labels, width tracking
  # support, span-scaled heights so many junctions nest rather than pile up,
  # and no PSI gutter unless one was computed.
  junction = list(
    palette = "omakase",
    arc_shape = "sine",
    arc_height_rule = "auto",
    arc_width_rule = "log",
    arc_width = 0.3,
    arc_label_format = "count",
    arc_label_position = "on",
    show_psi = FALSE,
    psi_pad = 0
  ),

  # Deliberately plain, for figures that will be recoloured downstream or
  # printed in one ink.
  minimal = list(
    palette = "mono",
    arc_shape = "sine",
    arc_width_rule = "constant",
    arc_width = 0.4,
    arc_label_position = "above",
    label_background = FALSE,
    show_psi = FALSE,
    psi_pad = 0,
    show_features = FALSE
  ),

  # Closest to what IGV draws: x-splines hugging the baseline, unboxed counts,
  # a tinted panel behind each track.
  igv = list(
    palette = "paired_dark",
    arc_shape = "xspline",
    arc_height_rule = "span",
    arc_width_rule = "log",
    arc_width = 0.35,
    arc_label_format = "count",
    arc_label_position = "above",
    background = "#F5F5F5",
    show_psi = FALSE,
    psi_pad = 0
  )
)

#' Plotting presets
#'
#' Named bundles of arguments for [plot_sashimi()]. A preset sets the house
#' style for a kind of figure; any argument you pass explicitly overrides it.
#'
#' \describe{
#'   \item{`default`}{ColorBrewer *Paired*, arcs above the coverage, counts on
#'     the curve.}
#'   \item{`tss`}{Start-site activity figures: *Set2* tracks, a slate/coral
#'     pair for the two model rows, staggered arc heights at 0.80 and 1.20 of
#'     the panel, adaptive-precision activity labels and a PSI gutter.}
#'   \item{`junction`}{Junction counts from alignments: integer labels, width
#'     tracking read support, span-scaled heights, no PSI gutter.}
#'   \item{`minimal`}{One ink, unboxed labels, no features.}
#'   \item{`igv`}{X-spline arcs on a tinted panel, in the manner of IGV.}
#' }
#'
#' @return A character vector of preset names.
#' @examples
#' presets()
#' @export
presets <- function() names(OM_PRESETS)

# Merge a preset under the caller's explicit arguments. `explicit` names the
# arguments the caller actually supplied, found with match.call(), so a default
# that happens to equal the preset value does not count as an override.
#' @noRd
apply_preset <- function(preset, args, explicit) {
  if (is.null(preset) || identical(preset, "default")) return(args)
  p <- OM_PRESETS[[rlang::arg_match(preset, presets())]]
  for (nm in names(p)) {
    if (!nm %in% explicit) args[[nm]] <- p[[nm]]
  }
  args
}
