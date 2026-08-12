# ---------------------------------------------------------------------------
# Colour.
#
# The defaults are house style: a muted qualitative set for sample groups, a
# dark slate / coral pair for the reference and alternative transcript models,
# and a single amber for annotation features.
# They are deliberately unsaturated, because a sashimi panel is mostly a large
# filled area and a saturated fill at that size overwhelms everything drawn on
# top of it.
# ---------------------------------------------------------------------------

OM_PALETTES <- list(
  # ColorBrewer "Paired". The default, and the palette the companion genome
  # tracks are drawn in, so a sashimi panel and a track view of the same locus
  # agree on colour. Its light/dark pairs also mean a two-condition figure gets
  # a related pair rather than two unrelated hues.
  omakase = c("#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99", "#E31A1C",
              "#FDBF6F", "#FF7F00", "#CAB2D6", "#6A3D9A", "#FFFF99", "#B15928"),
  # Just the saturated member of each Paired pair, for when adjacent groups are
  # unrelated and the light/dark pairing would imply a link that is not there.
  paired_dark = c("#1F78B4", "#33A02C", "#E31A1C", "#FF7F00", "#6A3D9A",
                  "#B15928"),
  # ColorBrewer "Set2". The muted trio the start-site stage tracks are drawn
  # in; `preset = "tss"` selects it.
  set2 = c("#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854",
           "#FFD92F", "#E5C494", "#B3B3B3"),
  # Cooler, for many groups.
  ocean   = c("#3A6B8F", "#4E9FBF", "#7FC5D6", "#A9DCE0", "#2C4A63",
              "#6BA292", "#9CC5A1", "#C7DFC5"),
  # Warm, for before/after contrasts.
  ember   = c("#B23A48", "#E4713A", "#F2A65A", "#F6D186", "#7C2E3C",
              "#C8553D", "#E8998D", "#F4C3B0"),
  # High-contrast, colour-vision-safe (Okabe-Ito).
  okabe   = c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2",
              "#D55E00", "#CC79A7", "#000000"),
  # Neutral, for print where colour is not available.
  mono    = c("#1A1A1A", "#4D4D4D", "#808080", "#B3B3B3", "#D9D9D9")
)

# Roles: which transcript is the reference and which the alternative. The blue
# and red are Paired's own, which is what the pyGenomeTracks configs for the
# same data use, so a sashimi model row and a genome-track model row are the
# same colour.
OM_ROLE_FILL <- c(
  main = "#1F78B4", alt = "#E31A1C", ATSS = "#E31A1C",
  inclusion = "#1F78B4", skipping = "#E31A1C",
  reference = "#1F78B4", alternative = "#E31A1C"
)

# The desaturated slate/coral pair the start-site sashimi figures use for
# their two model rows, selected by `preset = "tss"`.
OM_ROLE_FILL_TSS <- c(
  main = "#465674", alt = "#F18B89", ATSS = "#F18B89",
  inclusion = "#465674", skipping = "#F18B89",
  reference = "#465674", alternative = "#F18B89"
)

OM_MODEL_LINE <- "#2B3A55"
OM_FEATURE_FILL <- "#E08214"

# Consequence categories. Ordered by how much of the transcript the switch
# actually changes, so the ramp carries meaning: gold where only the 5' UTR
# moves, blue where the promoter changes but the protein does not, and a deep
# rust where the protein itself is altered. The two non-biological categories
# stay grey so they never compete with the three that matter.
OM_CONSEQUENCE_FILL <- c(
  "5'UTR change"   = "#FAC91E",
  "Promoter swap"  = "#1F95CE",
  "N-terminal/CDS" = "#781F0F",
  "Unclassified"   = "#CCCCCC",
  "Distal (excluded)" = "#E8E8E8"
)

#' Colour palettes
#'
#' @param name Palette name: `"omakase"` (the default, ColorBrewer *Paired*),
#'   `"paired_dark"` (its saturated members only), `"set2"` (ColorBrewer
#'   *Set2*), `"ocean"`, `"ember"`, `"okabe"` (colour-vision-safe), or
#'   `"mono"`.
#' @param n Number of colours needed. When `n` exceeds the palette length the
#'   palette is interpolated rather than recycled, so groups stay
#'   distinguishable.
#' @return A character vector of hex colours.
#' @examples
#' omakase_palette("omakase", 3)
#' omakase_palette("okabe", 12)
#' @export
omakase_palette <- function(name = "omakase", n = NULL) {
  if (is.character(name) && length(name) > 1L) {
    # Already an explicit vector of colours.
    pal <- name
  } else {
    name <- rlang::arg_match(name, names(OM_PALETTES))
    pal <- OM_PALETTES[[name]]
  }
  if (is.null(n)) return(pal)
  if (n <= length(pal)) return(pal[seq_len(n)])
  grDevices::colorRampPalette(pal)(n)
}

#' Palette names available
#' @return A character vector of palette names.
#' @examples
#' omakase_palettes()
#' @export
omakase_palettes <- function() names(OM_PALETTES)

#' Read a palette file
#'
#' Reads a one-colour-per-line file, the same format `ggsashimi` accepts, so an
#' existing palette file works unchanged. R colour names and hex values are
#' both valid, and only the first column is read.
#'
#' @param path Path to the file, or `NULL` to return `NULL`.
#' @return A character vector of colours, or `NULL`.
#' @examples
#' p <- tempfile()
#' writeLines(c("orange", "cornflowerblue", "#008000"), p)
#' read_palette(p)
#' @export
read_palette <- function(path) {
  if (is.null(path)) return(NULL)
  if (!file.exists(path)) om_abort("Palette file not found: {.path {path}}.")
  lines <- readLines(path, warn = FALSE)
  lines <- trimws(vapply(strsplit(lines, "\t"), function(x) x[1], character(1)))
  lines <- lines[nzchar(lines) & !startsWith(lines, "#!")]
  if (!length(lines)) om_abort("Palette file {.path {path}} contained no colours.")
  lines
}

# Resolve the `palette` argument, which accepts a palette name, an explicit
# vector of colours, a named vector keyed by group, or a file path.
#' @noRd
resolve_palette <- function(palette, groups) {
  n <- length(groups)
  if (is.null(palette)) {
    return(stats::setNames(omakase_palette("omakase", n), groups))
  }
  if (is.character(palette) && length(palette) == 1L &&
      palette %in% names(OM_PALETTES)) {
    return(stats::setNames(omakase_palette(palette, n), groups))
  }
  if (is.character(palette) && length(palette) == 1L && file.exists(palette)) {
    palette <- read_palette(palette)
  }
  if (!is.null(names(palette))) {
    missing <- setdiff(groups, names(palette))
    if (length(missing)) {
      om_abort("Palette has no colour for group{?s} {.val {missing}}.")
    }
    return(palette[groups])
  }
  cols <- if (length(palette) < n) {
    grDevices::colorRampPalette(palette)(n)
  } else {
    palette[seq_len(n)]
  }
  stats::setNames(cols, groups)
}

# Roles are a small closed vocabulary, so an unknown one gets a neutral grey
# rather than an error - a user may legitimately mark a third transcript.
#' @noRd
resolve_role_fill <- function(roles, override = NULL) {
  roles <- unique(roles[!is.na(roles)])
  if (!length(roles)) return(character(0))
  out <- OM_ROLE_FILL[roles]
  names(out) <- roles
  out[is.na(out)] <- "#9E9E9E"
  if (!is.null(override)) {
    for (nm in intersect(names(override), roles)) out[[nm]] <- override[[nm]]
  }
  out
}

#' Discrete colour and fill scales using the omakase palettes
#'
#' @param palette Palette name, see [omakase_palettes()].
#' @param ... Passed to [ggplot2::discrete_scale()].
#' @return A ggplot2 scale.
#' @examples
#' library(ggplot2)
#' ggplot(mtcars, aes(wt, mpg, colour = factor(cyl))) +
#'   geom_point() +
#'   scale_colour_omakase()
#' @export
scale_colour_omakase <- function(palette = "omakase", ...) {
  ggplot2::discrete_scale(
    "colour", palette = function(n) omakase_palette(palette, n), ...
  )
}

#' @rdname scale_colour_omakase
#' @export
scale_color_omakase <- scale_colour_omakase

#' @rdname scale_colour_omakase
#' @export
scale_fill_omakase <- function(palette = "omakase", ...) {
  ggplot2::discrete_scale(
    "fill", palette = function(n) omakase_palette(palette, n), ...
  )
}
