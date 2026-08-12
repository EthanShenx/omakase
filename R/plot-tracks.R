# ---------------------------------------------------------------------------
# Genome-track figures.
#
# The other way to look at a locus: stacked horizontal tracks, each one a
# transcript model, a coverage profile or a row of features, with titles down
# the right-hand side and a shared coordinate axis. This is the pyGenomeTracks
# / IGV idiom, and it answers a different question from a sashimi plot - not
# "how is this spliced" but "what is here, and where".
#
# Tracks are ordinary values: build a list of them, reorder it, drop one, or
# write your own. plot_tracks() stacks whatever it is given.
# ---------------------------------------------------------------------------

#' Build a transcript model track
#'
#' Draws transcript models the way a genome browser does: exons as blocks,
#' introns as a line, and the direction of transcription marked by chevrons
#' along that line. With `style = "UCSC"` the coding part of each exon is drawn
#' taller than the untranslated part.
#'
#' @param models A data frame of exons with `start`, `end` and `tx_id`, such as
#'   [read_bed12()] or [read_annotation()] returns; a `sashimi_data` object,
#'   whose `models` slot is used; or a path to a BED/GTF/GFF3 file.
#' @param title Track title, printed to the right of the panel.
#' @param color Fill for the blocks. `"bed_rgb"` uses each record's own colour
#'   when the source carried one; a named vector keyed by `role` colours by
#'   role.
#' @param style `"UCSC"` draws CDS taller than UTR; `"flat"` draws every block
#'   the same height.
#' @param labels Print each transcript's name beside its row.
#' @param chevrons Number of direction marks per intron; `0` for none.
#' @param height Relative height of this track.
#' @param row_by Which column separates rows: `"tx_id"` (one row per
#'   transcript, the default) or `"role"` (one row per role).
#'
#' @return A `omakase_track` object, for [plot_tracks()].
#'
#' @examples
#' m <- data.frame(tx_id = "tx1", start = c(1000, 3000), end = c(1500, 4000),
#'                 strand = "+", feature = "CDS")
#' track_models(m, title = "models")
#'
#' @export
track_models <- function(models, title = NULL, color = "#1F78B4",
                         style = c("UCSC", "flat"), labels = FALSE,
                         chevrons = 18, height = 1, row_by = c("tx_id", "role")) {
  style <- match.arg(style)
  row_by <- match.arg(row_by)
  d <- resolve_model_input(models)
  new_track("models", data = d, title = title, color = color, style = style,
            labels = labels, chevrons = chevrons, height = height,
            row_by = row_by)
}

#' Build a coverage track
#'
#' A filled profile with its value range printed as a bracket on the left, the
#' way pyGenomeTracks labels a bedGraph.
#'
#' @param coverage A data frame with `pos` and `value`, a `sashimi_data`
#'   object, or a path to a bedGraph/bigWig-style file.
#' @param title Track title.
#' @param color Fill colour.
#' @param group When `coverage` is a `sashimi_data`, which group to draw.
#' @param ymax Upper limit; `NULL` uses the track's own maximum.
#' @param show_range Print the `0`-to-`ymax` bracket at the left.
#' @param log_y Draw on a `log10(1 + value)` axis.
#' @param height Relative height of this track.
#' @return An `omakase_track` object.
#' @examples
#' cv <- data.frame(pos = 1:100, value = abs(sin(seq(0, pi, length.out = 100))))
#' track_coverage(cv, title = "RPM")
#' @export
track_coverage <- function(coverage, title = NULL, color = "#435469",
                           group = NULL, ymax = NULL, show_range = TRUE,
                           log_y = FALSE, height = 0.7) {
  d <- resolve_coverage_input(coverage, group)
  new_track("coverage", data = d, title = title, color = color, ymax = ymax,
            show_range = show_range, log_y = log_y, height = height)
}

#' Build a feature track
#'
#' A row of blocks - peaks, repeats, motifs, start sites - optionally labelled.
#' With `shape = "marker"` each feature is drawn as a small downward pointer
#' rather than a block, which is what a single-base position wants.
#'
#' @param features A data frame with `start`, `end` and optionally `name`, a
#'   `sashimi_data` object, or a path to a BED file.
#' @param title Track title.
#' @param color Fill; `"bed_rgb"` uses each record's own colour.
#' @param labels Print feature names.
#' @param shape `"block"` or `"marker"`.
#' @param collapse Draw every feature on one row rather than spreading
#'   overlapping ones across rows.
#' @param height Relative height of this track.
#' @return An `omakase_track` object.
#' @examples
#' f <- data.frame(start = c(100, 500), end = c(200, 560),
#'                 name = c("peak1", "peak2"))
#' track_features(f, title = "peaks")
#' @export
track_features <- function(features, title = NULL, color = OM_FEATURE_FILL,
                           labels = TRUE, shape = c("block", "marker"),
                           collapse = TRUE, height = 0.4) {
  shape <- match.arg(shape)
  d <- resolve_feature_input(features)
  new_track("features", data = d, title = title, color = color,
            labels = labels, shape = shape, collapse = collapse,
            height = height)
}

#' Build a coordinate axis track
#'
#' @param unit `"auto"`, `"bp"`, `"Kb"` or `"Mb"`.
#' @param n Target number of ticks.
#' @param show_chrom Print the contig name under the axis.
#' @param height Relative height of this track.
#' @return An `omakase_track` object.
#' @examples
#' track_axis()
#' @export
track_axis <- function(unit = c("auto", "bp", "Kb", "Mb"), n = 6,
                       show_chrom = TRUE, height = 0.42) {
  unit <- match.arg(unit)
  new_track("axis", data = NULL, title = NULL, unit = unit, n = n,
            show_chrom = show_chrom, height = height)
}

#' Build a spacer track
#'
#' @param height Relative height of the gap.
#' @return An `omakase_track` object.
#' @examples
#' track_spacer(0.1)
#' @export
track_spacer <- function(height = 0.15) {
  new_track("spacer", data = NULL, title = NULL, height = height)
}

#' @noRd
new_track <- function(type, ...) {
  structure(c(list(type = type), list(...)), class = "omakase_track")
}

#' @export
print.omakase_track <- function(x, ...) {
  n <- if (is.null(x$data)) 0L else nrow(x$data)
  cli::cli_text("<omakase track> {.field {x$type}}{if (is.null(x$title)) '' else paste0(' - ', x$title)}\\
                {if (n) paste0(' (', n, ' row', if (n > 1) 's' else '', ')') else ''}")
  invisible(x)
}

#' Draw a genome-track figure
#'
#' Stacks tracks over a shared coordinate range, with titles down the right and
#' an axis at top or bottom. This is the browser-style view of a locus, as
#' opposed to the splicing-focused view [plot_sashimi()] gives.
#'
#' @details
#' Called on a `sashimi_data` object with no `tracks` argument, it builds a
#' sensible default stack: one model track, one coverage track per group, a
#' feature track if the object has features, and an axis. Pass `tracks` to take
#' full control.
#'
#' @param x A `sashimi_data` object, or `NULL` when `tracks` is given.
#' @param tracks A list of tracks from [track_models()], [track_coverage()],
#'   [track_features()], [track_axis()] and [track_spacer()]. Built from `x`
#'   when `NULL`.
#' @param region The coordinate range, as a region string or a
#'   [parse_region()] object. Taken from `x` when `NULL`.
#' @param locus When `x` holds several loci, which one to draw.
#' @param axis `"bottom"`, `"top"` or `"none"`; ignored when `tracks` already
#'   contains an axis.
#' @param palette Palette for the per-group coverage tracks.
#' @param title Figure title.
#' @param title_width Fraction of the figure width reserved for track titles.
#' @param left_pad Fraction of the width reserved at the left for the coverage
#'   range labels.
#' @param base_size Base font size in points.
#' @param base_family Font family for every piece of text in the figure. The
#'   empty string uses the graphics device's default. Text geoms do not inherit
#'   a theme's family, so this is threaded to each of them explicitly.
#' @param hairline Line width for rules.
#' @param shrink Compress introns; `TRUE`, or an [intron_map()].
#' @param shrink_method,shrink_gamma Compression rule and its exponent.
#' @param chrom_style How the contig name is printed: `"keep"`, `"ucsc"` or
#'   `"ensembl"`.
#' @param theme A ggplot2 theme, or `NULL` for [theme_omakase()].
#'
#' @return A `patchwork` object.
#'
#' @examples
#' sd <- sashimi_data(
#'   loci = data.frame(locus_id = "a", gene_name = "DEMO", chrom = "chr1",
#'                     strand = "+", win_lo = 1000, win_hi = 5000),
#'   tracks = data.frame(locus_id = "a", group = "sample1",
#'                       pos = seq(1000, 5000, 40),
#'                       value = abs(sin(seq(0, pi, length.out = 101)))),
#'   models = data.frame(locus_id = "a", tx_id = "tx1", role = "main",
#'                       feature = "CDS",
#'                       start = c(1200, 3000), end = c(1600, 4200))
#' )
#' plot_tracks(sd)
#'
#' @export
plot_tracks <- function(x = NULL, tracks = NULL, region = NULL, locus = NULL,
                        axis = c("bottom", "top", "none"),
                        palette = NULL, title = NULL, title_width = 0.24,
                        left_pad = 0.05,
                        base_size = 9, base_family = "", hairline = 0.3,
                        shrink = FALSE, shrink_method = "power",
                        shrink_gamma = 0.7, chrom_style = "keep",
                        theme = NULL) {
  axis <- match.arg(axis)

  gi <- NULL
  if (!is.null(x)) {
    x <- as_sashimi_data(x)
    validate_sashimi_data(x)
    lid <- resolve_locus(x, locus)
    gi <- x$loci[x$loci$locus_id == lid, , drop = FALSE][1, ]
    x <- x[lid]
  }

  r <- if (!is.null(region)) {
    parse_region(region)
  } else if (!is.null(gi)) {
    new_region(gi$chrom, gi$win_lo, gi$win_hi, gi$strand)
  } else {
    om_abort("Give a {.arg region}, or an {.cls sashimi_data} object to take one from.")
  }

  map <- NULL
  if (isTRUE(shrink) && !is.null(x)) {
    map <- intron_map(introns_from_models(x$models), r$start, r$end,
                      method = shrink_method, gamma = shrink_gamma)
  } else if (inherits(shrink, "omakase_intron_map")) {
    map <- shrink
  }

  if (is.null(tracks)) {
    if (is.null(x)) om_abort("Give {.arg tracks}, or an {.cls sashimi_data} object to build them from.")
    tracks <- default_tracks(x, palette, axis)
  } else {
    if (!is.null(axis) && !identical(axis, "none") &&
        !any(vapply(tracks, function(t) identical(t$type, "axis"), logical(1)))) {
      ax <- track_axis()
      tracks <- if (identical(axis, "top")) c(list(ax), tracks) else c(tracks, list(ax))
    }
  }

  xr <- compress_coords(c(r$start, r$end), map)
  panels <- lapply(tracks, function(tr) {
    draw_track(tr, r = r, xr = xr, map = map, base_size = base_size,
               base_family = base_family, hairline = hairline,
               chrom_style = chrom_style, theme = theme,
               title_width = title_width, left_pad = left_pad)
  })
  heights <- vapply(tracks, function(t) t$height %||% 1, numeric(1))

  p <- patchwork::wrap_plots(panels, ncol = 1, heights = heights)
  if (!is.null(title)) {
    p <- p + patchwork::plot_annotation(
      title = title,
      theme = ggplot2::theme(
        plot.title = ggplot2::element_text(size = base_size, colour = "black",
                                           face = "plain", hjust = 0.5,
                                           family = base_family)
      )
    )
  }
  attr(p, "omakase_tracks") <- length(tracks)
  p
}

# The default stack for a sashimi_data: models, then one coverage row per
# group, then features, then the axis.
#' @noRd
default_tracks <- function(x, palette, axis) {
  out <- list()
  gi <- x$loci[1, ]
  if (nrow(x$models)) {
    out[[length(out) + 1]] <- track_models(
      x, title = paste(gi$gene_name, "models"),
      color = if (any(!is.na(x$models$role))) OM_ROLE_FILL else "#1F78B4",
      row_by = if (any(!is.na(x$models$role))) "role" else "tx_id"
    )
  }
  if (nrow(x$tracks)) {
    groups <- order_levels(unique(x$tracks$group), x$meta$group_order)
    cols <- resolve_palette(palette, groups)
    # One shared ceiling, so the rows can be read against each other - which is
    # the point of stacking them.
    ymax <- max(x$tracks$value, na.rm = TRUE)
    for (g in groups) {
      out[[length(out) + 1]] <- track_coverage(
        x, title = g, color = cols[[g]], group = g, ymax = ymax
      )
    }
  }
  if (nrow(x$features)) {
    out[[length(out) + 1]] <- track_features(x, title = "features")
  }
  if (!identical(axis, "none")) {
    ax <- track_axis()
    out <- if (identical(axis, "top")) c(list(ax), out) else c(out, list(ax))
  }
  out
}

# ---------------------------------------------------------------------------
# Input coercion. Each track builder accepts a data frame, a sashimi_data, or a
# path, so a caller never has to remember which.
# ---------------------------------------------------------------------------

#' @noRd
resolve_model_input <- function(m) {
  if (inherits(m, "sashimi_data")) {
    d <- m$models
    d$chrom <- m$loci$chrom[match(d$locus_id, m$loci$locus_id)]
    if (all(is.na(d$strand) | d$strand == "*")) {
      d$strand <- m$loci$strand[match(d$locus_id, m$loci$locus_id)]
    }
    return(d)
  }
  if (is.character(m) && length(m) == 1L) {
    d <- if (grepl("\\.bed(\\.gz)?$", m, ignore.case = TRUE)) {
      read_bed(m)
    } else {
      a <- read_annotation(m)
      names(a)[names(a) == "feature"] <- "feature"
      a
    }
    return(d)
  }
  d <- as_df(m)
  require_cols(d, c("start", "end"), "model table")
  if (!has_col(d, "tx_id")) d$tx_id <- "tx"
  if (!has_col(d, "feature")) d$feature <- "exon"
  if (!has_col(d, "strand")) d$strand <- "*"
  d
}

#' @noRd
resolve_coverage_input <- function(cv, group) {
  if (inherits(cv, "sashimi_data")) {
    d <- cv$tracks
    if (!is.null(group)) d <- d[d$group %in% group, , drop = FALSE]
    return(d)
  }
  if (is.character(cv) && length(cv) == 1L) {
    d <- utils::read.delim(cv, header = FALSE, sep = "\t",
                           stringsAsFactors = FALSE, comment.char = "#")
    names(d)[1:4] <- c("chrom", "start", "end", "value")
    # A bedGraph gives intervals; the midpoint is what a line wants.
    return(data.frame(pos = (as.numeric(d$start) + 1 + as.numeric(d$end)) / 2,
                      value = as.numeric(d$value),
                      stringsAsFactors = FALSE))
  }
  d <- as_df(cv)
  require_cols(d, c("pos", "value"), "coverage table")
  d
}

#' @noRd
resolve_feature_input <- function(f) {
  if (inherits(f, "sashimi_data")) return(f$features)
  if (is.character(f) && length(f) == 1L) {
    d <- read_bed(f, expand = FALSE)
    d$name <- d$tx_id
    return(d)
  }
  d <- as_df(f)
  require_cols(d, c("start", "end"), "feature table")
  if (!has_col(d, "name")) d$name <- NA_character_
  d
}

# ---------------------------------------------------------------------------
# Drawing
# ---------------------------------------------------------------------------

#' @noRd
draw_track <- function(tr, r, xr, map, base_size, base_family, hairline,
                       chrom_style, theme, title_width, left_pad = 0.05) {
  # Titles live in a reserved strip to the right of the data and the coverage
  # range labels in one to the left, so every panel shares a single x range and
  # the tracks stay in register.
  w <- diff(xr)
  xlim <- c(xr[1] - w * left_pad,
            xr[2] + w * title_width / (1 - title_width))
  p <- switch(tr$type,
    models   = draw_models_track(tr, xr, map, base_size, base_family, hairline),
    coverage = draw_coverage_track(tr, xr, map, base_size, base_family, hairline),
    features = draw_features_track(tr, xr, map, base_size, base_family, hairline),
    axis     = draw_axis_track(tr, r, xr, map, base_size, base_family,
                               hairline, chrom_style),
    spacer   = list(plot = ggplot2::ggplot(), ylim = c(0, 1))
  )
  g <- p$plot
  if (!is.null(tr$title) && !identical(tr$type, "spacer")) {
    g <- g + ggplot2::annotate(
      "text", x = xr[2] + diff(xr) * 0.02, y = mean(p$ylim),
      label = tr$title, hjust = 0, size = pt_to_mm(base_size), colour = "black",
      family = base_family
    )
  }
  g +
    ggplot2::coord_cartesian(xlim = xlim, ylim = p$ylim, expand = FALSE) +
    resolve_theme(theme, base_size, base_family, margin = c(0, 2, 0, 2))
}

#' @noRd
draw_models_track <- function(tr, xr, map, base_size, base_family, hairline) {
  d <- tr$data
  cx <- function(v) compress_coords(v, map)
  if (is.null(d) || !nrow(d)) {
    return(list(plot = ggplot2::ggplot(), ylim = c(0, 1)))
  }
  key <- if (identical(tr$row_by, "role") && has_col(d, "role")) {
    ifelse(is.na(d$role), "tx", as.character(d$role))
  } else {
    as.character(d$tx_id)
  }
  d$.key <- key
  lv <- unique(key)
  n <- length(lv)
  d$y <- n - match(d$.key, lv) + 1

  spans <- do.call(rbind, lapply(split(d, d$.key), function(z) {
    data.frame(y = z$y[1], x0 = min(z$start), x1 = max(z$end),
               strand = z$strand[1], key = z$.key[1],
               # The row may be keyed by role, but the label a reader wants is
               # still the transcript's own name.
               label = as.character(z$tx_id[1]), stringsAsFactors = FALSE)
  }))
  rownames(spans) <- NULL
  # A transcript usually runs past the window; clipping here rather than
  # dropping it keeps the model visible but stops the line ploughing through
  # the title strip.
  spans$x0c <- pmax(xr[1], cx(spans$x0))
  spans$x1c <- pmin(xr[2], cx(spans$x1))
  spans <- spans[spans$x1c > spans$x0c, , drop = FALSE]

  d$xmin <- pmax(xr[1], cx(d$start))
  d$xmax <- pmin(xr[2], cx(d$end))
  d <- d[d$xmax > d$xmin, , drop = FALSE]
  if (!nrow(d)) return(list(plot = ggplot2::ggplot(), ylim = c(0, 1)))
  # UCSC style: the coding part stands taller than the untranslated part, so
  # the reading frame is visible without a legend.
  d$half <- if (identical(tr$style, "UCSC")) {
    ifelse(toupper(d$feature) %in% c("CDS"), 0.34, 0.17)
  } else {
    0.30
  }
  d$fill <- track_fill(d, tr$color)

  p <- ggplot2::ggplot() +
    ggplot2::geom_segment(
      data = spans,
      ggplot2::aes(x = .data$x0c, xend = .data$x1c, y = .data$y, yend = .data$y),
      linewidth = hairline, colour = "black"
    )

  if (tr$chevrons > 0) {
    ch <- chevron_points(spans, tr$chevrons)
    if (nrow(ch)) {
      p <- p + ggplot2::geom_segment(
        data = ch,
        ggplot2::aes(x = .data$x, xend = .data$xend, y = .data$y,
                     yend = .data$y),
        linewidth = hairline * 0.8, colour = "black",
        # Absolute units: an npc-relative arrowhead grows with the panel and
        # swamps a short track.
        arrow = grid::arrow(length = grid::unit(1.1, "mm"), type = "open",
                            angle = 30)
      )
    }
  }

  p <- p + ggplot2::geom_rect(
    data = d,
    ggplot2::aes(xmin = .data$xmin, xmax = .data$xmax,
                 ymin = .data$y - .data$half, ymax = .data$y + .data$half,
                 fill = .data$fill),
    colour = NA
  ) + ggplot2::scale_fill_identity()

  if (isTRUE(tr$labels)) {
    p <- p + ggplot2::geom_text(
      data = spans,
      ggplot2::aes(x = .data$x1c, y = .data$y, label = .data$label),
      hjust = -0.06, size = pt_to_mm(base_size * 0.9), colour = "black",
      family = base_family
    )
  }
  list(plot = p, ylim = c(0.35, n + 0.65))
}

#' @noRd
draw_coverage_track <- function(tr, xr, map, base_size, base_family, hairline) {
  d <- tr$data
  cx <- function(v) compress_coords(v, map)
  if (is.null(d) || !nrow(d)) {
    return(list(plot = ggplot2::ggplot(), ylim = c(0, 1)))
  }
  d <- d[order(d$pos), , drop = FALSE]
  d$.x <- cx(d$pos)
  d$.y <- if (isTRUE(tr$log_y)) log10(1 + pmax(0, d$value)) else d$value

  ymax <- tr$ymax
  if (!is.null(ymax) && isTRUE(tr$log_y)) ymax <- log10(1 + ymax)
  ymax <- ymax %||% max(d$.y, na.rm = TRUE)
  if (!is.finite(ymax) || ymax <= 0) ymax <- 1

  p <- ggplot2::ggplot() +
    ggplot2::geom_area(data = d, ggplot2::aes(x = .data$.x, y = .data$.y),
                       fill = tr$color, colour = NA)

  if (isTRUE(tr$show_range)) {
    # The bracket pyGenomeTracks draws: a square hook at the left edge with the
    # ceiling and the floor written against it.
    tick <- diff(xr) * 0.012
    br <- data.frame(
      x = c(xr[1] + tick, xr[1], xr[1], xr[1] + tick),
      y = c(ymax, ymax, 0, 0)
    )
    p <- p +
      ggplot2::geom_path(data = br, ggplot2::aes(x = .data$x, y = .data$y),
                         linewidth = hairline, colour = "black") +
      ggplot2::annotate("text", x = xr[1] - tick * 0.6, y = ymax,
                        label = format_activity(tr$ymax %||% max(d$value, na.rm = TRUE)),
                        hjust = 1, vjust = 0.85, family = base_family,
                        size = pt_to_mm(base_size * 0.85), colour = "black") +
      ggplot2::annotate("text", x = xr[1] - tick * 0.6, y = 0, label = "0",
                        hjust = 1, vjust = 0.15, family = base_family,
                        size = pt_to_mm(base_size * 0.85), colour = "black")
  }
  list(plot = p, ylim = c(0, ymax * 1.05))
}

#' @noRd
draw_features_track <- function(tr, xr, map, base_size, base_family, hairline) {
  d <- tr$data
  cx <- function(v) compress_coords(v, map)
  if (is.null(d) || !nrow(d)) {
    return(list(plot = ggplot2::ggplot(), ylim = c(0, 1)))
  }
  d$xmin <- cx(d$start); d$xmax <- cx(d$end)
  d$y <- 1
  d$fill <- track_fill(d, tr$color)

  p <- ggplot2::ggplot()
  if (identical(tr$shape, "marker")) {
    # A single-base position has no width to draw, so it becomes a pointer.
    w <- diff(xr) * 0.004
    mk <- do.call(rbind, lapply(seq_len(nrow(d)), function(i) {
      cxm <- (d$xmin[i] + d$xmax[i]) / 2
      data.frame(x = c(cxm - w, cxm + w, cxm), y = c(1.35, 1.35, 0.75),
                 g = i, fill = d$fill[i], stringsAsFactors = FALSE)
    }))
    p <- p + ggplot2::geom_polygon(
      data = mk,
      ggplot2::aes(x = .data$x, y = .data$y, group = .data$g,
                   fill = .data$fill),
      colour = NA
    ) + ggplot2::scale_fill_identity()
  } else {
    p <- p + ggplot2::geom_rect(
      data = d,
      ggplot2::aes(xmin = .data$xmin, xmax = .data$xmax,
                   ymin = .data$y - 0.28, ymax = .data$y + 0.28,
                   fill = .data$fill),
      colour = NA
    ) + ggplot2::scale_fill_identity()
  }
  if (isTRUE(tr$labels) && has_col(d, "name") && any(!is.na(d$name))) {
    p <- p + ggplot2::geom_text(
      data = d[!is.na(d$name), , drop = FALSE],
      ggplot2::aes(x = (.data$xmin + .data$xmax) / 2, y = 0.42,
                   label = .data$name),
      size = pt_to_mm(base_size * 0.85), colour = "black", vjust = 1,
      family = base_family
    )
  }
  list(plot = p, ylim = c(0.05, 1.55))
}

#' @noRd
draw_axis_track <- function(tr, r, xr, map, base_size, base_family, hairline,
                            chrom_style) {
  # Ticks are chosen evenly in plot space and labelled with the genomic
  # coordinates they map back to, so a compressed axis still reads correctly.
  at <- seq(xr[1], xr[2], length.out = tr$n)
  gpos <- expand_coords(at, map)
  span <- diff(range(gpos))
  unit <- if (identical(tr$unit, "auto")) {
    if (span >= 2e6) "Mb" else if (span >= 2e3) "Kb" else "bp"
  } else {
    tr$unit
  }
  div <- switch(unit, bp = 1, Kb = 1e3, Mb = 1e6)
  # Enough decimals to tell adjacent ticks apart, and no more: three decimal
  # places on a kilobase axis is noise, not precision.
  step <- span / max(1, length(at) - 1) / div
  digits <- if (step >= 10) 0 else if (step >= 1) 1 else 2
  lab <- formatC(gpos / div, format = "f", digits = digits, big.mark = ",")
  lab[length(lab)] <- paste(lab[length(lab)], unit)
  # The end labels would otherwise be half outside the panel.
  hj <- c(0, rep(0.5, length(at) - 2), 1)

  p <- ggplot2::ggplot() +
    ggplot2::annotate("segment", x = xr[1], xend = xr[2], y = 1, yend = 1,
                      linewidth = hairline, colour = "black") +
    ggplot2::annotate("segment", x = at, xend = at, y = 1, yend = 1.22,
                      linewidth = hairline, colour = "black") +
    ggplot2::annotate("text", x = at, y = 0.80, label = lab, hjust = hj,
                      size = pt_to_mm(base_size * 0.9), colour = "black",
                      vjust = 1, family = base_family)
  if (isTRUE(tr$show_chrom)) {
    p <- p + ggplot2::annotate(
      "text", x = mean(xr), y = -0.30,
      label = harmonise_chrom(r$chrom, chrom_style),
      size = pt_to_mm(base_size * 0.9), colour = "black", vjust = 1,
      family = base_family
    )
  }
  list(plot = p, ylim = c(-0.95, 1.35))
}

# Resolve a track's `color`: a single colour, "bed_rgb" to use each record's
# own, or a named vector keyed by role.
#' @noRd
track_fill <- function(d, color) {
  if (identical(color, "bed_rgb")) {
    v <- if (has_col(d, "color")) d$color else NA_character_
    return(ifelse(is.na(v), "#1F78B4", v))
  }
  if (length(color) > 1L && !is.null(names(color))) {
    key <- if (has_col(d, "role")) as.character(d$role) else rep(NA_character_, nrow(d))
    v <- unname(color[key])
    return(ifelse(is.na(v), "#9E9E9E", v))
  }
  rep(color[1], nrow(d))
}

# Evenly spaced direction marks along each transcript's intron line.
#' @noRd
chevron_points <- function(spans, n) {
  out <- lapply(seq_len(nrow(spans)), function(i) {
    x0 <- spans$x0c[i]; x1 <- spans$x1c[i]
    if (!is.finite(x0) || !is.finite(x1) || x1 <= x0) return(NULL)
    at <- seq(x0, x1, length.out = n + 2)[-c(1, n + 2)]
    step <- (x1 - x0) * 0.004
    if (identical(spans$strand[i], "-")) step <- -step
    data.frame(x = at - step / 2, xend = at + step / 2, y = spans$y[i])
  })
  rbind_all(out) %||% data.frame(x = numeric(0), xend = numeric(0),
                                 y = numeric(0))
}

#' Save a genome-track figure
#'
#' @param plot A plot from [plot_tracks()].
#' @param file Output path; the extension selects the device.
#' @param width Width in inches.
#' @param height Height in inches, or `NULL` to size it from the track count.
#' @param dpi Resolution for raster devices.
#' @param ... Passed to [ggplot2::ggsave()].
#' @return `file`, invisibly.
#' @examples
#' sd <- sashimi_data(
#'   loci = data.frame(locus_id = "a", gene_name = "A", chrom = "chr1",
#'                     strand = "+", win_lo = 1, win_hi = 100),
#'   models = data.frame(locus_id = "a", tx_id = "t", start = 10, end = 90)
#' )
#' save_tracks(plot_tracks(sd), file.path(tempdir(), "tracks.pdf"))
#' @export
save_tracks <- function(plot, file, width = 8, height = NULL, dpi = 300, ...) {
  n <- attr(plot, "omakase_tracks") %||% 4
  height <- height %||% (0.42 * n + 0.6)
  ggplot2::ggsave(file, plot, width = width, height = height, units = "in",
                  dpi = dpi, ...)
  invisible(file)
}
