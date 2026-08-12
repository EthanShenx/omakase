# ---------------------------------------------------------------------------
# The annotation panel that sits under a stack of coverage tracks: a coordinate
# bar with the window's end points, one row per transcript model, optional
# start-site markers, and an optional row of annotation features (repeats,
# motifs, peaks) drawn at their true coordinates rather than snapped to the
# transcript they belong to.
# ---------------------------------------------------------------------------

#' Draw the annotation panel for one locus
#'
#' Returns the bottom panel of a sashimi figure as a standalone `ggplot`: the
#' coordinate bar, the transcript models, and any features. Called by
#' [plot_sashimi()], but exported so the panel can be used on its own or
#' restyled before being assembled.
#'
#' @param x A `sashimi_data` object.
#' @param locus A `locus_id` or gene name. Defaults to the first locus.
#' @param xlim Plot-space x limits. Defaults to the locus window, widened by
#'   `psi_pad`.
#' @param map An [intron_map()] to draw in compressed coordinates, or `NULL`.
#' @param base_size Base font size in points.
#' @param base_family Font family for every piece of text in the figure. The
#'   empty string uses the graphics device's default. Text geoms do not inherit
#'   a theme's family, so this is threaded to each of them explicitly.
#' @param hairline Line width for the coordinate bar and intron lines.
#' @param role_fill Named vector of fills keyed by model `role`.
#' @param feature_color Fill for the feature boxes.
#' @param chrom_style How the contig name is printed above the coordinate bar:
#'   `"keep"` (the default) prints it exactly as the data has it, `"ucsc"`
#'   ensures a `chr` prefix, `"ensembl"` strips one.
#' @param show_coord_bar Draw the coordinate bar and window end labels.
#' @param coord_ticks Intermediate tick marks on the coordinate bar, labelled
#'   with the genomic coordinate each plot position maps back to. `TRUE` (the
#'   default) draws them only when introns are compressed, where an axis with
#'   just two labelled ends would imply an even scale that is not there. Pass a
#'   number for that many ticks, or `FALSE` for none.
#' @param show_tx_label Print each transcript's identifier at the left of its
#'   row.
#' @param show_features Draw the `features` slot.
#' @param show_apex Mark `main_apex`/`alt_apex` from the `loci` slot with a
#'   downward triangle.
#' @param show_gene_label Print the gene name, in italic, under the models.
#' @param collapse Draw one merged row per `role` instead of one row per
#'   transcript. Useful when an annotation has thirty isoforms and only the
#'   reference/alternative distinction matters.
#' @param arrow_bins Number of strand arrowheads to place along each intron
#'   line. `0` disables them.
#' @param theme A ggplot2 theme, or `NULL` for [theme_omakase()].
#'
#' @return A `ggplot` object.
#'
#' @examples
#' sd <- sashimi_data(
#'   loci = data.frame(locus_id = "a", gene_name = "DEMO", chrom = "chr1",
#'                     strand = "+", win_lo = 1000, win_hi = 5000),
#'   models = data.frame(locus_id = "a", tx_id = "tx1", role = "main",
#'                       start = c(1200, 3000), end = c(1600, 4200))
#' )
#' sashimi_annotation(sd)
#'
#' @export
sashimi_annotation <- function(x, locus = NULL, xlim = NULL, map = NULL,
                               base_size = 9, base_family = "",
                               hairline = 0.3,
                               role_fill = NULL,
                               feature_color = OM_FEATURE_FILL,
                               chrom_style = c("keep", "ucsc", "ensembl"),
                               show_coord_bar = TRUE, coord_ticks = TRUE,
                               show_tx_label = FALSE,
                               show_features = TRUE, show_apex = TRUE,
                               show_gene_label = TRUE, collapse = FALSE,
                               arrow_bins = 0, theme = NULL) {
  chrom_style <- match.arg(chrom_style)
  validate_sashimi_data(x)
  lid <- resolve_locus(x, locus)
  gi <- x$loci[x$loci$locus_id == lid, , drop = FALSE][1, ]
  md <- x$models[x$models$locus_id == lid, , drop = FALSE]
  ft <- x$features[x$features$locus_id == lid, , drop = FALSE]

  txt <- pt_to_mm(base_size)
  cx <- function(v) compress_coords(v, map)
  xr <- cx(c(gi$win_lo, gi$win_hi))
  xlim <- xlim %||% xr

  rows <- model_rows(md, collapse = collapse)
  n_row <- max(1L, nrow(rows$index))
  # The coordinate bar sits one row above the topmost transcript; the gene name
  # goes below the bottom one. The floor drops when features are drawn, because
  # a feature box and its label hang below their transcript's row.
  y_bar <- n_row + 1.15
  y_top <- n_row + 1.70
  # A feature box and its label hang below their transcript's row, but the
  # lowest row is at y = 1 and its label lands at 0.38, so the gene name at
  # 0.05 clears both. One pair of constants therefore serves either case.
  y_gene <- 0.05
  y_bot <- -0.05

  p <- ggplot2::ggplot()

  if (show_coord_bar) {
    # With introns compressed the axis is no longer linear, so the two end
    # labels alone would imply an even scale that is not there. Intermediate
    # ticks are placed at evenly spaced *plot* positions and labelled with the
    # genomic coordinates they map back to, which is exactly what the inverse
    # of the compression map is for.
    n_ticks <- if (isTRUE(coord_ticks)) {
      if (is.null(map) || identical(map$method, "none")) 0L else 3L
    } else {
      as.integer(coord_ticks %||% 0L)
    }
    if (n_ticks > 0) {
      at <- seq(xr[1], xr[2], length.out = n_ticks + 2)
      at <- at[-c(1, n_ticks + 2)]
      tk <- data.frame(x = at, lab = format_coord(expand_coords(at, map)))
      p <- p +
        ggplot2::geom_segment(
          data = tk,
          ggplot2::aes(x = .data$x, xend = .data$x, y = y_bar,
                       yend = y_bar + 0.09),
          linewidth = hairline, colour = "black"
        ) +
        ggplot2::geom_text(
          data = tk,
          ggplot2::aes(x = .data$x, y = y_bar - 0.30, label = .data$lab),
          size = txt * 0.9, colour = "black", family = base_family
        )
    }
    p <- p +
      ggplot2::annotate("segment", x = xr[1], xend = xr[2], y = y_bar,
                        yend = y_bar, linewidth = hairline, colour = "black") +
      ggplot2::annotate("text", x = mean(xr), y = y_bar + 0.35,
                        label = sprintf("%s : %s",
                                        harmonise_chrom(gi$chrom, chrom_style),
                                        gi$strand),
                        size = txt, colour = "black", family = base_family) +
      ggplot2::annotate("text", x = xr[1], y = y_bar - 0.30,
                        label = format_coord(gi$win_lo), hjust = 0,
                        size = txt, colour = "black", family = base_family) +
      ggplot2::annotate("text", x = xr[2], y = y_bar - 0.30,
                        label = format_coord(gi$win_hi), hjust = 1,
                        size = txt, colour = "black", family = base_family)
  }

  if (nrow(rows$exons)) {
    ex <- rows$exons
    ex$xmin <- cx(ex$start)
    ex$xmax <- cx(ex$end)
    sp <- rows$spans
    sp$x0 <- cx(sp$x0)
    sp$x1 <- cx(sp$x1)

    fills <- resolve_role_fill(ex$role, role_fill)
    # A CDS is drawn taller than a UTR, the way genome browsers do it, so the
    # coding part of a transcript reads at a glance.
    ex$half <- ifelse(ex$feature %in% c("UTR", "utr", "five_prime_UTR",
                                        "three_prime_UTR"), 0.09, 0.16)

    p <- p +
      ggplot2::geom_segment(
        data = sp, ggplot2::aes(x = .data$x0, xend = .data$x1,
                                y = .data$y, yend = .data$y),
        linewidth = hairline, colour = OM_MODEL_LINE
      ) +
      ggplot2::geom_rect(
        data = ex,
        ggplot2::aes(xmin = .data$xmin, xmax = .data$xmax,
                     ymin = .data$y - .data$half, ymax = .data$y + .data$half,
                     fill = .data$role),
        colour = NA
      ) +
      ggplot2::scale_fill_manual(values = fills, na.value = "#9E9E9E")

    if (arrow_bins > 0) {
      ar <- strand_arrows(sp, gi$strand, arrow_bins)
      if (nrow(ar)) {
        p <- p + ggplot2::geom_segment(
          data = ar,
          ggplot2::aes(x = .data$x, xend = .data$xend, y = .data$y,
                       yend = .data$y),
          linewidth = hairline, colour = OM_MODEL_LINE,
          arrow = grid::arrow(length = grid::unit(0.035, "npc"), type = "open")
        )
      }
    }
    if (show_tx_label) {
      p <- p + ggplot2::geom_text(
        data = rows$index,
        ggplot2::aes(x = xlim[1], y = .data$y + 0.30, label = .data$label),
        hjust = 0, size = txt, colour = "black", family = base_family
      )
    }
  }

  if (show_apex) {
    ap <- apex_points(gi, rows$index)
    if (nrow(ap)) {
      ap$x <- cx(ap$x)
      p <- p + ggplot2::geom_point(
        data = ap, ggplot2::aes(x = .data$x, y = .data$y + 0.30),
        shape = 25, size = 1.3, fill = "black", colour = "black"
      )
    }
  }

  if (show_features && nrow(ft)) {
    ft$y <- rows$y_for(ft$role)
    ft$xmin <- cx(ft$start)
    ft$xmax <- cx(ft$end)
    fcol <- ifelse(is.na(ft$color), feature_color, ft$color)
    p <- p +
      ggplot2::geom_rect(
        data = ft,
        ggplot2::aes(xmin = .data$xmin, xmax = .data$xmax,
                     ymin = .data$y - 0.46, ymax = .data$y - 0.28),
        fill = fcol, colour = NA
      ) +
      ggplot2::geom_text(
        data = ft,
        ggplot2::aes(x = (.data$xmin + .data$xmax) / 2, y = .data$y - 0.62,
                     label = .data$name),
        size = txt, colour = "black", family = base_family
      )
  }

  if (show_gene_label && !is.na(gi$gene_name)) {
    p <- p + ggplot2::annotate("text", x = mean(xr), y = y_gene,
                               label = gi$gene_name, size = txt,
                               colour = "black", fontface = "italic",
                               family = base_family)
  }

  p +
    ggplot2::coord_cartesian(xlim = xlim, ylim = c(y_bot, y_top),
                             expand = FALSE) +
    resolve_theme(theme, base_size, base_family)
}

# Lay transcripts out on rows. Reference roles go on top, and rows are numbered
# from the bottom so that adding a transcript does not move the others.
#' @noRd
model_rows <- function(md, collapse = FALSE) {
  empty <- list(
    exons = data.frame(), spans = data.frame(),
    index = data.frame(y = numeric(0), label = character(0),
                       role = character(0)),
    y_for = function(role) rep(1, length(role))
  )
  if (is.null(md) || !nrow(md)) return(empty)

  md$role <- ifelse(is.na(md$role), "main", md$role)
  key <- if (collapse) md$role else md$tx_id

  role_rank <- function(r) {
    pref <- c("main", "reference", "inclusion")
    ifelse(r %in% pref, 0L, 1L)
  }
  idx <- unique(data.frame(key = key, role = md$role,
                           stringsAsFactors = FALSE))
  idx <- idx[order(role_rank(idx$role), idx$role, idx$key), , drop = FALSE]
  n <- nrow(idx)
  idx$y <- n:1
  idx$label <- idx$key

  md$.key <- key
  md$y <- idx$y[match(md$.key, idx$key)]

  spans <- do.call(rbind, lapply(split(md, md$.key), function(d) {
    data.frame(y = d$y[1], x0 = min(d$start), x1 = max(d$end),
               role = d$role[1], stringsAsFactors = FALSE)
  }))
  rownames(spans) <- NULL

  list(
    exons = md,
    spans = spans,
    index = idx[, c("y", "label", "role", "key")],
    # Features name the role they belong to; anything unrecognised lands on the
    # bottom row rather than being dropped.
    y_for = function(role) {
      y <- idx$y[match(role, idx$role)]
      ifelse(is.na(y), min(idx$y), y)
    }
  )
}

# Start-site triangles. Uses the loci slot's apex columns when present, and
# otherwise falls back to the transcription start of each drawn model.
#' @noRd
apex_points <- function(gi, index) {
  pts <- list()
  if (!is.na(col(gi, "main_apex") %||% NA)) {
    y <- index$y[match("main", index$role)]
    pts[[length(pts) + 1]] <- data.frame(
      x = gi$main_apex, y = if (is.na(y)) max(index$y, 1) else y
    )
  }
  if (!is.na(col(gi, "alt_apex") %||% NA)) {
    y <- index$y[!index$role %in% c("main", "reference", "inclusion")][1]
    pts[[length(pts) + 1]] <- data.frame(
      x = gi$alt_apex, y = if (is.na(y)) min(index$y, 1) else y
    )
  }
  out <- rbind_all(pts)
  out %||% data.frame(x = numeric(0), y = numeric(0))
}

# Evenly spaced arrowheads showing the direction of transcription, drawn as
# very short segments so the head is all that shows.
#' @noRd
strand_arrows <- function(spans, strand, bins) {
  if (!nrow(spans) || bins <= 0 || identical(strand, "*")) {
    return(data.frame(x = numeric(0), xend = numeric(0), y = numeric(0)))
  }
  out <- lapply(seq_len(nrow(spans)), function(i) {
    x0 <- spans$x0[i]; x1 <- spans$x1[i]
    if (!is.finite(x0) || !is.finite(x1) || x1 <= x0) return(NULL)
    at <- seq(x0, x1, length.out = bins + 2)[-c(1, bins + 2)]
    step <- (x1 - x0) / (bins + 1) * 0.02
    if (identical(strand, "-")) step <- -step
    data.frame(x = at, xend = at + step, y = spans$y[i])
  })
  rbind_all(out) %||% data.frame(x = numeric(0), xend = numeric(0),
                                 y = numeric(0))
}

#' @noRd
resolve_locus <- function(x, locus) {
  if (!nrow(x$loci)) om_abort("This {.cls sashimi_data} object has no loci.")
  if (is.null(locus)) return(x$loci$locus_id[1])
  hit <- x$loci$locus_id == locus | x$loci$gene_name == locus
  if (!any(hit)) {
    om_abort(c("No locus matches {.val {locus}}.",
               "i" = "Available: {.val {utils::head(x$loci$locus_id, 10)}}."))
  }
  x$loci$locus_id[which(hit)[1]]
}
