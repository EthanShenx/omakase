# ---------------------------------------------------------------------------
# The sashimi figure: a stack of coverage panels with junction arcs, over an
# annotation panel. Assembled with patchwork and returned as a plot object, so
# a user can restyle any panel afterwards instead of re-running the whole
# pipeline with different flags.
# ---------------------------------------------------------------------------

#' Draw a single coverage panel with its junction arcs
#'
#' One row of a sashimi figure. Exported so a panel can be built, inspected or
#' restyled on its own; [plot_sashimi()] calls this once per group and stacks
#' the results.
#'
#' @param x A `sashimi_data` object.
#' @param locus A `locus_id` or gene name. Defaults to the first locus.
#' @param group The group (sample, stage, condition) to draw. Defaults to the
#'   first. Several groups may be given, in which case their coverage is
#'   overlaid in this one panel, each in its own colour from `fill`.
#' @param xlim Plot-space x limits; defaults to the locus window widened by
#'   `psi_pad`.
#' @param ymax Upper limit of the coverage axis. `NULL` takes the panel's own
#'   maximum, which is what makes each panel use its full height; pass a shared
#'   value (as [plot_sashimi()] does when `fix_y_scale = TRUE`) to make panels
#'   quantitatively comparable.
#' @param map An [intron_map()] for compressed coordinates, or `NULL`.
#' @param fill,arc_color Colours for the coverage area and the arcs. `arc_color`
#'   defaults to `fill`.
#' @param alpha Opacity of the coverage area.
#' @param base_size Base font size in points.
#' @param base_family Font family for every piece of text in the figure. The
#'   empty string uses the graphics device's default. Text geoms do not inherit
#'   a theme's family, so this is threaded to each of them explicitly.
#' @param arc_shape One of [arc_shapes()].
#' @param arc_height_rule One of [arc_height_rules()]. `"auto"` (the default)
#'   staggers a handful of arcs at the `arc_height_frac` heights and switches to
#'   `"span"` once there are enough junctions that equal-height arcs would
#'   cross; `"constant"` always staggers; `"span"` scales height with junction
#'   width so nested junctions nest; the rest scale height with count.
#' @param arc_height_frac Arc apex heights as a fraction of `ymax`. The values
#'   are cycled across arcs - arcs sharing an end point first - so their labels
#'   cannot collide.
#' @param arc_width_rule,arc_width Junction line width: a rule from
#'   [arc_widths()] and its scale factor.
#' @param arc_side `"above"` or `"below"` the coverage.
#' @param arc_n Points per arc path.
#' @param show_arc_label Print the count on each arc.
#' @param arc_label_format A formatter: a function, or one of `"activity"`,
#'   `"count"`, `"coord"`, `"none"`.
#' @param arc_label_position Where the count sits relative to its arc.
#'   `"on"` (the default) puts it astride the apex, its opaque box
#'   interrupting the curve, which keeps the digits readable where arcs cross;
#'   `"above"` floats it clear of an unbroken curve; `"none"` is the same as
#'   `show_arc_label = FALSE`.
#' @param label_background Draw the arc label on an opaque box, so an arc
#'   passing underneath does not strike through the digits. Only applies when
#'   `arc_label_position = "on"`.
#' @param label_padding Padding inside the label box, in points.
#' @param label_offset Gap between the arc apex and a label placed `"above"`,
#'   as a fraction of the panel height.
#' @param background Fill for the coverage panel. `NA` (the default) leaves it
#'   transparent; a colour tints the whole panel, which separates stacked
#'   tracks without drawing rules between them.
#' @param background_alpha Opacity of `background`.
#' @param min_count Junctions with a count below this are not drawn.
#' @param group_label Text printed at the top left of the panel. `NULL` uses
#'   the group name; `NA` prints nothing.
#' @param psi_label Text printed in the right-hand gutter, or `NULL` to take it
#'   from the `psi` slot.
#' @param psi_pad Width of the right-hand gutter, as a fraction of the window.
#' @param log_y Draw the coverage on a `log10(1 + value)` axis.
#' @param overlay_junction_fun How a junction shared by several overlaid groups
#'   is combined into one arc: `"mean"` (the default), `"median"`, `"sum"` or
#'   `"max"`. Only used when `group` names more than one group.
#' @param theme A ggplot2 theme, or `NULL` for [theme_omakase()].
#'
#' @return A `ggplot` object.
#'
#' @examples
#' sd <- sashimi_data(
#'   loci = data.frame(locus_id = "a", gene_name = "A", chrom = "chr1",
#'                     strand = "+", win_lo = 1000, win_hi = 2000),
#'   tracks = data.frame(locus_id = "a", group = "g", pos = seq(1000, 2000, 10),
#'                       value = abs(sin(seq(0, pi, length.out = 101))) * 20),
#'   junctions = data.frame(locus_id = "a", group = "g", x0 = 1200, x1 = 1800,
#'                          count = 42)
#' )
#' sashimi_track(sd)
#'
#' @export
sashimi_track <- function(x, locus = NULL, group = NULL, xlim = NULL,
                          ymax = NULL, map = NULL,
                          fill = "#66C2A5", arc_color = NULL, alpha = 1,
                          base_size = 9, base_family = "",
                          arc_shape = "sine", arc_height_rule = "auto",
                          arc_height_frac = c(0.80, 1.20),
                          arc_width_rule = "constant", arc_width = 0.5,
                          arc_side = "above", arc_n = 121,
                          show_arc_label = TRUE,
                          arc_label_format = "activity",
                          arc_label_position = c("on", "above", "none"),
                          label_background = TRUE, label_padding = 0.6,
                          label_offset = 0.03,
                          background = NA, background_alpha = 1,
                          min_count = 0,
                          group_label = NULL, psi_label = NULL,
                          psi_pad = 0.30, log_y = FALSE,
                          overlay_junction_fun = "mean", theme = NULL) {
  validate_sashimi_data(x)
  lid <- resolve_locus(x, locus)
  gi <- x$loci[x$loci$locus_id == lid, , drop = FALSE][1, ]
  grp <- group %||% x$tracks$group[x$tracks$locus_id == lid][1] %||% "1"

  cv <- x$tracks[x$tracks$locus_id == lid & x$tracks$group %in% grp, , drop = FALSE]
  jn <- x$junctions[x$junctions$locus_id == lid & x$junctions$group %in% grp, ,
                    drop = FALSE]
  if (nrow(jn) && min_count > 0) jn <- jn[jn$count >= min_count, , drop = FALSE]
  if (length(grp) > 1L && nrow(jn)) {
    # An overlay panel stands for several groups at once, so a junction they
    # share must become one arc. Drawing it once per group stacks identical
    # curves and piles their labels on top of each other.
    jn <- stats::aggregate(
      jn$count,
      by = list(x0 = jn$x0, x1 = jn$x1,
                role = ifelse(is.na(jn$role), "", jn$role)),
      FUN = function(v) do.call(overlay_junction_fun, list(v, na.rm = TRUE))
    )
    names(jn)[ncol(jn)] <- "count"
    jn$role[jn$role == ""] <- NA_character_
    jn$group <- grp[1]
    jn$label <- NA_character_
  }

  arc_color <- arc_color %||% fill[1]
  txt <- pt_to_mm(base_size)
  fmt <- resolve_formatter(arc_label_format)
  cx <- function(v) compress_coords(v, map)

  tf <- if (log_y) function(v) log10(1 + pmax(0, v)) else identity
  if (nrow(cv)) cv$.y <- tf(cv$value)

  # A panel with no signal still needs a height, or the arcs and labels have
  # nothing to be positioned against.
  ymax <- ymax %||% if (nrow(cv)) max(cv$.y, na.rm = TRUE) else 1
  if (!is.finite(ymax) || ymax <= 0) ymax <- 1

  xr <- cx(c(gi$win_lo, gi$win_hi))
  xlim <- xlim %||% c(xr[1], xr[2] + diff(xr) * psi_pad)

  p <- ggplot2::ggplot()

  if (!is.na(background)) {
    # An annotate() rectangle rather than a theme panel fill, so it spans
    # exactly the coordinate range and so patchwork cannot inset it.
    p <- p + ggplot2::annotate(
      "rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf,
      fill = background, alpha = background_alpha, colour = NA
    )
  }

  if (nrow(cv)) {
    cv <- cv[order(cv$pos), , drop = FALSE]
    cv$.x <- cx(cv$pos)
    if (length(grp) > 1L) {
      # Overlay: one filled area per group in the same panel. Drawing the
      # largest first means a smaller profile is never hidden behind a bigger
      # one, which is the failure mode that makes overlays unreadable.
      cols <- if (length(fill) == length(grp)) {
        stats::setNames(fill, grp)
      } else {
        resolve_palette(NULL, grp)
      }
      ord <- names(sort(tapply(cv$.y, cv$group, max, na.rm = TRUE),
                        decreasing = TRUE))
      for (g in ord) {
        one <- cv[cv$group == g, , drop = FALSE]
        p <- p + ggplot2::geom_area(
          data = one, ggplot2::aes(x = .data$.x, y = .data$.y),
          fill = cols[[g]], colour = NA, alpha = alpha
        )
      }
    } else {
      p <- p + ggplot2::geom_area(
        data = cv, ggplot2::aes(x = .data$.x, y = .data$.y),
        fill = fill[1], colour = NA, alpha = alpha
      )
    }
  }

  if (nrow(jn)) {
    jn$x0 <- cx(jn$x0)
    jn$x1 <- cx(jn$x1)
    arcs <- build_arcs(jn, ymax = ymax, shape = arc_shape,
                       height_rule = arc_height_rule,
                       height_frac = arc_height_frac, n = arc_n,
                       side = arc_side)
    if (!is.null(arcs$paths) && nrow(arcs$paths)) {
      # Widths are computed on the sorted table build_arcs returns, so .id
      # indexes them directly.
      lw <- arc_widths(arcs$junctions$count, rule = arc_width_rule,
                       w0 = arc_width)
      arcs$paths$.lw <- lw[arcs$paths$.id]
      p <- p + ggplot2::geom_path(
        data = arcs$paths,
        ggplot2::aes(x = .data$x, y = .data$y, group = .data$.id,
                     linewidth = .data$.lw),
        colour = arc_color, lineend = "round"
      ) + ggplot2::scale_linewidth_identity()
    }
    pos <- rlang::arg_match(arc_label_position, c("on", "above", "none"))
    if (show_arc_label && !identical(pos, "none") &&
        !is.null(arcs$labels) && nrow(arcs$labels)) {
      lab <- arcs$labels
      lab$txt <- ifelse(is.na(lab$label), fmt(lab$count), lab$label)

      if (identical(pos, "on")) {
        # The label sits astride the apex and its opaque box interrupts the
        # curve, which is what keeps the digits readable when several arcs
        # cross. This is the default and the reference style.
        p <- p + if (label_background) {
          ggplot2::geom_label(
            data = lab,
            ggplot2::aes(x = .data$x, y = .data$y, label = .data$txt),
            size = txt, colour = "black", family = base_family,
            fill = "white", linewidth = 0,
            label.padding = grid::unit(label_padding, "pt")
          )
        } else {
          ggplot2::geom_text(
            data = lab,
            ggplot2::aes(x = .data$x, y = .data$y, label = .data$txt),
            size = txt, colour = "black", family = base_family
          )
        }
      } else {
        # Clear of the curve: the arc is drawn unbroken and the number floats
        # just above its apex. Preferable when arc shape is what matters, or
        # when a boxed label would punch a hole in a dense figure.
        lab$y <- lab$y + sign_of(arc_side) * ymax * label_offset
        p <- p + ggplot2::geom_text(
          data = lab, ggplot2::aes(x = .data$x, y = .data$y, label = .data$txt),
          size = txt, colour = "black", family = base_family,
          vjust = if (identical(arc_side, "below")) 1 else 0
        )
      }
    }
  }

  if (!identical(group_label, NA) && !is.null(group_label %||% grp)) {
    p <- p + ggplot2::annotate(
      "text", x = xr[1], y = ymax * 1.45, label = group_label %||% grp,
      hjust = 0, size = txt, colour = "black", family = base_family
    )
  }

  psi_txt <- psi_label %||% psi_text(x, lid, grp)
  if (!is.null(psi_txt) && !is.na(psi_txt)) {
    p <- p + ggplot2::annotate(
      "text", x = xr[2] + diff(xr) * 0.03, y = ymax * 0.75, label = psi_txt,
      hjust = 0, size = txt, colour = "black", family = base_family
    )
  }

  # Head-room above the arcs: enough for the tallest arc plus its label.
  head <- max(1.62, max(arc_height_frac) * 1.30)
  ylim <- if (identical(arc_side, "below")) {
    c(-ymax * head, ymax * head)
  } else {
    c(0, ymax * head)
  }

  p +
    ggplot2::coord_cartesian(xlim = xlim, ylim = ylim, expand = FALSE) +
    resolve_theme(theme, base_size, base_family)
}

# Turn the `overlay` argument into a named list of panel -> groups. NULL gives
# one panel per group, which is the ordinary case.
#' @noRd
resolve_overlay <- function(overlay, groups, tracks) {
  if (is.null(overlay)) {
    return(stats::setNames(as.list(groups), groups))
  }
  if (is.list(overlay)) {
    unknown <- setdiff(unlist(overlay), groups)
    if (length(unknown)) {
      om_abort("{.arg overlay} names group{?s} not present: {.val {unknown}}.")
    }
    nm <- names(overlay) %||% vapply(overlay, paste, character(1), collapse = "+")
    return(stats::setNames(overlay, nm))
  }
  if (is.character(overlay) && length(overlay) == 1L) {
    if (!has_col(tracks, overlay)) {
      om_abort(c("{.arg overlay} column {.field {overlay}} is not in the tracks table.",
                 "i" = "Columns present: {.field {names(tracks)}}."))
    }
    key <- as.character(tracks[[overlay]])
    m <- unique(data.frame(g = tracks$group, k = key, stringsAsFactors = FALSE))
    m <- m[m$g %in% groups, , drop = FALSE]
    out <- split(m$g, m$k)
    return(out[order_levels(names(out))])
  }
  om_abort("{.arg overlay} must be a named list of groups or a column name.")
}

#' @noRd
psi_text <- function(x, lid, grp, digits = 3) {
  if (!nrow(x$psi)) return(NULL)
  v <- x$psi$psi[x$psi$locus_id == lid & x$psi$group %in% grp]
  if (!length(v)) return(NULL)
  if (is.na(v[1])) return("PSI = n/a")
  sprintf("PSI = %.*f", digits, v[1])
}

#' Draw a sashimi figure
#'
#' Stacks one coverage panel per group over a shared annotation panel. The
#' result is a `patchwork` object: print it, modify a panel with `[[`, or save
#' it with [save_sashimi()].
#'
#' @section Reading the figure:
#' Each panel is one group - a sample, a condition, a developmental stage. The
#' filled area is coverage; the arcs are junctions, drawn from donor to
#' acceptor with the supporting count printed at the apex. Where a `psi` slot
#' is present, the inclusion value for that group is printed in the right-hand
#' gutter. Beneath the stack sit the coordinate bar, the transcript models and
#' any annotation features.
#'
#' By default each panel is scaled to its own maximum, so a lowly expressed
#' group is still readable; set `fix_y_scale = TRUE` when relative depth
#' between panels is the point of the figure.
#'
#' @param x A `sashimi_data` object, or anything [as_sashimi_data()] accepts.
#' @param preset A named bundle of arguments setting the house style for a kind
#'   of figure; see [presets()]. Anything you pass explicitly overrides it, so
#'   `preset = "tss"` gets you the start-site style and
#'   `preset = "tss", arc_shape = "bezier"` gets you that style with a
#'   different curve.
#' @param locus A `locus_id` or gene name. Defaults to the first locus.
#' @param groups Groups to draw, in order. Defaults to every group present.
#' @param overlay Combine groups into shared panels instead of one panel each.
#'   A named list maps panel names to the groups they hold
#'   (`list(Ctrl = c("c1", "c2"), Treated = c("t1", "t2"))`); a single column
#'   name in the object's `tracks` slot groups by that column. This is
#'   `ggsashimi`'s `--overlay`, and pairs naturally with `alpha` below 1.
#' @param reverse_minus Draw a minus-strand locus 5\' to 3\', left to right, by
#'   flipping the x axis. MISO calls this `reverse_minus`; it makes the start
#'   site of a minus-strand gene appear on the left, where a reader expects it.
#' @param palette Palette name, colour vector, named colour vector, or the path
#'   to a one-colour-per-line palette file. See [omakase_palette()].
#' @param alpha Opacity of the coverage areas.
#' @param normalize Coverage normalisation, see [normalize_methods()].
#' @param library_sizes Named vector of library sizes for `normalize`.
#' @param aggregate Collapse replicates within each group; see
#'   [aggregate_tracks()].
#' @param fix_y_scale Give every panel the same y limit.
#' @param ymax Explicit y limit; overrides `fix_y_scale`.
#' @param log_y Draw coverage on a `log10(1 + value)` axis.
#' @param min_count Drop junctions supported by fewer than this many reads.
#' @param shrink Compress introns. See [intron_map()].
#' @param shrink_method,shrink_gamma,shrink_min Compression rule, its exponent,
#'   and the shortest intron worth compressing.
#' @param base_family Font family for every piece of text in the figure. The
#'   empty string uses the graphics device's default. Text geoms do not inherit
#'   a theme's family, so this is threaded to each of them explicitly.
#' @param arc_shape,arc_height_rule,arc_height_frac,arc_width_rule,arc_width,arc_side,arc_n
#'   Arc geometry; see [sashimi_track()].
#' @param show_arc_label,arc_label_format,arc_label_position,label_background,label_padding,label_offset
#'   Arc labelling: whether to print the count, how to format it, whether it
#'   sits `"on"` the arc (interrupting it) or `"above"` it, and the box padding
#'   and clearance. See [sashimi_track()].
#' @param background,background_alpha Fill and opacity for the coverage panels.
#'   `NA` leaves them transparent.
#' @param overlay_junction_fun How a junction shared by several overlaid groups
#'   is combined: `"mean"`, `"median"`, `"sum"` or `"max"`.
#' @param show_psi,psi_pad Whether to print PSI, and how much room to leave for
#'   it.
#' @param show_model,show_features,show_apex,show_tx_label,show_coord_bar,show_gene_label
#'   Which parts of the annotation panel to draw.
#' @param chrom_style How the contig name is printed: `"keep"`, `"ucsc"`
#'   (ensure a `chr` prefix) or `"ensembl"` (strip one).
#' @param coord_ticks Intermediate ticks on the coordinate bar; see
#'   [sashimi_annotation()]. Drawn by default whenever introns are compressed.
#' @param collapse_models Draw one row per role rather than one per transcript.
#' @param arrow_bins Strand arrowheads per transcript; `0` for none.
#' @param role_fill,feature_color Colours for transcript models and features.
#' @param base_size Base font size in points.
#' @param hairline Line width for rules.
#' @param panel_height,ann_height Relative heights of a coverage panel and the
#'   annotation panel.
#' @param title Figure title, or `NULL`.
#' @param theme A ggplot2 theme, or `NULL` for [theme_omakase()].
#'
#' @return A `patchwork` object.
#'
#' @examples
#' sd <- sashimi_data(
#'   loci = data.frame(locus_id = "a", gene_name = "DEMO", chrom = "chr1",
#'                     strand = "+", win_lo = 1000, win_hi = 2000),
#'   tracks = rbind(
#'     data.frame(locus_id = "a", group = "early", pos = seq(1000, 2000, 10),
#'                value = abs(sin(seq(0, pi, length.out = 101))) * 20),
#'     data.frame(locus_id = "a", group = "late", pos = seq(1000, 2000, 10),
#'                value = abs(sin(seq(0, pi, length.out = 101))) * 8)
#'   ),
#'   junctions = data.frame(locus_id = "a", group = c("early", "late"),
#'                          x0 = 1200, x1 = 1800, count = c(120, 30)),
#'   models = data.frame(locus_id = "a", tx_id = "tx1", role = "main",
#'                       start = c(1100, 1700), end = c(1300, 1900))
#' )
#' plot_sashimi(sd)
#'
#' @export
plot_sashimi <- function(x, preset = NULL, locus = NULL, groups = NULL,
                         overlay = NULL, overlay_junction_fun = "mean",
                         reverse_minus = FALSE,
                         palette = NULL, alpha = 1,
                         normalize = "none", library_sizes = NULL,
                         aggregate = "none",
                         fix_y_scale = FALSE, ymax = NULL, log_y = FALSE,
                         min_count = 0,
                         shrink = FALSE, shrink_method = "power",
                         shrink_gamma = 0.7, shrink_min = 100,
                         arc_shape = "sine", arc_height_rule = "auto",
                         arc_height_frac = c(0.80, 1.20),
                         arc_width_rule = "constant", arc_width = 0.5,
                         arc_side = "above", arc_n = 121,
                         show_arc_label = TRUE, arc_label_format = "activity",
                         arc_label_position = c("on", "above", "none"),
                         label_background = TRUE, label_padding = 0.6,
                         label_offset = 0.03,
                         background = NA, background_alpha = 1,
                         show_psi = TRUE, psi_pad = 0.30,
                         show_model = TRUE, show_features = TRUE,
                         show_apex = TRUE, show_tx_label = FALSE,
                         show_coord_bar = TRUE, show_gene_label = TRUE,
                         coord_ticks = TRUE, chrom_style = "keep",
                         collapse_models = FALSE, arrow_bins = 0,
                         role_fill = NULL, feature_color = OM_FEATURE_FILL,
                         base_size = 9, base_family = "", hairline = 0.3,
                         panel_height = 1, ann_height = 1.5,
                         title = NULL, theme = NULL) {
  # Merge the preset under whatever the caller actually named, then re-bind, so
  # every line below sees the resolved value and nothing has to test twice.
  if (!is.null(preset)) {
    explicit <- setdiff(names(as.list(match.call())[-1]), "preset")
    args <- apply_preset(preset, as.list(environment()), explicit)
    for (nm in names(args)) assign(nm, args[[nm]])
  }
  arc_label_position <- rlang::arg_match(arc_label_position,
                                         c("on", "above", "none"))

  x <- as_sashimi_data(x)
  validate_sashimi_data(x)

  if (!identical(normalize, "none")) {
    x <- normalize_tracks(x, normalize, library_sizes = library_sizes)
  }
  if (!identical(aggregate, "none")) {
    x <- aggregate_tracks(x, aggregate)
  }

  lid <- resolve_locus(x, locus)
  gi <- x$loci[x$loci$locus_id == lid, , drop = FALSE][1, ]
  cv <- x$tracks[x$tracks$locus_id == lid, , drop = FALSE]

  groups <- groups %||% order_levels(
    unique(c(cv$group, x$junctions$group[x$junctions$locus_id == lid])),
    x$meta$group_order
  )
  groups <- groups[!is.na(groups)]
  if (!length(groups)) om_abort("No groups to draw for locus {.val {lid}}.")

  map <- NULL
  if (isTRUE(shrink)) {
    iv <- introns_from_models(x$models, lid)
    map <- intron_map(iv, gi$win_lo, gi$win_hi, method = shrink_method,
                      gamma = shrink_gamma, min_intron = shrink_min)
  } else if (inherits(shrink, "omakase_intron_map")) {
    map <- shrink
  }

  cols <- resolve_palette(palette, groups)

  tf <- if (log_y) function(v) log10(1 + pmax(0, v)) else identity
  shared_ymax <- ymax
  if (is.null(shared_ymax) && fix_y_scale && nrow(cv)) {
    shared_ymax <- max(tf(cv$value), na.rm = TRUE)
  }

  xr <- compress_coords(c(gi$win_lo, gi$win_hi), map)
  pad <- if (show_psi && nrow(x$psi)) psi_pad else 0
  xlim <- c(xr[1], xr[2] + diff(xr) * pad)

  # Resolve the overlay spec into a list of panels, each a character vector of
  # the groups it holds. Without one, every group gets its own panel.
  panel_groups <- resolve_overlay(overlay, groups, x$tracks)

  panels <- lapply(names(panel_groups), function(pname) {
    g <- panel_groups[[pname]]
    sashimi_track(
      x, locus = lid, group = g, xlim = xlim, ymax = shared_ymax, map = map,
      fill = unname(cols[g]), alpha = alpha, base_size = base_size,
      base_family = base_family,
      group_label = pname,
      arc_shape = arc_shape, arc_height_rule = arc_height_rule,
      arc_height_frac = arc_height_frac, arc_width_rule = arc_width_rule,
      arc_width = arc_width, arc_side = arc_side, arc_n = arc_n,
      show_arc_label = show_arc_label, arc_label_format = arc_label_format,
      arc_label_position = arc_label_position,
      label_background = label_background, label_padding = label_padding,
      label_offset = label_offset,
      background = background, background_alpha = background_alpha,
      overlay_junction_fun = overlay_junction_fun, min_count = min_count,
      psi_label = if (show_psi) NULL else NA, psi_pad = pad,
      log_y = log_y, theme = theme
    )
  })

  parts <- panels
  heights <- rep(panel_height, length(panels))

  if (show_model || show_coord_bar || show_features) {
    ann <- sashimi_annotation(
      x, locus = lid, xlim = xlim, map = map, base_size = base_size,
      base_family = base_family, hairline = hairline, role_fill = role_fill,
      feature_color = feature_color, show_coord_bar = show_coord_bar,
      coord_ticks = coord_ticks, chrom_style = chrom_style,
      show_tx_label = show_tx_label,
      show_features = show_features && show_model,
      show_apex = show_apex, show_gene_label = show_gene_label,
      collapse = collapse_models, arrow_bins = arrow_bins, theme = theme
    )
    parts <- c(parts, list(ann))
    heights <- c(heights, ann_height)
  }

  if (isTRUE(reverse_minus) && identical(gi$strand, "-")) {
    # Reverse every panel's x axis so the locus reads 5' to 3' left to right.
    # coord_cartesian is already set on each panel, so the flip has to replace
    # it rather than be added alongside.
    parts <- lapply(parts, function(g) {
      g$coordinates <- ggplot2::coord_cartesian(
        xlim = rev(g$coordinates$limits$x),
        ylim = g$coordinates$limits$y, expand = FALSE
      )
      g
    })
  }

  p <- patchwork::wrap_plots(parts, ncol = 1, heights = heights)
  if (!is.null(title)) {
    p <- p + patchwork::plot_annotation(
      title = title,
      theme = ggplot2::theme(
        plot.title = ggplot2::element_text(size = base_size, colour = "black",
                                           face = "plain", hjust = 0,
                                           family = base_family)
      )
    )
  }
  # Carry the panel count so save_sashimi() can size the device without the
  # caller having to know how many groups were drawn.
  attr(p, "omakase_panels") <- length(panels)
  attr(p, "omakase_locus") <- lid
  p
}

#' Draw a sashimi figure for every locus
#'
#' Loops [plot_sashimi()] over the loci in an object, optionally writing one
#' file per locus.
#'
#' @param x A `sashimi_data` object.
#' @param dir Directory to write into, or `NULL` to return the plots without
#'   writing.
#' @param device File extension: `"pdf"`, `"png"`, `"svg"`, `"tiff"`, `"jpeg"`.
#' @param width Figure width in inches.
#' @param height Figure height in inches, or `NULL` to size it from the number
#'   of panels.
#' @param dpi Resolution for raster devices.
#' @param quiet Suppress the per-file progress message.
#' @param ... Passed to [plot_sashimi()].
#' @return A named list of plot objects, invisibly when writing to disk.
#' @examples
#' sd <- sashimi_data(
#'   loci = data.frame(locus_id = c("a", "b"), gene_name = c("A", "B"),
#'                     chrom = "chr1", strand = "+",
#'                     win_lo = 1, win_hi = 100),
#'   tracks = data.frame(locus_id = rep(c("a", "b"), each = 10), group = "g",
#'                       pos = rep(seq(1, 100, length.out = 10), 2),
#'                       value = runif(20))
#' )
#' plots <- plot_sashimi_all(sd)
#' names(plots)
#' @export
plot_sashimi_all <- function(x, dir = NULL, device = "pdf", width = 5.9,
                             height = NULL, dpi = 300, quiet = FALSE, ...) {
  x <- as_sashimi_data(x)
  ids <- loci(x)
  out <- stats::setNames(vector("list", length(ids)), ids)
  if (!is.null(dir)) dir.create(dir, showWarnings = FALSE, recursive = TRUE)

  for (id in ids) {
    p <- plot_sashimi(x, locus = id, ...)
    out[[id]] <- p
    if (!is.null(dir)) {
      nm <- x$loci$gene_name[x$loci$locus_id == id][1]
      nm <- if (is.na(nm) || !nzchar(nm)) id else nm
      f <- file.path(dir, paste0(sanitise_filename(nm), ".", device))
      save_sashimi(p, f, width = width, height = height, dpi = dpi)
    }
  }
  if (!is.null(dir)) {
    if (!quiet) om_inform("Wrote {length(ids)} figure{?s} to {.path {dir}}.")
    return(invisible(out))
  }
  out
}

#' Save a sashimi figure
#'
#' A thin wrapper over [ggplot2::ggsave()] that defaults the height to the one
#' the panel stack wants - `0.95` inches per coverage panel plus `2.0` for the
#' annotation - and picks a device from the file extension.
#'
#' @param plot A plot from [plot_sashimi()].
#' @param file Output path; the extension selects the device.
#' @param width Width in inches.
#' @param height Height in inches, or `NULL` to compute it.
#' @param dpi Resolution for raster devices.
#' @param ... Passed to [ggplot2::ggsave()].
#' @return `file`, invisibly.
#' @examples
#' sd <- sashimi_data(
#'   loci = data.frame(locus_id = "a", gene_name = "A", chrom = "chr1",
#'                     strand = "+", win_lo = 1, win_hi = 10),
#'   tracks = data.frame(locus_id = "a", group = "g", pos = 1:10, value = 1:10)
#' )
#' f <- file.path(tempdir(), "demo.pdf")
#' save_sashimi(plot_sashimi(sd), f)
#' @export
save_sashimi <- function(plot, file, width = 5.9, height = NULL, dpi = 300,
                         ...) {
  n <- attr(plot, "omakase_panels") %||% 3
  height <- height %||% (0.95 * n + 2.0)
  ggplot2::ggsave(file, plot, width = width, height = height, units = "in",
                  dpi = dpi, ...)
  invisible(file)
}

#' @noRd
sanitise_filename <- function(x) {
  x <- gsub("[/\\\\:*?\"<>|]+", "_", x)
  gsub("^[.]+", "", x)
}
