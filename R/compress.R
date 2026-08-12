# ---------------------------------------------------------------------------
# Intron compression.
#
# A gene is mostly intron, so a linear axis spends most of its width on empty
# sequence. Every sashimi tool shrinks introns somehow; ggsashimi raises the
# length to the power 0.7 and then rewrites the axis breaks by generating R
# source, and rmats2sashimiplot divides by an integer. omakase builds an
# explicit, invertible piecewise-linear map phi: genome -> plot, exposes it as
# an object, and hands ggplot2 a proper transform so the axis labels come out
# right by construction rather than by string surgery.
# ---------------------------------------------------------------------------

#' Intron compression methods
#'
#' @return A character vector of method names accepted by [intron_map()] and by
#'   the `shrink_method` argument of [plot_sashimi()].
#' @examples
#' shrink_methods()
#' @export
shrink_methods <- function() {
  c("none", "power", "log", "fixed", "scale")
}

# The length rule g(L): how long a real intron of length L is drawn.
#' @noRd
shrink_g <- function(L, method = "power", gamma = 0.7, cap = 200, scale = 5,
                     coef = 20) {
  method <- rlang::arg_match(method, shrink_methods())
  L <- pmax(0, as.numeric(L))
  out <- switch(method,
    none  = L,
    power = L^gamma,
    log   = coef * log1p(L),
    fixed = pmin(L, cap),
    scale = L / scale
  )
  # Never lengthen an intron, and never collapse it to nothing - a zero-width
  # intron would put a donor and an acceptor at the same plot coordinate and
  # the arc between them would vanish.
  pmin(L, pmax(out, pmin(L, 1)))
}

#' Build a genome-to-plot coordinate map that compresses introns
#'
#' Given a window and a set of intronic intervals, constructs a monotone
#' piecewise-linear map \eqn{\phi} that leaves exonic sequence at scale 1 and
#' draws each intron of length \eqn{L = b - a} at a reduced length \eqn{g(L)}:
#'
#' \describe{
#'   \item{`none`}{\eqn{g(L) = L}}
#'   \item{`power`}{\eqn{g(L) = L^{\gamma}}, default \eqn{\gamma = 0.7}
#'     (the `ggsashimi` rule)}
#'   \item{`log`}{\eqn{g(L) = c \log(1 + L)}}
#'   \item{`fixed`}{\eqn{g(L) = \min(L, k)}}
#'   \item{`scale`}{\eqn{g(L) = L / s\,} (the `rmats2sashimiplot --intron_s`
#'     rule)}
#' }
#'
#' Writing the cumulative shift after the \eqn{j}th intron as
#' \eqn{\Delta_j = \sum_{i \le j} (L_i - g(L_i))}, the map is
#' \eqn{\phi(x) = x - \Delta_{j(x)}} on exonic stretches, and interpolates
#' linearly across each compressed intron. Because \eqn{\phi} is strictly
#' increasing it has an exact inverse, which is what [plot_sashimi()] uses to
#' label the axis with true genomic coordinates.
#'
#' @param introns A two-column data frame or matrix of intron `start`/`end`
#'   coordinates, or `NULL` for no compression. Overlapping intervals are
#'   merged; intervals are sorted internally.
#' @param win_lo,win_hi The plotted window. Introns are clipped to it.
#' @param method One of [shrink_methods()].
#' @param gamma Exponent for `method = "power"`.
#' @param cap Maximum drawn intron length for `method = "fixed"`.
#' @param scale Divisor for `method = "scale"`.
#' @param coef Multiplier \eqn{c} for `method = "log"`.
#' @param min_intron Introns shorter than this are left uncompressed; shrinking
#'   a 60 bp intron buys no width and distorts short-intron genes.
#' @param exon_scale Divisor applied to exonic stretches, the counterpart of
#'   `rmats2sashimiplot`'s `--exon_s`. The default of `1` leaves exons at true
#'   scale; larger values shrink them too, which is occasionally wanted when a
#'   single huge exon dominates the view.
#'
#' @return An object of class `omakase_intron_map` with elements `genome` and
#'   `plot` (the knot coordinates), plus the settings used. Apply it with
#'   [compress_coords()] and invert it with [expand_coords()].
#'
#' @examples
#' m <- intron_map(data.frame(start = c(120, 400), end = c(300, 900)),
#'                 win_lo = 100, win_hi = 1000)
#' compress_coords(c(100, 300, 1000), m)
#' # phi is exactly invertible
#' expand_coords(compress_coords(c(150, 950), m), m)
#'
#' @export
intron_map <- function(introns, win_lo, win_hi, method = "power", gamma = 0.7,
                       cap = 200, scale = 5, coef = 20, min_intron = 100,
                       exon_scale = 1) {
  method <- rlang::arg_match(method, shrink_methods())

  identity_map <- function() {
    # With exons scaled the map is still linear, just not at scale 1, so it
    # remains exactly invertible.
    structure(
      list(genome = c(win_lo, win_hi),
           plot = c(win_lo, win_lo + (win_hi - win_lo) / exon_scale),
           method = if (exon_scale == 1) "none" else "exon_only",
           exon_scale = exon_scale, introns = NULL),
      class = "omakase_intron_map"
    )
  }
  if (identical(method, "none") || is.null(introns) || !NROW(introns)) {
    return(identity_map())
  }

  iv <- as_df(introns)
  names(iv)[1:2] <- c("start", "end")
  iv <- iv[stats::complete.cases(iv$start, iv$end), c("start", "end"), drop = FALSE]
  # Clip to the window, then merge overlaps: two transcripts usually share most
  # introns, and compressing the same sequence twice would break monotonicity.
  iv$start <- pmax(iv$start, win_lo)
  iv$end <- pmin(iv$end, win_hi)
  iv <- iv[iv$end - iv$start >= min_intron, , drop = FALSE]
  if (!nrow(iv)) return(identity_map())

  iv <- iv[order(iv$start), , drop = FALSE]
  merged <- list(c(iv$start[1], iv$end[1]))
  if (nrow(iv) > 1) {
    for (i in 2:nrow(iv)) {
      last <- merged[[length(merged)]]
      if (iv$start[i] <= last[2]) {
        merged[[length(merged)]] <- c(last[1], max(last[2], iv$end[i]))
      } else {
        merged[[length(merged)]] <- last
        merged <- c(merged, list(c(iv$start[i], iv$end[i])))
      }
    }
  }
  iv <- do.call(rbind, merged)
  colnames(iv) <- c("start", "end")

  # Walk the window laying down knots: exonic stretches keep their length,
  # intronic ones are drawn at g(L).
  gpos <- win_lo
  ppos <- win_lo
  gk <- gpos
  pk <- ppos
  for (i in seq_len(nrow(iv))) {
    a <- iv[i, "start"]
    b <- iv[i, "end"]
    if (a > gpos) {
      ppos <- ppos + (a - gpos) / exon_scale
      gpos <- a
      gk <- c(gk, gpos); pk <- c(pk, ppos)
    }
    L <- b - a
    ppos <- ppos + shrink_g(L, method, gamma = gamma, cap = cap,
                            scale = scale, coef = coef)
    gpos <- b
    gk <- c(gk, gpos); pk <- c(pk, ppos)
  }
  if (win_hi > gpos) {
    pk <- c(pk, ppos + (win_hi - gpos) / exon_scale)
    gk <- c(gk, win_hi)
  }

  structure(
    list(genome = gk, plot = pk, method = method, gamma = gamma, cap = cap,
         scale = scale, coef = coef, exon_scale = exon_scale,
         introns = as.data.frame(iv)),
    class = "omakase_intron_map"
  )
}

#' Apply or invert an intron compression map
#'
#' `compress_coords()` is \eqn{\phi}: genome coordinates to plot coordinates.
#' `expand_coords()` is \eqn{\phi^{-1}}. Values outside the mapped window are
#' extrapolated at scale 1, so an arc anchored just beyond the window still
#' lands somewhere sensible.
#'
#' @param x Numeric vector of coordinates.
#' @param map An `omakase_intron_map` from [intron_map()].
#' @return A numeric vector the same length as `x`.
#' @examples
#' m <- intron_map(data.frame(start = 200, end = 800), 100, 1000)
#' compress_coords(c(100, 500, 1000), m)
#' expand_coords(compress_coords(500, m), m)
#' @export
compress_coords <- function(x, map) {
  if (is.null(map) || identical(map$method, "none")) return(as.numeric(x))
  if (identical(map$method, "exon_only")) {
    return(map$plot[1] + (as.numeric(x) - map$genome[1]) / map$exon_scale)
  }
  interp(as.numeric(x), map$genome, map$plot)
}

#' @rdname compress_coords
#' @export
expand_coords <- function(x, map) {
  if (is.null(map) || identical(map$method, "none")) return(as.numeric(x))
  if (identical(map$method, "exon_only")) {
    return(map$genome[1] + (as.numeric(x) - map$plot[1]) * map$exon_scale)
  }
  interp(as.numeric(x), map$plot, map$genome)
}

# Monotone linear interpolation with linear (slope-1 at the ends)
# extrapolation. stats::approx returns NA outside the range, which would drop
# arcs anchored beyond the window.
#' @noRd
interp <- function(x, from, to) {
  out <- stats::approx(from, to, xout = x, rule = 1)$y
  lo <- !is.na(x) & x < from[1]
  hi <- !is.na(x) & x > from[length(from)]
  if (any(lo)) out[lo] <- to[1] - (from[1] - x[lo])
  if (any(hi)) {
    n <- length(from)
    out[hi] <- to[n] + (x[hi] - from[n])
  }
  out
}

#' @export
print.omakase_intron_map <- function(x, ...) {
  if (identical(x$method, "none")) {
    cli::cli_text("<omakase intron map> identity (no compression)")
    return(invisible(x))
  }
  gw <- diff(range(x$genome))
  pw <- diff(range(x$plot))
  cli::cli_text(
    "<omakase intron map> method {.val {x$method}}, {nrow(x$introns)} intron{?s}, \\
     {round(gw)} bp drawn as {round(pw)} ({round(100 * pw / gw)}%)"
  )
  invisible(x)
}

#' A ggplot2 axis transform for compressed genomic coordinates
#'
#' Wraps an [intron_map()] as a `scales` transform so that
#' `scale_x_continuous(transform = ...)` places tick marks at round *genomic*
#' coordinates while drawing the data in compressed space. This is what keeps
#' the axis honest when introns are shrunk.
#'
#' @param map An `omakase_intron_map`.
#' @param n Target number of axis breaks.
#' @return A transform object from [scales::new_transform()].
#' @examples
#' m <- intron_map(data.frame(start = 200, end = 8000), 100, 10000)
#' tr <- intron_trans(m)
#' tr$breaks(c(100, 10000))
#' @export
intron_trans <- function(map, n = 5) {
  force(map); force(n)
  scales::new_transform(
    name = "omakase-intron",
    transform = function(x) compress_coords(x, map),
    inverse = function(x) expand_coords(x, map),
    breaks = function(limits, ...) {
      # Choose pretty breaks in *plot* space, so they are visually even, then
      # report them in genome space; this is the same idea ggsashimi implements
      # by post-processing generated R code.
      pl <- compress_coords(range(limits, na.rm = TRUE), map)
      b <- labeling::extended(pl[1], pl[2], m = n, only.loose = FALSE)
      g <- expand_coords(b, map)
      g[g >= min(limits, na.rm = TRUE) & g <= max(limits, na.rm = TRUE)]
    },
    format = function(x, ...) format_coord(x)
  )
}

# Collect the introns implied by a set of exon models, per locus. Used when the
# caller asks for shrinking but supplies models rather than explicit introns.
#' @noRd
introns_from_models <- function(models, locus_id = NULL) {
  if (is.null(models) || !nrow(models)) return(NULL)
  if (!is.null(locus_id)) {
    models <- models[models$locus_id %in% locus_id, , drop = FALSE]
  }
  if (!nrow(models)) return(NULL)
  parts <- lapply(split(models, models$tx_id), function(d) {
    d <- d[order(d$start), , drop = FALSE]
    if (nrow(d) < 2) return(NULL)
    data.frame(start = utils::head(d$end, -1) + 1,
               end = utils::tail(d$start, -1) - 1)
  })
  out <- rbind_all(parts)
  if (is.null(out) || !nrow(out)) return(NULL)
  out[out$end > out$start, , drop = FALSE]
}
