# ---------------------------------------------------------------------------
# Junction arc geometry.
#
# Every sashimi tool draws its arcs differently, and the choice is visible: the
# sine arc reads as a clean pointer, the x-spline hugs the baseline and looks
# more like IGV, the Bezier is what MISO drew. omakase exposes all of them
# through one function that returns a plain path, so arcs are ordinary ggplot
# layers rather than grid grobs bolted on with annotation_custom().
# ---------------------------------------------------------------------------

#' Arc shapes available for junction curves
#'
#' @return A character vector of shape names accepted by [arc_path()] and by
#'   the `arc_shape` argument of [plot_sashimi()].
#' @examples
#' arc_shapes()
#' @export
arc_shapes <- function() {
  c("sine", "parabola", "bezier", "xspline", "elbow", "arch")
}

#' Build the path of a junction arc
#'
#' Returns the polyline for a single arc running from `x0` to `x1` and peaking
#' at height `h`. All shapes start and end at `y = 0` and reach exactly `h` at
#' their apex, so arcs of different shapes are directly comparable.
#'
#' The shapes, parameterised by \eqn{t \in [0, 1]} with
#' \eqn{x(t) = x_0 + (x_1 - x_0)t}:
#'
#' \describe{
#'   \item{`sine`}{\eqn{y = h \sin(\pi t)}. The default: symmetric, with a
#'     gentle take-off from the baseline.}
#'   \item{`parabola`}{\eqn{y = 4 h\, t (1 - t)}. Slightly fuller shoulders
#'     than the sine.}
#'   \item{`bezier`}{Cubic Bezier
#'     \eqn{B(t) = (1-t)^3 P_0 + 3t(1-t)^2 P_1 + 3t^2(1-t) P_2 + t^3 P_3}
#'     with control points lifted to \eqn{4h/3} so the apex lands on \eqn{h};
#'     this is the curve MISO's `sashimi_plot` draws.}
#'   \item{`xspline`}{An approximation of the X-spline used by `ggsashimi`,
#'     which rises steeply from the donor and flattens across the middle.}
#'   \item{`elbow`}{Straight risers joined by a flat top - useful when many
#'     arcs overlap and curvature becomes noise.}
#'   \item{`arch`}{A semi-ellipse, \eqn{y = h\sqrt{1 - (2t - 1)^2}}, the
#'     roundest of the set.}
#' }
#'
#' @param x0,x1 Genomic start and end of the junction.
#' @param h Apex height, in the units of the coverage axis.
#' @param shape One of [arc_shapes()].
#' @param n Number of points in the returned path. Higher is smoother; 120 is
#'   plenty at print size, and `elbow` ignores it.
#' @param y0 Baseline the arc springs from. Non-zero when arcs are stacked
#'   above a coverage track that does not start at zero.
#'
#' @return A data frame with columns `x` and `y`.
#'
#' @examples
#' head(arc_path(100, 200, h = 10))
#' # every shape peaks at h
#' vapply(arc_shapes(), function(s) max(arc_path(0, 1, 5, shape = s)$y), numeric(1))
#'
#' @export
arc_path <- function(x0, x1, h, shape = "sine", n = 120, y0 = 0) {
  shape <- rlang::arg_match(shape, arc_shapes())
  if (!is.finite(x0) || !is.finite(x1) || !is.finite(h)) {
    return(data.frame(x = numeric(0), y = numeric(0)))
  }
  if (identical(shape, "elbow")) {
    # Four corners: up, across, down. Drawn without interpolation so the
    # verticals stay vertical at any device size.
    return(data.frame(
      x = c(x0, x0, x1, x1),
      y = c(y0, y0 + h, y0 + h, y0)
    ))
  }

  # An even number of samples straddles t = 0.5 and so never evaluates the
  # apex, leaving the drawn arc a hair short of h. Forcing an odd count puts a
  # point exactly on the peak.
  n <- max(3L, as.integer(n))
  if (n %% 2L == 0L) n <- n + 1L
  t <- seq(0, 1, length.out = n)
  x <- x0 + (x1 - x0) * t
  y <- switch(shape,
    sine     = h * sin(pi * t),
    parabola = 4 * h * t * (1 - t),
    arch     = h * sqrt(pmax(0, 1 - (2 * t - 1)^2)),
    bezier   = {
      # Control points at 4h/3 make the cubic peak at exactly h when t = 0.5,
      # since B(0.5) = (P0 + 3P1 + 3P2 + P3)/8 = (0 + 3k + 3k + 0)/8 = 3k/4.
      k <- 4 * h / 3
      3 * t * (1 - t)^2 * k + 3 * t^2 * (1 - t) * k
    },
    xspline  = {
      # ggsashimi joins two xsplineGrob halves that leave the baseline almost
      # vertically and flatten across the top. A raised-cosine on a
      # square-rooted parameter reproduces that profile closely without
      # dropping to grid.
      u <- sqrt(sin(pi * t))
      h * u
    }
  )
  # Snap the ends to the baseline: sin(pi) and the xspline approximation both
  # leave a float's worth of residue, which shows up as a hairline stub where
  # the arc meets the coverage.
  y[c(1L, length(y))] <- 0
  data.frame(x = x, y = y0 + y)
}

#' Arc height rules
#'
#' How a junction's count is turned into the height of its arc.
#'
#' @return A character vector of rule names.
#' @examples
#' arc_height_rules()
#' @export
arc_height_rules <- function() {
  c("auto", "constant", "span", "linear", "sqrt", "log")
}

#' Compute arc apex heights from junction counts
#'
#' @param count Numeric vector of junction counts or activities.
#' @param rule One of [arc_height_rules()]. `"constant"` gives every arc the
#'   same height, which keeps the picture legible when counts span orders of
#'   magnitude; `"span"` scales height with the width of the junction, so
#'   nested junctions nest visually instead of crossing; the rest scale height
#'   with count.
#' @param min_h,max_h The height range to map onto, in coverage-axis units.
#' @param base Logarithm base used by `rule = "log"`.
#' @param span Junction widths, required for `rule = "span"`.
#'
#' @return A numeric vector of heights, the same length as `count`.
#'
#' @examples
#' arc_heights(c(1, 10, 100), rule = "log", min_h = 1, max_h = 5)
#'
#' @export
arc_heights <- function(count, rule = "constant", min_h = 0, max_h = 1,
                        base = 10, span = NULL) {
  rule <- rlang::arg_match(rule, arc_height_rules())
  count <- as.numeric(count)
  if (!length(count)) return(numeric(0))
  if (rule %in% c("constant", "auto")) return(rep(max_h, length(count)))
  if (identical(rule, "span") && is.null(span)) {
    om_abort('{.arg span} is required when {.code rule = "span"}.')
  }

  z <- switch(rule,
    span   = sqrt(pmax(0, as.numeric(span))),
    linear = count,
    sqrt   = sqrt(pmax(0, count)),
    log    = log(pmax(0, count) + 1, base = base)
  )
  rng <- range(z, na.rm = TRUE)
  if (!all(is.finite(rng)) || diff(rng) == 0) {
    return(rep((min_h + max_h) / 2, length(count)))
  }
  min_h + (z - rng[1]) / diff(rng) * (max_h - min_h)
}

#' Compute arc line widths from junction counts
#'
#' The classic sashimi convention is that a junction supported by more reads is
#' drawn with a thicker curve. The default here is
#' \deqn{w = w_0 \left(\log_b (n + 1)\right)^{\alpha}}
#' clamped to `range`. Passing `rule = "miso"` instead reproduces MISO's
#' \eqn{w = \log(b)\,(n+1)^{0.33} \times 0.1}, and `rule = "constant"` disables
#' the effect entirely.
#'
#' @param count Numeric vector of junction counts.
#' @param rule `"log"` (the default), `"miso"`, `"linear"`, or `"constant"`.
#' @param w0 Scale factor \eqn{w_0}.
#' @param base Logarithm base \eqn{b}.
#' @param alpha Exponent \eqn{\alpha} applied to the log term.
#' @param range Numeric length-2 clamp applied to the result.
#'
#' @return A numeric vector of line widths.
#'
#' @examples
#' arc_widths(c(1, 5, 50, 500))
#' arc_widths(c(1, 5, 50, 500), rule = "constant", w0 = 0.5)
#'
#' @export
arc_widths <- function(count, rule = c("log", "miso", "linear", "constant"),
                       w0 = 0.5, base = 10, alpha = 1, range = c(0.15, 2)) {
  rule <- match.arg(rule)
  count <- as.numeric(count)
  if (!length(count)) return(numeric(0))
  if (identical(rule, "constant")) return(rep(w0, length(count)))

  w <- switch(rule,
    log    = w0 * (log(pmax(0, count) + 1, base = base))^alpha,
    miso   = log(base) * (pmax(0, count) + 1)^0.33 * 0.1,
    linear = {
      m <- max(count, na.rm = TRUE)
      if (!is.finite(m) || m == 0) rep(w0, length(count)) else w0 * count / m
    }
  )
  w[!is.finite(w)] <- w0
  pmin(pmax(w, range[1]), range[2])
}

# Build the full set of arc paths for one panel, staggering the heights of arcs
# that share an endpoint so their labels cannot collide. This is the behaviour
# a start-site figure relies on when a main TSS and an alternative TSS both
# point at the same downstream anchor.
#' @noRd
build_arcs <- function(j, ymax, shape = "sine", height_rule = "constant",
                       height_frac = c(0.80, 1.20), n = 120, side = "above",
                       stagger = TRUE) {
  if (!nrow(j)) {
    return(list(paths = NULL, labels = NULL))
  }
  j <- j[order(j$x0, j$x1), , drop = FALSE]
  j$.id <- seq_len(nrow(j))
  # Deliberately not called `n`: that is the per-arc point count parameter, and
  # shadowing it draws every arc as a three-point triangle.
  n_arc <- nrow(j)

  # "auto" keeps the staggered constant height that reads best for the two or
  # three arcs of a designed figure, and switches to span-scaling once there
  # are enough junctions that equal-height arcs would cross into a thicket.
  if (identical(height_rule, "auto")) {
    height_rule <- if (n_arc > length(height_frac) + 1L) "span" else "constant"
  }

  if (identical(height_rule, "constant")) {
    if (!is.null(names(height_frac)) && any(nzchar(names(height_frac)))) {
      # Named fractions are keyed by role, so a given transcript always draws
      # at the same height across every panel of a figure. That is what makes a
      # stack of stages readable: the reader learns "the tall arc is the
      # alternative site" once, rather than per panel.
      frac <- unname(height_frac[as.character(j$role)])
      frac[is.na(frac)] <- height_frac[[1]]
    } else if (stagger && n_arc > 1L) {
      # Otherwise cycle the fractions so neighbours differ in height. Arcs
      # sharing an end point go first in the cycle, since those are the ones
      # whose labels would otherwise land on top of each other. Either end may
      # be the shared one - junctions fanning out from a common donor and
      # junctions converging on a common acceptor are both routine.
      key0 <- paste(j$x0, j$group)
      key1 <- paste(j$x1, j$group)
      dup <- function(k) stats::ave(rep(1L, n_arc), k, FUN = sum) > 1L
      key <- ifelse(dup(key0), key0, key1)
      rank_in_key <- stats::ave(seq_len(n_arc), key, FUN = seq_along)
      shared <- dup(key0) | dup(key1)
      pos <- ifelse(shared, rank_in_key, seq_len(n_arc))
      frac <- height_frac[((pos - 1L) %% length(height_frac)) + 1L]
    } else {
      frac <- rep(height_frac[1], n_arc)
    }
    h <- ymax * frac
  } else {
    h <- arc_heights(j$count, rule = height_rule,
                     min_h = ymax * min(height_frac),
                     max_h = ymax * max(height_frac),
                     span = j$x1 - j$x0)
  }
  j$.h <- h

  sgn <- if (identical(side, "below")) -1 else 1
  paths <- rbind_all(lapply(seq_len(nrow(j)), function(i) {
    p <- arc_path(j$x0[i], j$x1[i], sgn * j$.h[i], shape = shape, n = n)
    if (!nrow(p)) return(NULL)
    p$.id <- j$.id[i]
    p$role <- j$role[i]
    p$group <- j$group[i]
    p$count <- j$count[i]
    p
  }))
  labels <- data.frame(
    x = (j$x0 + j$x1) / 2,
    y = sgn * j$.h * 1.02,
    count = j$count,
    role = j$role,
    group = j$group,
    # Carried through rather than looked up by the caller: this function sorts
    # its input, so an index into the original table no longer lines up.
    label = col(j, "label") %||% NA_character_,
    stringsAsFactors = FALSE
  )
  list(paths = paths, labels = labels, heights = j$.h, junctions = j)
}
