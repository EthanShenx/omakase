# ---------------------------------------------------------------------------
# Coverage normalisation.
#
# Neither ggsashimi nor rmats2sashimiplot normalises: they plot raw depth, so a
# deeply sequenced library looks like a highly expressed gene. omakase makes
# the choice explicit and offers the scalings people actually use, including
# DESeq2's median-of-ratios, which is what start-site activity tracks are
# usually built on.
# ---------------------------------------------------------------------------

#' Coverage normalisation methods
#'
#' @return A character vector of method names accepted by [normalize_tracks()]
#'   and the `normalize` argument of [plot_sashimi()].
#' @examples
#' normalize_methods()
#' @export
normalize_methods <- function() {
  c("none", "cpm", "rpm", "rpkm", "size_factor", "manual", "max", "sum")
}

#' Normalise coverage tracks
#'
#' Rescales the `value` column of a `sashimi_data` object's `tracks` slot (and,
#' where it makes sense, the matching junction counts) so panels from libraries
#' of different depth are comparable.
#'
#' With \eqn{c_i} the raw signal in bin \eqn{i} of sample \eqn{j}, \eqn{N_j}
#' the library size and \eqn{L} the bin width in kilobases:
#'
#' \describe{
#'   \item{`none`}{\eqn{c_i}}
#'   \item{`cpm`, `rpm`}{\eqn{c_i \times 10^6 / N_j}}
#'   \item{`rpkm`}{\eqn{c_i \times 10^9 / (N_j L)}}
#'   \item{`size_factor`}{\eqn{c_i / s_j}, with DESeq2's median-of-ratios
#'     factor \eqn{s_j = \mathrm{median}_i \left( k_{ij} / (\prod_v k_{iv})^{1/m} \right)}
#'     over the \eqn{m} samples, taken across bins with non-zero signal in
#'     every sample}
#'   \item{`manual`}{\eqn{c_i / f_j} for user-supplied factors `factors`}
#'   \item{`max`}{each sample scaled to a maximum of 1}
#'   \item{`sum`}{each sample scaled to sum to 1}
#' }
#'
#' @param x A `sashimi_data` object.
#' @param method One of [normalize_methods()].
#' @param library_sizes Named numeric vector of \eqn{N_j}, named by sample. If
#'   absent for `cpm`/`rpm`/`rpkm`, the sum of the plotted signal is used
#'   instead, and a note is emitted - a window-local total is not a library
#'   size, and only makes panels comparable in a rough sense.
#' @param factors Named numeric vector of divisors for `method = "manual"`.
#' @param bin_width Bin width in base pairs, needed for `rpkm`. Taken from the
#'   track spacing when `NULL`.
#' @param scale_junctions If `TRUE`, junction counts are divided by the same
#'   per-sample factor as the coverage, so arc labels stay on the same scale as
#'   the track beneath them.
#'
#' @return The `sashimi_data` object with rescaled values and a record of the
#'   scaling in `meta`.
#'
#' @examples
#' sd <- sashimi_data(
#'   loci = data.frame(locus_id = "a", gene_name = "A", chrom = "chr1",
#'                     strand = "+", win_lo = 1, win_hi = 100),
#'   tracks = data.frame(locus_id = "a", group = "g", sample = "s1",
#'                       pos = 1:100, value = 1)
#' )
#' normalize_tracks(sd, "cpm")$tracks$value[1]
#'
#' @export
normalize_tracks <- function(x, method = "none", library_sizes = NULL,
                             factors = NULL, bin_width = NULL,
                             scale_junctions = TRUE) {
  method <- rlang::arg_match(method, normalize_methods())
  if (identical(method, "none") || !nrow(x$tracks)) {
    return(x)
  }
  tr <- x$tracks
  # Where samples are not distinguished, everything shares one factor.
  key <- if (all(is.na(tr$sample))) tr$group else tr$sample
  key <- as.character(key)
  samples <- unique(key)

  bw <- bin_width %||% infer_bin_width(tr)

  div <- switch(method,
    manual = {
      if (is.null(factors)) om_abort('{.arg factors} is required when {.code method = "manual"}.')
      lookup_factors(factors, samples, "factors")
    },
    size_factor = deseq_size_factors(tr, key),
    max = vapply(samples, function(s) max(tr$value[key == s], na.rm = TRUE), numeric(1)),
    sum = vapply(samples, function(s) sum(tr$value[key == s], na.rm = TRUE), numeric(1)),
    {
      # cpm / rpm / rpkm all need a library size.
      N <- if (!is.null(library_sizes)) {
        lookup_factors(library_sizes, samples, "library_sizes")
      } else {
        om_inform(c(
          "!" = "No {.arg library_sizes} given; using the total signal in the plotted window.",
          "i" = "Pass real library sizes for a normalisation that is comparable between genes."
        ))
        vapply(samples, function(s) sum(tr$value[key == s], na.rm = TRUE), numeric(1))
      }
      if (identical(method, "rpkm")) N * bw / 1e3 / 1e6 else N / 1e6
    }
  )
  div[!is.finite(div) | div == 0] <- 1

  f <- div[match(key, samples)]
  x$tracks$value <- tr$value / f

  if (scale_junctions && nrow(x$junctions) &&
      method %in% c("cpm", "rpm", "rpkm", "size_factor", "manual")) {
    jkey <- as.character(if (all(is.na(x$junctions$sample))) x$junctions$group
                         else x$junctions$sample)
    jf <- div[match(jkey, samples)]
    jf[is.na(jf)] <- 1
    x$junctions$count <- x$junctions$count / jf
  }

  x$meta$normalize <- method
  x$meta$normalize_factors <- stats::setNames(round(div, 6), samples)
  x
}

#' @noRd
lookup_factors <- function(v, samples, what) {
  if (is.null(names(v))) {
    v <- recycle_to(as.numeric(v), length(samples), what)
    return(stats::setNames(v, samples))
  }
  missing <- setdiff(samples, names(v))
  if (length(missing)) {
    om_abort("{what} is missing entr{?y/ies} for {.val {missing}}.")
  }
  as.numeric(v[samples])
}

# Median spacing between consecutive bin positions; robust to a locus whose
# window is shorter than the rest.
#' @noRd
infer_bin_width <- function(tr) {
  d <- unlist(lapply(split(tr$pos, paste(tr$locus_id, tr$group)), function(p) {
    p <- sort(unique(p))
    if (length(p) < 2) return(NULL)
    diff(p)
  }), use.names = FALSE)
  if (!length(d)) return(1)
  stats::median(d, na.rm = TRUE)
}

# DESeq2 median-of-ratios, computed on the binned coverage matrix. Bins with a
# zero in any sample drop out (the geometric mean would be zero), matching the
# reference implementation.
#' @noRd
deseq_size_factors <- function(tr, key) {
  samples <- unique(key)
  if (length(samples) < 2) return(stats::setNames(rep(1, length(samples)), samples))

  bin_id <- paste(tr$locus_id, tr$pos, sep = ":")
  mat <- tapply(tr$value, list(bin_id, key), sum, default = NA_real_)
  mat <- mat[, samples, drop = FALSE]
  ok <- stats::complete.cases(mat) & apply(mat > 0, 1, all)
  if (sum(ok) < 2) {
    om_warn("Too few bins with signal in every sample for size factors; falling back to no scaling.")
    return(stats::setNames(rep(1, length(samples)), samples))
  }
  m <- mat[ok, , drop = FALSE]
  gm <- exp(rowMeans(log(m)))
  sf <- apply(m / gm, 2, stats::median, na.rm = TRUE)
  sf[!is.finite(sf) | sf <= 0] <- 1
  stats::setNames(as.numeric(sf), samples)
}

#' Aggregate replicate tracks into one track per group
#'
#' Replicates are usually plotted one panel each, but with many samples that is
#' unreadable. Aggregation collapses each group to a single track, the way
#' `ggsashimi`'s `--aggr` does, with the option to keep junction counts summed
#' while the coverage is averaged.
#'
#' @param x A `sashimi_data` object.
#' @param fun One of `"mean"`, `"median"`, `"sum"`, `"max"`, or `"none"` to
#'   leave the object untouched.
#' @param junction_fun Aggregation applied to junction counts; defaults to the
#'   same as `fun`. `"sum"` is the honest choice when arcs should reflect total
#'   support across replicates.
#' @return The aggregated `sashimi_data` object.
#' @examples
#' sd <- sashimi_data(
#'   loci = data.frame(locus_id = "a", gene_name = "A", chrom = "chr1",
#'                     strand = "+", win_lo = 1, win_hi = 3),
#'   tracks = data.frame(locus_id = "a", group = "g",
#'                       sample = rep(c("s1", "s2"), each = 3),
#'                       pos = rep(1:3, 2), value = c(1, 2, 3, 3, 4, 5))
#' )
#' aggregate_tracks(sd, "mean")$tracks
#' @export
aggregate_tracks <- function(x, fun = c("mean", "median", "sum", "max", "none"),
                             junction_fun = NULL) {
  fun <- match.arg(fun)
  if (identical(fun, "none")) return(x)
  f <- match.fun(fun)
  jf <- match.fun(junction_fun %||% fun)

  if (nrow(x$tracks)) {
    tr <- x$tracks
    agg <- stats::aggregate(
      tr$value,
      by = list(locus_id = tr$locus_id, group = tr$group, pos = tr$pos,
                strand = tr$strand),
      FUN = function(v) f(v, na.rm = TRUE)
    )
    names(agg)[ncol(agg)] <- "value"
    agg$sample <- NA_character_
    agg$bin <- NA_real_
    x$tracks <- normalise_slot(agg, "tracks")
  }
  if (nrow(x$junctions)) {
    jn <- x$junctions
    agg <- stats::aggregate(
      jn$count,
      by = list(locus_id = jn$locus_id, group = jn$group, x0 = jn$x0,
                x1 = jn$x1, role = ifelse(is.na(jn$role), "", jn$role),
                strand = jn$strand),
      FUN = function(v) jf(v, na.rm = TRUE)
    )
    names(agg)[ncol(agg)] <- "count"
    agg$role[agg$role == ""] <- NA_character_
    agg$sample <- NA_character_
    agg$label <- NA_character_
    x$junctions <- normalise_slot(agg, "junctions")
  }
  x$meta$aggregate <- fun
  x
}
