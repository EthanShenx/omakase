# ---------------------------------------------------------------------------
# 5'-end tag data.
#
# CAGE, STRT, CamoTSS and related protocols report a single base per molecule -
# the transcription start site - rather than a read that covers an exon. A raw
# tag track is therefore a row of spikes, not something that reads as coverage,
# and none of the existing sashimi tools will draw it at all.
#
# The recipe here is the one behind 5'-tag activity tracks: bin the tags,
# scale each library, extend each tag toward the 3' end so the track has body
# without the start site moving, then average across the samples of a group.
# ---------------------------------------------------------------------------

#' Read a 5'-tag BED file
#'
#' Reads single-base transcription start site tags. Any BED-like file works:
#' columns 1-3 are contig, start, end, column 5 is taken as a count if present,
#' and column 6 as strand.
#'
#' @param path Path to the BED file, optionally gzipped.
#' @param region Optional region to restrict to.
#' @param strand Keep only tags on this strand (`"+"`, `"-"`, or `NULL` for
#'   both).
#' @param count_col Column holding a per-tag count. `NULL` treats every row as
#'   one tag.
#'
#' @return A data frame with `chrom`, `pos`, `count`, `strand`.
#'
#' @examples
#' f <- tempfile()
#' writeLines(c("chr1\t999\t1000\ttag1\t3\t+"), f)
#' read_tag_bed(f)
#'
#' @export
read_tag_bed <- function(path, region = NULL, strand = NULL, count_col = 5) {
  d <- utils::read.delim(path, header = FALSE, sep = "\t",
                         stringsAsFactors = FALSE, comment.char = "#")
  if (ncol(d) < 3) om_abort("{.path {path}} does not look like a BED file.")

  out <- data.frame(
    chrom = as.character(d[[1]]),
    # BED starts are 0-based half-open; the tag's base is start + 1.
    pos = as.numeric(d[[2]]) + 1,
    stringsAsFactors = FALSE
  )
  out$count <- if (!is.null(count_col) && ncol(d) >= count_col) {
    v <- suppressWarnings(as.numeric(d[[count_col]]))
    ifelse(is.na(v), 1, v)
  } else {
    1
  }
  out$strand <- if (ncol(d) >= 6) as.character(d[[6]]) else "*"

  if (!is.null(strand)) out <- out[out$strand == strand, , drop = FALSE]
  if (!is.null(region)) {
    r <- parse_region(region)
    out <- out[out$chrom == r$chrom & out$pos >= r$start & out$pos <= r$end, ,
               drop = FALSE]
  }
  rownames(out) <- NULL
  out
}

#' Build a sashimi data object from 5'-tag data
#'
#' Turns single-base 5'-end tags (CAGE, STRT, CamoTSS) into the binned tracks a
#' sashimi panel draws, and optionally into arcs that run from each start site
#' to a shared downstream anchor.
#'
#' @details
#' For each sample the tags on the requested strand inside the window are
#' binned at `bin` base pairs and scaled to tags per million,
#' \eqn{t_i \times 10^6 / N_j}, using that library's total tag count. Each tag
#' is then extended `footprint` base pairs toward the 3' end, which gives the
#' track the body of a coverage profile without moving the start site. Finally
#' the samples of a group are averaged.
#'
#' Arcs in this mode are not junctions - a 5'-tag protocol produces no spliced
#' reads to count. An arc here is a pointer from a start site to the body of
#' the transcript it starts, and its label is that site's activity. Supply
#' those numbers through `activity`; the resulting figure is the one the
#' `omakase` sashimi style was designed for.
#'
#' @param tags A manifest read by [read_manifest()], a path to one, a vector of
#'   BED paths, or a data frame of tags with `chrom`, `pos`, `count`,
#'   `strand`.
#' @param region A region string, or a data frame of regions with `chrom`,
#'   `start`, `end` and optionally `name` and `strand`.
#' @param annotation Optional GTF/GFF3 for transcript models.
#' @param bin Bin width in base pairs.
#' @param footprint Distance each tag is extended toward the 3' end.
#' @param normalize `"tpm"` scales each library to tags per million; `"none"`
#'   leaves raw counts.
#' @param library_sizes Named vector of total tag counts per sample. Computed
#'   from the file when absent, which counts only tags inside the window unless
#'   the file is region-restricted already.
#' @param aggregate How to combine the samples of a group: `"mean"`,
#'   `"median"`, `"sum"`, or `"none"` to keep one track per sample.
#' @param strand Keep only tags on this strand. `"gene"` uses each region's own
#'   strand, which is usually what you want.
#' @param activity A data frame of start-site activities to draw as arcs, with
#'   columns `locus_id` (or `name`), `group`, `x0` (the start site), `x1` (the
#'   anchor), `count`, and optionally `role`.
#' @param group_col,label_col Manifest columns for grouping and labels.
#' @param gene Restrict the annotation to this gene.
#'
#' @return A `sashimi_data` object.
#'
#' @examples
#' tags <- data.frame(chrom = "chr1", pos = c(1200, 1200, 1800),
#'                    count = 1, strand = "+", sample = "s1", group = "early")
#' sashimi_from_tags(tags, "chr1:1000-2000", bin = 25, footprint = 100)
#'
#' @export
sashimi_from_tags <- function(tags, region, annotation = NULL, bin = 25,
                              footprint = 250, normalize = c("tpm", "none"),
                              library_sizes = NULL,
                              aggregate = c("mean", "median", "sum", "none"),
                              strand = "gene", activity = NULL,
                              group_col = NULL, label_col = NULL, gene = NULL) {
  normalize <- match.arg(normalize)
  aggregate <- match.arg(aggregate)

  man <- NULL
  tag_df <- NULL
  if (is.data.frame(tags)) {
    tag_df <- as_df(tags)
    require_cols(tag_df, c("chrom", "pos"), "tag table")
    tag_df$count <- as.numeric(col(tag_df, "count") %||% 1)
    tag_df$strand <- as.character(col(tag_df, "strand") %||% "*")
    tag_df$sample <- as.character(col(tag_df, "sample") %||% "all")
    tag_df$group <- as.character(col(tag_df, "group") %||% tag_df$sample)
  } else {
    man <- read_manifest(tags, group_col, label_col)
  }

  regions <- as_region_table(region)
  ann <- resolve_annotation(annotation, gene)

  loci <- list(); tracks <- list(); models <- list()

  for (ri in seq_len(nrow(regions))) {
    r <- new_region(regions$chrom[ri], regions$start[ri], regions$end[ri],
                    regions$strand[ri] %||% "*")
    lid <- regions$locus_id[ri]
    hint <- annotation_hint(ann, r)
    named <- !identical(regions$name[ri], format(r))
    gstrand <- if (!identical(r$strand, "*")) r$strand else hint$strand %||% "+"

    loci[[length(loci) + 1]] <- data.frame(
      locus_id = lid,
      gene_name = if (named) regions$name[ri] else hint$gene %||% regions$name[ri],
      chrom = r$chrom, strand = gstrand,
      win_lo = r$start, win_hi = r$end, stringsAsFactors = FALSE
    )

    want_strand <- if (identical(strand, "gene")) gstrand else strand
    if (identical(want_strand, "*")) want_strand <- NULL
    # An explicit strand also settles which way the footprint extends. Without
    # this, a window with no annotation falls back to "+" and a minus-strand
    # tag set would be smeared upstream of its own start sites.
    if (!identical(strand, "gene") && strand %in% c("+", "-")) gstrand <- strand

    per_sample <- if (!is.null(tag_df)) {
      split(tag_df, tag_df$sample)
    } else {
      stats::setNames(lapply(seq_len(nrow(man)), function(si) {
        d <- read_tag_bed(man$path[si], region = r, strand = want_strand)
        d$sample <- man$sample[si]
        d$group <- man$group[si]
        d
      }), man$sample)
    }

    for (nm in names(per_sample)) {
      d <- per_sample[[nm]]
      if (!is.null(want_strand) && has_col(d, "strand")) {
        d <- d[d$strand %in% c(want_strand, "*"), , drop = FALSE]
      }
      d <- d[d$chrom == r$chrom & d$pos >= r$start & d$pos <= r$end, ,
             drop = FALSE]
      N <- library_sizes[[nm]] %||% sum(d$count)
      tr <- bin_tags(d, r, bin = bin, footprint = footprint,
                     strand = gstrand,
                     scale = if (identical(normalize, "tpm") && N > 0) 1e6 / N else 1)
      if (is.null(tr)) next
      tr$locus_id <- lid
      tr$sample <- nm
      tr$group <- (col(d, "group") %||% nm)[1]
      tracks[[length(tracks) + 1]] <- tr
    }

    if (!is.null(ann)) {
      m <- models_in_region(ann, r)
      if (nrow(m)) {
        m$locus_id <- lid
        models[[length(models) + 1]] <- m
      }
    }
  }

  x <- sashimi_data(
    loci = rbind_all(loci), tracks = rbind_all(tracks),
    models = rbind_all(models),
    junctions = resolve_activity(activity, rbind_all(loci)),
    meta = list(source = "tags", bin = bin, footprint = footprint,
                normalize = normalize,
                group_order = if (is.null(man)) NULL else unique(man$group))
  )
  if (!identical(aggregate, "none")) {
    # Junctions here are supplied activities, not counted reads, so they must
    # not be re-aggregated along with the tracks.
    keep <- x$junctions
    x <- aggregate_tracks(x, aggregate)
    x$junctions <- keep
  }
  x
}

# Bin tags, then extend each bin's signal `footprint` bp downstream. The
# extension is a running sum over the bins ahead (or behind, on the minus
# strand), which is the cheap equivalent of widening every tag individually.
#' @noRd
bin_tags <- function(d, r, bin, footprint, strand, scale = 1) {
  width <- region_width(r)
  nb <- ceiling(width / bin)
  if (nb < 1) return(NULL)

  v <- numeric(nb)
  if (nrow(d)) {
    idx <- pmin(nb, pmax(1, floor((d$pos - r$start) / bin) + 1))
    agg <- tapply(d$count, idx, sum)
    v[as.integer(names(agg))] <- as.numeric(agg)
  }
  v <- v * scale

  fp_bins <- max(0, floor(footprint / bin))
  if (fp_bins > 0) {
    v <- extend_3prime(v, strand, fp_bins)
  }
  data.frame(
    bin = seq_len(nb) - 1L,
    pos = r$start + (seq_len(nb) - 1L) * bin,
    value = v,
    strand = strand,
    stringsAsFactors = FALSE
  )
}

# Spread each bin's signal over the `n` bins in the 3' direction, so a
# single-base tag becomes a plateau that starts at the tag.
#' @noRd
extend_3prime <- function(v, strand, n) {
  if (n <= 0) return(v)
  k <- length(v)
  out <- numeric(k)
  if (identical(strand, "-")) {
    # Transcription runs toward lower coordinates, so extend leftward.
    cs <- cumsum(c(0, v))
    for (i in seq_len(k)) {
      hi <- min(k, i + n)
      out[i] <- cs[hi + 1] - cs[i]
    }
  } else {
    cs <- cumsum(c(0, v))
    for (i in seq_len(k)) {
      lo <- max(1, i - n)
      out[i] <- cs[i + 1] - cs[lo]
    }
  }
  out
}

# Activities supplied for the arcs: accept a data frame keyed by locus_id or by
# gene name.
#' @noRd
resolve_activity <- function(activity, loci) {
  if (is.null(activity)) return(NULL)
  a <- as_df(activity)
  if (!has_col(a, "locus_id") && has_col(a, "name") && !is.null(loci)) {
    a$locus_id <- loci$locus_id[match(a$name, loci$gene_name)]
  }
  if (!has_col(a, "locus_id") && has_col(a, "gene_name") && !is.null(loci)) {
    a$locus_id <- loci$locus_id[match(a$gene_name, loci$gene_name)]
  }
  require_cols(a, c("locus_id", "group", "x0", "x1", "count"), "activity table")
  a
}
