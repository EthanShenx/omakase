# ---------------------------------------------------------------------------
# BED12.
#
# The format genome browsers use for transcript models: twelve columns, where
# the blocks are exons and the thick range is the coding sequence. Splitting
# each exon against the thick range is what lets a track draw coding and
# untranslated parts at different heights, the way UCSC and pyGenomeTracks do.
# ---------------------------------------------------------------------------

BED12_COLS <- c("chrom", "start", "end", "name", "score", "strand",
                "thick_start", "thick_end", "item_rgb", "block_count",
                "block_sizes", "block_starts")

#' Read a BED file
#'
#' Reads BED3 through BED12. For BED12 the blocks are expanded into exons and
#' each exon is split against the thick range, so the result carries a
#' `feature` column of `"CDS"` and `"UTR"` that a genome track can draw at
#' different heights.
#'
#' @param path Path to the BED file, optionally gzipped.
#' @param region Optional region to restrict to.
#' @param expand For BED12, expand blocks into one row per exon (the default).
#'   `FALSE` keeps one row per transcript.
#' @param name_col Column holding the transcript name; defaults to column 4.
#'
#' @return A data frame with `chrom`, `start`, `end`, `tx_id`, `strand`,
#'   `feature`, and `color` when the file carries an `itemRgb` column.
#'   Coordinates are 1-based inclusive.
#'
#' @examples
#' f <- tempfile()
#' writeLines(paste("chr1", 1000, 5000, "tx1", 0, "+", 1200, 4000,
#'                  "31,120,180", 2, "500,800,", "0,3200,", sep = "\t"), f)
#' read_bed(f)
#'
#' @export
read_bed <- function(path, region = NULL, expand = TRUE, name_col = 4) {
  if (!file.exists(path)) om_abort("BED file not found: {.path {path}}.")
  d <- utils::read.delim(path, header = FALSE, sep = "\t",
                         stringsAsFactors = FALSE, comment.char = "#")
  if (ncol(d) < 3) om_abort("{.path {path}} does not look like a BED file.")
  names(d)[seq_len(min(ncol(d), length(BED12_COLS)))] <-
    BED12_COLS[seq_len(min(ncol(d), length(BED12_COLS)))]

  # BED is 0-based half-open; the package is 1-based inclusive throughout.
  d$start <- as.numeric(d$start) + 1
  d$end <- as.numeric(d$end)
  d$chrom <- as.character(d$chrom)
  d$tx_id <- if (ncol(d) >= name_col) {
    as.character(d[[name_col]])
  } else {
    sprintf("%s:%s-%s", d$chrom, d$start, d$end)
  }
  d$strand <- if (has_col(d, "strand")) as.character(d$strand) else "*"
  d$color <- if (has_col(d, "item_rgb")) rgb_string(d$item_rgb) else NA_character_

  if (!is.null(region)) {
    r <- parse_region(region)
    if (!is.na(r$start)) {
      d <- d[d$chrom == r$chrom & d$start <= r$end & d$end >= r$start, ,
             drop = FALSE]
    }
  }
  if (!nrow(d)) return(empty_bed())

  has_blocks <- all(c("block_count", "block_sizes", "block_starts") %in% names(d))
  if (!expand || !has_blocks) {
    out <- d[, c("chrom", "start", "end", "tx_id", "strand", "color")]
    out$feature <- "exon"
    out$thick_start <- if (has_col(d, "thick_start")) as.numeric(d$thick_start) + 1 else NA_real_
    out$thick_end <- if (has_col(d, "thick_end")) as.numeric(d$thick_end) else NA_real_
    rownames(out) <- NULL
    return(out)
  }
  bed12_to_exons(d)
}

#' @rdname read_bed
#' @export
read_bed12 <- function(path, region = NULL, expand = TRUE, name_col = 4) {
  read_bed(path, region = region, expand = expand, name_col = name_col)
}

#' @noRd
empty_bed <- function() {
  data.frame(chrom = character(0), start = numeric(0), end = numeric(0),
             tx_id = character(0), strand = character(0),
             color = character(0), feature = character(0),
             thick_start = numeric(0), thick_end = numeric(0),
             stringsAsFactors = FALSE)
}

# "31,120,180" -> "#1F78B4". A single "0" means "no colour given", which is
# what most tools write when they do not use itemRgb.
#' @noRd
rgb_string <- function(x) {
  x <- as.character(x)
  vapply(x, function(v) {
    if (is.na(v) || !nzchar(v) || identical(v, "0")) return(NA_character_)
    p <- suppressWarnings(as.numeric(strsplit(v, ",", fixed = TRUE)[[1]]))
    if (length(p) != 3 || anyNA(p)) return(NA_character_)
    grDevices::rgb(p[1], p[2], p[3], maxColorValue = 255)
  }, character(1), USE.NAMES = FALSE)
}

# Expand BED12 blocks into exons, splitting each against the thick range so the
# coding and untranslated parts become separate rows.
#' @noRd
bed12_to_exons <- function(d) {
  parts <- lapply(seq_len(nrow(d)), function(i) {
    sizes <- as.numeric(strsplit(trimws(as.character(d$block_sizes[i])), ",")[[1]])
    offs <- as.numeric(strsplit(trimws(as.character(d$block_starts[i])), ",")[[1]])
    sizes <- sizes[!is.na(sizes)]
    offs <- offs[!is.na(offs)]
    if (!length(sizes) || length(sizes) != length(offs)) {
      return(data.frame(chrom = d$chrom[i], start = d$start[i], end = d$end[i],
                        tx_id = d$tx_id[i], strand = d$strand[i],
                        color = d$color[i], feature = "exon",
                        stringsAsFactors = FALSE))
    }
    ex_start <- d$start[i] + offs
    ex_end <- ex_start + sizes - 1

    ts <- as.numeric(d$thick_start[i]) + 1
    te <- as.numeric(d$thick_end[i])
    # A thick range of zero width means the transcript is non-coding, so every
    # block is untranslated.
    coding <- is.finite(ts) && is.finite(te) && te >= ts

    rows <- lapply(seq_along(ex_start), function(k) {
      s <- ex_start[k]; e <- ex_end[k]
      if (!coding) {
        return(data.frame(start = s, end = e, feature = "UTR"))
      }
      seg <- list()
      if (s < ts) seg[[length(seg) + 1]] <- data.frame(
        start = s, end = min(e, ts - 1), feature = "UTR")
      cs <- max(s, ts); ce <- min(e, te)
      if (ce >= cs) seg[[length(seg) + 1]] <- data.frame(
        start = cs, end = ce, feature = "CDS")
      if (e > te) seg[[length(seg) + 1]] <- data.frame(
        start = max(s, te + 1), end = e, feature = "UTR")
      rbind_all(seg)
    })
    out <- rbind_all(rows)
    if (is.null(out)) return(NULL)
    out$chrom <- d$chrom[i]
    out$tx_id <- d$tx_id[i]
    out$strand <- d$strand[i]
    out$color <- d$color[i]
    out[out$end >= out$start, , drop = FALSE]
  })
  out <- rbind_all(parts)
  if (is.null(out)) return(empty_bed())
  out <- out[, c("chrom", "start", "end", "tx_id", "strand", "color", "feature")]
  rownames(out) <- NULL
  out
}

#' Write junctions to a BED file
#'
#' Writes the `junctions` slot as a BED12 file in the layout regtools and
#' TopHat use: two anchor blocks either side of the intron, the score column
#' carrying the supporting read count. This is the counterpart of `ggsashimi`'s
#' `--junctions-bed`, and lets counted junctions be loaded into a genome
#' browser or fed back in through [read_junctions()].
#'
#' @param x A `sashimi_data` object, or a junction data frame with `x0`, `x1`
#'   and `count`.
#' @param file Output path.
#' @param anchor Width of the flanking blocks, in base pairs.
#' @param name_prefix Prefix for the junction names.
#' @return `file`, invisibly.
#' @examples
#' j <- data.frame(locus_id = "a", group = "g", x0 = 1000, x1 = 2000,
#'                 count = 42, chrom = "chr1")
#' f <- tempfile(fileext = ".bed")
#' write_junctions(j, f)
#' readLines(f)
#' @export
write_junctions <- function(x, file, anchor = 25, name_prefix = "JUNC") {
  j <- if (inherits(x, "sashimi_data")) {
    d <- x$junctions
    if (nrow(d)) {
      d$chrom <- x$loci$chrom[match(d$locus_id, x$loci$locus_id)]
      d$strand <- ifelse(is.na(d$strand) | d$strand == "*",
                         x$loci$strand[match(d$locus_id, x$loci$locus_id)],
                         d$strand)
    }
    d
  } else {
    as_df(x)
  }
  require_cols(j, c("x0", "x1", "count"), "junction table")
  if (!has_col(j, "chrom")) om_abort("Junctions need a {.field chrom} column to be written as BED.")
  # A plain junction table need not carry a strand; BED writes "." for unknown.
  if (!has_col(j, "strand")) j$strand <- NA_character_
  if (!nrow(j)) {
    writeLines(character(0), file)
    return(invisible(file))
  }

  # `x0` is the last exonic base before the intron and `x1` the first one
  # after, so the two anchor blocks must END at x0 and START at x1. In 0-based
  # half-open BED terms the record therefore runs from `x0 - anchor` to
  # `x1 + anchor - 1`, and the second block begins at `x1 - 1`.
  bstart <- j$x0 - anchor
  bend <- j$x1 + anchor - 1
  bed <- data.frame(
    chrom = j$chrom,
    start = bstart,
    end = bend,
    name = sprintf("%s%05d", name_prefix, seq_len(nrow(j))),
    score = round(j$count),
    strand = ifelse(is.na(j$strand), ".", j$strand),
    thick_start = bstart,
    thick_end = bend,
    item_rgb = "255,0,0",
    block_count = 2L,
    block_sizes = sprintf("%d,%d,", anchor, anchor),
    block_starts = sprintf("0,%d,", (j$x1 - 1) - bstart),
    stringsAsFactors = FALSE
  )
  utils::write.table(bed, file, sep = "\t", quote = FALSE, row.names = FALSE,
                     col.names = FALSE)
  invisible(file)
}
