# ---------------------------------------------------------------------------
# Reading coverage and junctions out of alignment files.
#
# Junctions are the N operations in the CIGAR string; coverage is everything
# else that consumes reference. Both are computed in one pass over the reads,
# because a second pass over a large BAM is the slowest thing this package
# does.
# ---------------------------------------------------------------------------

#' Strand specificity settings
#'
#' How the strand of a read is inferred from the flag, for stranded libraries.
#'
#' @return A character vector of setting names.
#' @examples
#' strand_modes()
#' @export
strand_modes <- function() {
  c("none", "sense", "antisense", "mate1_sense", "mate2_sense")
}

# Which reads get their apparent strand flipped, following the convention
# ggsashimi uses. For an unstranded library everything is called "+" and the
# two strands are never separated.
#' @noRd
flip_strand <- function(mode, flag) {
  is_paired <- bitwAnd(flag, 1L) > 0L
  is_mate1 <- bitwAnd(flag, 64L) > 0L
  switch(mode,
    none = rep(FALSE, length(flag)),
    sense = rep(FALSE, length(flag)),
    antisense = rep(TRUE, length(flag)),
    # For a mate1-sense protocol, mate 2 reports the opposite strand.
    mate1_sense = ifelse(is_paired, !is_mate1, FALSE),
    mate2_sense = ifelse(is_paired, is_mate1, TRUE)
  )
}

#' Read coverage and junctions from alignment files
#'
#' Reads one or more BAM/CRAM files over a region and returns a
#' [sashimi_data()] object holding binned coverage and spliced-junction counts.
#' This is the main entry point when starting from alignments.
#'
#' @details
#' Coverage is the depth of reference-consuming CIGAR operations (`M`, `D`,
#' `=`, `X`), binned to `bin` base pairs by taking the mean depth within each
#' bin. Junctions are `N` operations, counted once per read and keyed by their
#' donor and acceptor coordinates; a junction supported by fewer than
#' `min_count` reads is dropped.
#'
#' For a stranded library, set `strand` so that reads are assigned to the
#' correct strand, and `keep_strand` to draw only one of them. For an
#' unstranded library leave `strand = "none"` and everything is pooled.
#'
#' @param bam A path to a BAM/CRAM/SAM file, a vector of paths, or a manifest
#'   read by [read_manifest()] (a path to a TSV works directly). A `.sam` file
#'   is converted to a sorted, indexed BAM in the session's temporary directory
#'   the first time it is read, which is what `rmats2sashimiplot` does for its
#'   `--s1`/`--s2` inputs.
#' @param region A region string such as `"chr10:27035000-27050000"`, a
#'   [parse_region()] object, or a data frame of regions with columns `chrom`,
#'   `start`, `end` and optionally `name`.
#' @param annotation Optional GTF/GFF3 path, or the result of
#'   [read_annotation()], used to draw transcript models under the tracks.
#' @param bin Bin width in base pairs. `1` gives per-base coverage, which is
#'   exact but slow to draw over a wide window; the default of `NULL` picks a
#'   bin that gives roughly `n_bins` points.
#' @param n_bins Target number of bins when `bin` is `NULL`.
#' @param group_col,label_col Manifest columns used for grouping and panel
#'   labels; see [read_manifest()].
#' @param strand Strand specificity, see [strand_modes()].
#' @param keep_strand Which strand to keep once reads are assigned: `"both"`,
#'   `"+"` or `"-"`. Only meaningful when `strand` is not `"none"`.
#' @param min_count Drop junctions with fewer supporting reads.
#' @param junction_overlap Which junctions to keep: `"within"` (the default)
#'   keeps only junctions whose donor and acceptor both fall inside the region,
#'   `"any"` also keeps those with one end beyond it, whose arcs then run off
#'   the edge of the panel.
#' @param min_mapq Ignore alignments below this mapping quality.
#' @param flags A [Rsamtools::scanBamFlag()] value. The default drops
#'   unmapped reads, secondary alignments, duplicates and QC failures.
#' @param per_sample Keep one track per sample. When `FALSE`, samples sharing a
#'   group are summed as they are read, which is cheaper in memory.
#' @param gene Restrict `annotation` to this gene name.
#' @param flank Widen the region by this many base pairs on each side.
#'
#' @return A `sashimi_data` object.
#'
#' @examples
#' bams <- system.file("extdata", "samples.tsv", package = "omakase")
#' gtf <- system.file("extdata", "annotation.gtf", package = "omakase")
#' if (nzchar(bams)) {
#'   sd <- sashimi_from_bam(bams, "chr10:27040584-27048100", annotation = gtf,
#'                          min_count = 10)
#'   plot_sashimi(sd, aggregate = "mean")
#' }
#'
#' @export
sashimi_from_bam <- function(bam, region, annotation = NULL, bin = NULL,
                             n_bins = 800, group_col = NULL, label_col = NULL,
                             strand = "none", keep_strand = "both",
                             min_count = 1, min_mapq = 0, flags = NULL,
                             per_sample = TRUE, gene = NULL, flank = 0,
                             junction_overlap = c("within", "any")) {
  need_pkg("Rsamtools", "reading alignment files")
  need_pkg("GenomicAlignments", "reading alignment files")
  strand <- rlang::arg_match(strand, strand_modes())
  junction_overlap <- match.arg(junction_overlap)

  man <- if (is.data.frame(bam)) bam else read_manifest(bam, group_col, label_col)
  missing <- man$path[!file.exists(man$path)]
  if (length(missing)) om_abort("Alignment file{?s} not found: {.path {missing}}.")
  man$path <- vapply(man$path, ensure_bam, character(1), USE.NAMES = FALSE)

  regions <- as_region_table(region)
  ann <- resolve_annotation(annotation, gene)

  flags <- flags %||% Rsamtools::scanBamFlag(
    isUnmappedQuery = FALSE, isSecondaryAlignment = FALSE,
    isDuplicate = FALSE, isNotPassingQualityControls = FALSE
  )

  loci <- vector("list", nrow(regions))
  tracks <- list(); juncs <- list(); models <- list()

  for (ri in seq_len(nrow(regions))) {
    r <- expand_region(
      new_region(regions$chrom[ri], regions$start[ri], regions$end[ri],
                 regions$strand[ri] %||% "*"),
      flank
    )
    lid <- regions$locus_id[ri]
    width <- region_width(r)
    bw <- bin %||% max(1, round(width / n_bins))

    # When an annotation is supplied it knows the gene's name and strand better
    # than a coordinate string does, so let it fill both in unless the caller
    # was explicit.
    named <- !is.na(regions$name[ri]) && !identical(regions$name[ri], format(r))
    hint <- annotation_hint(ann, r)
    loci[[ri]] <- data.frame(
      locus_id = lid,
      gene_name = if (named) regions$name[ri] else hint$gene %||% regions$name[ri],
      chrom = r$chrom,
      strand = if (!identical(r$strand, "*")) r$strand else hint$strand %||% "*",
      win_lo = r$start, win_hi = r$end,
      stringsAsFactors = FALSE
    )

    for (si in seq_len(nrow(man))) {
      one <- read_one_bam(man$path[si], r, strand = strand,
                          keep_strand = keep_strand, min_mapq = min_mapq,
                          flags = flags, bin = bw,
                          junction_overlap = junction_overlap)
      if (!is.null(one$cov) && nrow(one$cov)) {
        one$cov$locus_id <- lid
        one$cov$group <- man$group[si]
        one$cov$sample <- if (per_sample) man$sample[si] else NA_character_
        tracks[[length(tracks) + 1]] <- one$cov
      }
      if (!is.null(one$jun) && nrow(one$jun)) {
        one$jun$locus_id <- lid
        one$jun$group <- man$group[si]
        one$jun$sample <- if (per_sample) man$sample[si] else NA_character_
        juncs[[length(juncs) + 1]] <- one$jun
      }
    }

    if (!is.null(ann)) {
      m <- models_in_region(ann, r)
      if (nrow(m)) {
        m$locus_id <- lid
        models[[length(models) + 1]] <- m
      }
    }
  }

  tracks <- rbind_all(tracks)
  juncs <- rbind_all(juncs)

  if (!per_sample && !is.null(tracks)) {
    tracks <- stats::aggregate(
      tracks$value,
      by = list(locus_id = tracks$locus_id, group = tracks$group,
                pos = tracks$pos, strand = tracks$strand),
      FUN = sum
    )
    names(tracks)[ncol(tracks)] <- "value"
  }
  if (!is.null(juncs)) {
    juncs <- stats::aggregate(
      juncs$count,
      by = list(locus_id = juncs$locus_id, group = juncs$group,
                sample = ifelse(is.na(juncs$sample), "", juncs$sample),
                x0 = juncs$x0, x1 = juncs$x1, strand = juncs$strand),
      FUN = sum
    )
    names(juncs)[ncol(juncs)] <- "count"
    juncs$sample[juncs$sample == ""] <- NA_character_
    if (min_count > 0) juncs <- juncs[juncs$count >= min_count, , drop = FALSE]
  }

  sashimi_data(
    loci = rbind_all(loci), tracks = tracks, junctions = juncs,
    models = rbind_all(models),
    meta = list(
      source = "bam", n_files = nrow(man), strand = strand,
      junction_overlap = junction_overlap,
      keep_strand = keep_strand, min_count = min_count, min_mapq = min_mapq,
      group_order = unique(man$group)
    )
  )
}

# One file, one region. Returns binned coverage and a junction table.
# Rsamtools reads BAM and CRAM but not SAM, so a SAM input is converted once
# and cached for the session. Sorting and indexing are required for the
# region query that follows.
#' @noRd
ensure_bam <- function(path) {
  if (!grepl("\\.sam$", path, ignore.case = TRUE)) return(path)
  need_pkg("Rsamtools", "converting SAM to BAM")
  key <- paste0("omakase_sam_", tools::file_path_sans_ext(basename(path)))
  dest <- file.path(tempdir(), key)
  if (file.exists(paste0(dest, ".bam"))) return(paste0(dest, ".bam"))
  om_inform("Converting {.path {basename(path)}} to an indexed BAM (once per session).")
  raw <- Rsamtools::asBam(path, destination = paste0(dest, "_unsorted"),
                          overwrite = TRUE, indexDestination = FALSE)
  sorted <- Rsamtools::sortBam(raw, destination = dest)
  Rsamtools::indexBam(sorted)
  sorted
}

#' @noRd
read_one_bam <- function(path, r, strand, keep_strand, min_mapq, flags, bin,
                         junction_overlap = "within") {
  which <- region_to_granges(new_region(r$chrom, r$start, r$end, "*"))
  param <- Rsamtools::ScanBamParam(which = which, flag = flags,
                                   what = c("flag", "mapq"))

  ga <- tryCatch(
    GenomicAlignments::readGAlignments(path, param = param, use.names = FALSE),
    error = function(e) {
      om_abort(c("Could not read {.path {path}}.",
                 "x" = conditionMessage(e),
                 "i" = "Is the file indexed? Try {.code Rsamtools::indexBam()}."))
    }
  )
  if (!length(ga)) return(list(cov = NULL, jun = NULL))

  meta <- S4Vectors::mcols(ga)
  keep <- rep(TRUE, length(ga))
  if (min_mapq > 0 && has_col(meta, "mapq")) {
    keep <- keep & !is.na(meta$mapq) & meta$mapq >= min_mapq
  }

  # Strand assignment. The aligner's own strand is flipped for the protocols
  # where the read reports the opposite of the transcript it came from.
  rs <- rep("+", length(ga))
  if (!identical(strand, "none")) {
    flag <- col(meta, "flag") %||% rep(0L, length(ga))
    rev <- as.character(GenomicAlignments::strand(ga)) == "-"
    rs <- ifelse(xor(flip_strand(strand, flag), rev), "-", "+")
    if (!identical(keep_strand, "both")) keep <- keep & rs == keep_strand
  }
  ga <- ga[keep]
  rs <- rs[keep]
  if (!length(ga)) return(list(cov = NULL, jun = NULL))

  strands <- if (identical(strand, "none")) "*" else unique(rs)
  cov_parts <- list(); jun_parts <- list()

  for (s in strands) {
    sel <- if (identical(s, "*")) rep(TRUE, length(ga)) else rs == s
    g <- ga[sel]
    if (!length(g)) next

    cov <- binned_coverage(g, r, bin)
    if (!is.null(cov)) {
      cov$strand <- s
      cov_parts[[length(cov_parts) + 1]] <- cov
    }
    jun <- junction_counts(g, r, overlap = junction_overlap)
    if (!is.null(jun)) {
      jun$strand <- s
      jun_parts[[length(jun_parts) + 1]] <- jun
    }
  }
  list(cov = rbind_all(cov_parts), jun = rbind_all(jun_parts))
}

# Mean depth per bin over the window. grglist() expands each alignment into its
# aligned blocks, so N gaps correctly contribute no coverage.
#' @noRd
binned_coverage <- function(ga, r, bin) {
  blocks <- unlist(GenomicAlignments::grglist(ga), use.names = FALSE)
  blocks <- blocks[GenomicRanges::seqnames(blocks) == r$chrom]
  if (!length(blocks)) return(NULL)

  cvg <- GenomicRanges::coverage(blocks)[[r$chrom]]
  win <- IRanges::IRanges(start = r$start, end = min(r$end, length(cvg)))
  if (IRanges::width(win) <= 0) return(NULL)
  v <- as.numeric(cvg[win])
  # Reads may stop short of the window's right edge; pad so every sample has
  # the same bin grid and panels line up.
  full <- region_width(r)
  if (length(v) < full) v <- c(v, rep(0, full - length(v)))

  starts <- seq(1, full, by = bin)
  idx <- rep(seq_along(starts), each = bin, length.out = full)
  value <- as.numeric(tapply(v, idx, mean))
  data.frame(
    bin = seq_along(starts) - 1L,
    pos = r$start + starts - 1L,
    value = value,
    stringsAsFactors = FALSE
  )
}

# Spliced junctions: the N operations. junctions() gives the intron ranges
# directly, so no CIGAR parsing is needed here.
#' @noRd
junction_counts <- function(ga, r, overlap = "within") {
  jx <- GenomicAlignments::junctions(ga, use.mcols = FALSE)
  jx <- unlist(jx, use.names = FALSE)
  if (!length(jx)) return(NULL)
  jx <- jx[GenomicRanges::seqnames(jx) == r$chrom]
  if (!length(jx)) return(NULL)

  st <- GenomicRanges::start(jx)
  en <- GenomicRanges::end(jx)
  # An arc with one end outside the window has its apex off-screen, so what
  # shows is a near-flat line entering from the edge - visually it reads as a
  # junction that is not there. Keeping only fully-contained junctions is both
  # clearer and what ggsashimi does; "any" restores the partial arcs for
  # callers who want to see that signal leaves the window.
  keep <- if (identical(overlap, "any")) {
    st <= r$end & en >= r$start
  } else {
    st >= r$start & en <= r$end
  }
  st <- st[keep]; en <- en[keep]
  if (!length(st)) return(NULL)

  tab <- table(paste(st, en, sep = "-"))
  parts <- do.call(rbind, strsplit(names(tab), "-", fixed = TRUE))
  data.frame(
    x0 = as.numeric(parts[, 1]) - 1,
    x1 = as.numeric(parts[, 2]) + 1,
    count = as.numeric(tab),
    stringsAsFactors = FALSE
  )
}

# Normalise the `region` argument into a table, so one code path handles a
# single region string and a batch of them.
#' @noRd
as_region_table <- function(region) {
  if (is.data.frame(region)) {
    d <- as_df(region)
    require_cols(d, c("chrom", "start", "end"), "region table")
    d$name <- as.character(col(d, "name") %||% col(d, "locus_id") %||%
                             sprintf("%s:%s-%s", d$chrom, d$start, d$end))
    d$locus_id <- as.character(col(d, "locus_id") %||% d$name)
    d$strand <- as.character(col(d, "strand") %||% "*")
    return(d)
  }
  if (is.character(region) && length(region) > 1L) {
    return(rbind_all(lapply(region, as_region_table)))
  }
  r <- parse_region(region)
  if (is.na(r$start)) {
    om_abort("Region {.val {r$chrom}} has no coordinates; give {.val chr:start-end}.")
  }
  data.frame(
    locus_id = format(r), name = format(r), chrom = r$chrom,
    start = r$start, end = r$end, strand = r$strand,
    stringsAsFactors = FALSE
  )
}
