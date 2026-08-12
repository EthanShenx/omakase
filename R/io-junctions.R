# ---------------------------------------------------------------------------
# Junction files.
#
# Not every project keeps its BAMs. STAR's SJ.out.tab and regtools' junction
# BED are both small, both survive when the alignments do not, and both carry
# everything an arc needs.
# ---------------------------------------------------------------------------

#' Read a splice junction file
#'
#' Reads STAR `SJ.out.tab`, a regtools/TopHat junction BED, or a plain BED-like
#' table of junctions. The format is detected from the columns.
#'
#' @details
#' STAR's `SJ.out.tab` has no header and nine columns: contig, first intron
#' base, last intron base, strand code (`0` undefined, `1` `+`, `2` `-`),
#' intron motif, annotated flag, uniquely-mapping reads, multi-mapping reads,
#' maximum overhang. The uniquely-mapping count is used unless
#' `include_multimappers` is set.
#'
#' A regtools junction BED has twelve columns, where the thick start/end and
#' block sizes describe the anchors flanking the intron; the intron itself is
#' recovered from them.
#'
#' @param path Path to the junction file.
#' @param format `"auto"`, `"star"`, `"bed"`, or `"plain"`.
#' @param region Optional region to restrict to.
#' @param min_count Drop junctions with fewer supporting reads.
#' @param include_multimappers For STAR input, add the multi-mapping read count
#'   to the unique count.
#' @param annotated_only For STAR input, keep only junctions the annotation
#'   already contains.
#'
#' @return A data frame with `chrom`, `x0`, `x1`, `count`, `strand`.
#'
#' @examples
#' f <- tempfile()
#' writeLines("chr1\t100\t200\t1\t1\t1\t42\t3\t30", f)
#' read_junctions(f)
#'
#' @export
read_junctions <- function(path, format = c("auto", "star", "bed", "plain"),
                           region = NULL, min_count = 1,
                           include_multimappers = FALSE,
                           annotated_only = FALSE) {
  format <- match.arg(format)
  raw <- utils::read.delim(path, header = FALSE, sep = "\t",
                           stringsAsFactors = FALSE, comment.char = "#")
  # A header row would leave the coordinate columns non-numeric.
  if (ncol(raw) >= 3 && anyNA(suppressWarnings(as.numeric(raw[[2]][1])))) {
    raw <- utils::read.delim(path, header = TRUE, sep = "\t",
                             stringsAsFactors = FALSE, comment.char = "#")
  }
  if (identical(format, "auto")) {
    format <- if (ncol(raw) >= 9 && all(raw[[4]] %in% c(0, 1, 2))) {
      "star"
    } else if (ncol(raw) >= 12) {
      "bed"
    } else {
      "plain"
    }
  }

  out <- switch(format,
    star = {
      cnt <- as.numeric(raw[[7]])
      if (include_multimappers) cnt <- cnt + as.numeric(raw[[8]])
      d <- data.frame(
        chrom = as.character(raw[[1]]),
        x0 = as.numeric(raw[[2]]) - 1,
        x1 = as.numeric(raw[[3]]) + 1,
        count = cnt,
        strand = c("*", "+", "-")[as.numeric(raw[[4]]) + 1],
        annotated = as.numeric(raw[[6]]),
        stringsAsFactors = FALSE
      )
      if (annotated_only) d <- d[d$annotated > 0, , drop = FALSE]
      d[, c("chrom", "x0", "x1", "count", "strand")]
    },
    bed = {
      # Blocks are the anchors either side of the intron; the intron runs from
      # the end of the first block to the start of the second.
      sizes <- lapply(strsplit(as.character(raw[[11]]), ","), as.numeric)
      starts <- lapply(strsplit(as.character(raw[[12]]), ","), as.numeric)
      bstart <- as.numeric(raw[[2]])
      x0 <- bstart + vapply(sizes, function(s) s[1], numeric(1))
      x1 <- bstart + vapply(starts, function(s) s[2], numeric(1)) + 1
      data.frame(
        chrom = as.character(raw[[1]]), x0 = x0, x1 = x1,
        count = as.numeric(raw[[5]]),
        strand = as.character(raw[[6]]),
        stringsAsFactors = FALSE
      )
    },
    plain = {
      d <- as_df(raw)
      if (!has_col(d, "chrom")) {
        names(d)[seq_len(min(4, ncol(d)))] <-
          c("chrom", "x0", "x1", "count")[seq_len(min(4, ncol(d)))]
      }
      d <- rename_cols(d, list(chrom = "chr", x0 = "start", x1 = "end",
                               count = "score"))
      require_cols(d, c("chrom", "x0", "x1"), "junction table")
      d$count <- as.numeric(col(d, "count") %||% 1)
      d$strand <- as.character(col(d, "strand") %||% "*")
      d[, c("chrom", "x0", "x1", "count", "strand")]
    }
  )
  out$x0 <- as.numeric(out$x0)
  out$x1 <- as.numeric(out$x1)
  out <- out[out$count >= min_count, , drop = FALSE]

  if (!is.null(region)) {
    r <- parse_region(region)
    out <- out[out$chrom == r$chrom & out$x0 >= r$start & out$x1 <= r$end, ,
               drop = FALSE]
  }
  rownames(out) <- NULL
  out
}

#' Build a sashimi data object from junction files
#'
#' For projects that keep junction tables rather than alignments. Produces
#' arcs but no coverage, which draws a perfectly readable figure - the arcs
#' carry the splicing signal, and the filled area is context.
#'
#' @param junctions A path, a vector of paths, a manifest read by
#'   [read_manifest()], or a data frame of junctions.
#' @param region A region string, or a data frame of regions.
#' @param annotation Optional GTF/GFF3 for transcript models.
#' @param group_col,label_col Manifest columns for grouping and labels.
#' @param min_count Drop junctions with fewer supporting reads.
#' @param format Junction file format, see [read_junctions()].
#' @param gene Restrict the annotation to this gene.
#' @param ... Passed to [read_junctions()].
#'
#' @return A `sashimi_data` object with an empty `tracks` slot.
#'
#' @examples
#' f <- tempfile()
#' writeLines(c("chr1\t1200\t1800\t1\t1\t1\t42\t0\t30"), f)
#' sashimi_from_junctions(f, "chr1:1000-2000")
#'
#' @export
sashimi_from_junctions <- function(junctions, region, annotation = NULL,
                                   group_col = NULL, label_col = NULL,
                                   min_count = 1, format = "auto",
                                   gene = NULL, ...) {
  man <- if (is.data.frame(junctions) && has_col(junctions, "path")) {
    junctions
  } else if (is.data.frame(junctions)) {
    NULL
  } else {
    m <- read_manifest(junctions, group_col, label_col)
    m
  }
  regions <- as_region_table(region)
  ann <- resolve_annotation(annotation, gene)

  loci <- list(); juncs <- list(); models <- list()
  for (ri in seq_len(nrow(regions))) {
    r <- new_region(regions$chrom[ri], regions$start[ri], regions$end[ri],
                    regions$strand[ri] %||% "*")
    lid <- regions$locus_id[ri]
    hint <- annotation_hint(ann, r)
    # `format` is this function's own argument, so the region formatter has to
    # be reached explicitly.
    named <- !identical(regions$name[ri], format.omakase_region(r))
    loci[[length(loci) + 1]] <- data.frame(
      locus_id = lid,
      gene_name = if (named) regions$name[ri] else hint$gene %||% regions$name[ri],
      chrom = r$chrom,
      strand = if (!identical(r$strand, "*")) r$strand else hint$strand %||% "*",
      win_lo = r$start, win_hi = r$end, stringsAsFactors = FALSE
    )

    if (is.null(man)) {
      d <- as_df(junctions)
      d <- d[d$chrom == r$chrom & d$x0 >= r$start & d$x1 <= r$end, , drop = FALSE]
      if (nrow(d)) {
        d$locus_id <- lid
        d$group <- col(d, "group") %||% "all"
        juncs[[length(juncs) + 1]] <- d
      }
    } else {
      for (si in seq_len(nrow(man))) {
        d <- read_junctions(man$path[si], format = format, region = r,
                            min_count = min_count, ...)
        if (!nrow(d)) next
        d$locus_id <- lid
        d$group <- man$group[si]
        d$sample <- man$sample[si]
        juncs[[length(juncs) + 1]] <- d
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

  sashimi_data(
    loci = rbind_all(loci), junctions = rbind_all(juncs),
    models = rbind_all(models),
    meta = list(source = "junctions", min_count = min_count,
                group_order = if (is.null(man)) NULL else unique(man$group))
  )
}
