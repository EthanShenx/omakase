# ---------------------------------------------------------------------------
# rMATS event tables.
#
# rMATS writes one file per event type, and each type lays its exon coordinates
# out in different columns. The tables below record that layout once so the
# rest of the code can treat every event type the same way.
# ---------------------------------------------------------------------------

#' Alternative splicing event types
#'
#' The five event classes rMATS reports.
#'
#' @return A character vector: `SE`, `A5SS`, `A3SS`, `MXE`, `RI`.
#' @examples
#' event_types()
#' @export
event_types <- function() c("SE", "A5SS", "A3SS", "MXE", "RI")

# For each event type: the columns holding each exon, and which of them is the
# alternative feature whose inclusion PSI is being measured.
RMATS_LAYOUT <- list(
  SE = list(
    exons = list(
      upstream = c("upstreamES", "upstreamEE"),
      alt = c("exonStart_0base", "exonEnd"),
      downstream = c("downstreamES", "downstreamEE")
    ),
    alt_exon = "alt", label = "skipped exon"
  ),
  MXE = list(
    exons = list(
      upstream = c("upstreamES", "upstreamEE"),
      alt = c("1stExonStart_0base", "1stExonEnd"),
      alt2 = c("2ndExonStart_0base", "2ndExonEnd"),
      downstream = c("downstreamES", "downstreamEE")
    ),
    alt_exon = "alt", label = "mutually exclusive exons"
  ),
  A5SS = list(
    exons = list(
      alt = c("longExonStart_0base", "longExonEnd"),
      short = c("shortES", "shortEE"),
      flanking = c("flankingES", "flankingEE")
    ),
    alt_exon = "alt", label = "alternative 5' splice site"
  ),
  A3SS = list(
    exons = list(
      alt = c("longExonStart_0base", "longExonEnd"),
      short = c("shortES", "shortEE"),
      flanking = c("flankingES", "flankingEE")
    ),
    alt_exon = "alt", label = "alternative 3' splice site"
  ),
  RI = list(
    exons = list(
      alt = c("riExonStart_0base", "riExonEnd"),
      upstream = c("upstreamES", "upstreamEE"),
      downstream = c("downstreamES", "downstreamEE")
    ),
    alt_exon = "alt", label = "retained intron"
  )
)

#' Read an rMATS event table
#'
#' Reads a `*.MATS.JC.txt` or `*.MATS.JCEC.txt` file and returns one row per
#' event with the coordinates needed to draw it, plus rMATS's own inclusion
#' levels and significance.
#'
#' @param path Path to the rMATS output file.
#' @param event_type One of [event_types()]. Inferred from the file name when
#'   `NULL`.
#' @param fdr Keep only events with `FDR` at or below this.
#' @param min_dpsi Keep only events whose absolute `IncLevelDifference` is at
#'   least this.
#' @param top Keep only the `top` most significant events.
#' @param genes Restrict to these gene names or IDs.
#' @param chrom_style Rewrite contig names: `"keep"`, `"ucsc"` (ensure a `chr`
#'   prefix) or `"ensembl"` (strip it), for when the event file and the BAM
#'   header disagree.
#' @param flank Padding added around the event when computing its window.
#'
#' @return A data frame with one row per event: `event_id`, `event_type`,
#'   `gene_id`, `gene_name`, `chrom`, `strand`, the exon coordinates, the
#'   window (`win_lo`, `win_hi`), and rMATS's `psi_1`, `psi_2`, `dpsi`,
#'   `pvalue`, `fdr`.
#'
#' @examples
#' f <- system.file("extdata", "SE.MATS.JC.txt", package = "omakase")
#' if (nzchar(f)) read_rmats(f, "SE", fdr = 0.05, top = 20)
#'
#' @export
read_rmats <- function(path, event_type = NULL, fdr = NULL, min_dpsi = NULL,
                       top = NULL, genes = NULL, chrom_style = "keep",
                       flank = 300) {
  d <- read_table_any(path)
  event_type <- event_type %||% infer_event_type(path, d)
  event_type <- rlang::arg_match(event_type, event_types())
  layout <- RMATS_LAYOUT[[event_type]]

  need <- unique(unlist(layout$exons))
  require_cols(d, c("chr", "strand", need), paste0(event_type, " event file"))

  gene_name <- col(d, "geneSymbol") %||% col(d, "GeneID") %||% NA_character_
  gene_name <- gsub('^"|"$', "", as.character(gene_name))

  out <- data.frame(
    event_id = paste0(event_type, ".", col(d, "ID") %||% seq_len(nrow(d))),
    event_type = event_type,
    gene_id = as.character(col(d, "GeneID") %||% NA_character_),
    gene_name = gene_name,
    chrom = harmonise_chrom(as.character(d$chr), chrom_style),
    strand = as.character(d$strand),
    stringsAsFactors = FALSE
  )
  # rMATS writes exon starts 0-based and ends 1-based; shifting the starts here
  # puts everything on the package's 1-based inclusive footing.
  for (nm in names(layout$exons)) {
    cols <- layout$exons[[nm]]
    out[[paste0(nm, "_start")]] <- as.numeric(d[[cols[1]]]) + 1
    out[[paste0(nm, "_end")]] <- as.numeric(d[[cols[2]]])
  }

  starts <- out[, grep("_start$", names(out)), drop = FALSE]
  ends <- out[, grep("_end$", names(out)), drop = FALSE]
  out$win_lo <- apply(starts, 1, min, na.rm = TRUE) - flank
  out$win_hi <- apply(ends, 1, max, na.rm = TRUE) + flank

  out$psi_1 <- mean_inc_level(col(d, "IncLevel1"))
  out$psi_2 <- mean_inc_level(col(d, "IncLevel2"))
  out$dpsi <- suppressWarnings(as.numeric(d$IncLevelDifference %||% NA))
  out$pvalue <- suppressWarnings(as.numeric(d$PValue %||% NA))
  out$fdr <- suppressWarnings(as.numeric(d$FDR %||% NA))
  out$inc_counts_1 <- as.character(d$IJC_SAMPLE_1 %||% NA)
  out$skip_counts_1 <- as.character(d$SJC_SAMPLE_1 %||% NA)
  out$inc_counts_2 <- as.character(d$IJC_SAMPLE_2 %||% NA)
  out$skip_counts_2 <- as.character(d$SJC_SAMPLE_2 %||% NA)
  out$inc_len <- suppressWarnings(as.numeric(d$IncFormLen %||% NA))
  out$skip_len <- suppressWarnings(as.numeric(d$SkipFormLen %||% NA))

  if (!is.null(genes)) {
    out <- out[out$gene_name %in% genes | out$gene_id %in% genes, , drop = FALSE]
  }
  if (!is.null(fdr)) out <- out[!is.na(out$fdr) & out$fdr <= fdr, , drop = FALSE]
  if (!is.null(min_dpsi)) {
    out <- out[!is.na(out$dpsi) & abs(out$dpsi) >= min_dpsi, , drop = FALSE]
  }
  if (!is.null(top) && nrow(out) > top) {
    ord <- order(out$fdr, -abs(out$dpsi), na.last = TRUE)
    out <- out[utils::head(ord, top), , drop = FALSE]
  }
  rownames(out) <- NULL
  out
}

#' @noRd
infer_event_type <- function(path, d) {
  bn <- basename(path)
  for (et in event_types()) {
    if (grepl(paste0("(^|[^A-Z0-9])", et, "([^A-Z0-9]|$)"), bn)) return(et)
  }
  # Fall back on the columns: each layout has at least one column no other
  # layout uses.
  for (et in event_types()) {
    need <- unique(unlist(RMATS_LAYOUT[[et]]$exons))
    if (all(need %in% names(d))) return(et)
  }
  om_abort(c("Could not tell which event type {.path {path}} holds.",
             "i" = "Pass {.arg event_type}, one of {.val {event_types()}}."))
}

# rMATS writes per-replicate inclusion levels as a comma-separated string.
#' @noRd
mean_inc_level <- function(x) {
  if (is.null(x)) return(NA_real_)
  vapply(strsplit(as.character(x), ",", fixed = TRUE), function(v) {
    v <- suppressWarnings(as.numeric(v))
    if (!length(v) || all(is.na(v))) NA_real_ else mean(v, na.rm = TRUE)
  }, numeric(1))
}

#' Build a sashimi data object from rMATS events and alignments
#'
#' Reads an rMATS event file, turns each event into a plotting window with the
#' event's exons as transcript models, and reads coverage and junctions from
#' the supplied BAM files over those windows. The result is one locus per
#' event, ready for [plot_sashimi_all()].
#'
#' @details
#' The two isoforms of an event are drawn as two model rows: the inclusion form
#' (role `"main"`) carries the alternative exon, the skipping form (role
#' `"alt"`) does not. That is the same reading as `rmats2sashimiplot`'s
#' two-row layout, but as data rather than as pixels, so the rows can be
#' recoloured or relabelled afterwards.
#'
#' @param events Path to an rMATS output file, or a data frame from
#'   [read_rmats()].
#' @param bam A manifest or vector of BAM paths, as for [sashimi_from_bam()].
#' @param event_type One of [event_types()]; inferred when `NULL`.
#' @param group_col,label_col Manifest columns for grouping and labels.
#' @param psi_from `"rmats"` takes the inclusion levels rMATS reported;
#'   `"junctions"` recomputes PSI from the junctions read out of the BAMs;
#'   `"none"` leaves the slot empty.
#' @param ... Passed to [sashimi_from_bam()] (`bin`, `strand`, `min_count`,
#'   `min_mapq`, and so on).
#' @inheritParams read_rmats
#'
#' @return A `sashimi_data` object with one locus per event.
#'
#' @examples
#' ev <- system.file("extdata", "SE.MATS.JC.txt", package = "omakase")
#' bams <- system.file("extdata", "samples.tsv", package = "omakase")
#' if (nzchar(ev) && nzchar(bams)) {
#'   sd <- sashimi_from_rmats(ev, bams, min_count = 5)
#'   plot_sashimi(sd)
#' }
#'
#' @export
sashimi_from_rmats <- function(events, bam, event_type = NULL, fdr = NULL,
                               min_dpsi = NULL, top = NULL, genes = NULL,
                               chrom_style = "keep", flank = 300,
                               group_col = NULL, label_col = NULL,
                               psi_from = c("rmats", "junctions", "none"),
                               ...) {
  psi_from <- match.arg(psi_from)
  ev <- if (is.data.frame(events)) {
    as_df(events)
  } else {
    read_rmats(events, event_type, fdr = fdr, min_dpsi = min_dpsi, top = top,
               genes = genes, chrom_style = chrom_style, flank = flank)
  }
  if (!nrow(ev)) om_abort("No events left after filtering.")

  regions <- data.frame(
    locus_id = ev$event_id,
    name = ifelse(is.na(ev$gene_name), ev$event_id, ev$gene_name),
    chrom = ev$chrom, start = ev$win_lo, end = ev$win_hi, strand = ev$strand,
    stringsAsFactors = FALSE
  )
  x <- sashimi_from_bam(bam, regions, group_col = group_col,
                        label_col = label_col, ...)

  # Replace whatever models the annotation supplied with the event's own
  # two-isoform picture, which is the thing the figure is about.
  x$models <- normalise_slot(rmats_models(ev), "models")
  x$meta$source <- "rmats"
  x$meta$event_type <- unique(ev$event_type)

  if (identical(psi_from, "rmats")) {
    groups <- order_levels(unique(x$tracks$group), x$meta$group_order)
    psi <- rmats_psi(ev, groups)
    if (!is.null(psi)) x$psi <- normalise_slot(psi, "psi")
  } else if (identical(psi_from, "junctions")) {
    x <- label_event_junctions(x, ev)
    x <- compute_psi(x, main = "main", alt = "alt")
  }
  attr(x, "events") <- ev
  x
}

# Two model rows per event: the inclusion isoform keeps the alternative exon,
# the skipping isoform drops it.
#' @noRd
rmats_models <- function(ev) {
  parts <- lapply(seq_len(nrow(ev)), function(i) {
    e <- ev[i, ]
    layout <- RMATS_LAYOUT[[e$event_type]]
    nms <- names(layout$exons)
    box <- function(nm) data.frame(
      start = e[[paste0(nm, "_start")]], end = e[[paste0(nm, "_end")]]
    )
    inc <- rbind_all(lapply(nms, box))
    # The skipping form is everything except the alternative exon; for MXE the
    # two forms differ in which of the two mutually exclusive exons they carry.
    skip_nms <- if (identical(e$event_type, "MXE")) {
      setdiff(nms, "alt")
    } else {
      setdiff(nms, layout$alt_exon)
    }
    skp <- rbind_all(lapply(skip_nms, box))

    rbind(
      data.frame(locus_id = e$event_id, tx_id = paste0(e$event_id, ".inclusion"),
                 role = "main", feature = "exon", inc, strand = e$strand,
                 stringsAsFactors = FALSE),
      data.frame(locus_id = e$event_id, tx_id = paste0(e$event_id, ".skipping"),
                 role = "alt", feature = "exon", skp, strand = e$strand,
                 stringsAsFactors = FALSE)
    )
  })
  rbind_all(parts)
}

# rMATS reports PSI for exactly two sample groups; map them onto the first two
# groups present in the manifest.
#' @noRd
rmats_psi <- function(ev, groups) {
  if (length(groups) < 2) return(NULL)
  rbind(
    data.frame(locus_id = ev$event_id, group = groups[1], psi = ev$psi_1,
               stringsAsFactors = FALSE),
    data.frame(locus_id = ev$event_id, group = groups[2], psi = ev$psi_2,
               stringsAsFactors = FALSE)
  )
}

# Mark junctions read from the BAMs as inclusion or skipping by matching their
# coordinates against the event's exon boundaries, so PSI can be recomputed
# from the alignments rather than taken on trust.
#' @noRd
label_event_junctions <- function(x, ev, tol = 5) {
  if (!nrow(x$junctions)) return(x)
  j <- x$junctions
  j$role <- NA_character_

  for (i in seq_len(nrow(ev))) {
    e <- ev[i, ]
    sel <- j$locus_id == e$event_id
    if (!any(sel)) next
    a_s <- e$alt_start; a_e <- e$alt_end
    up_e <- e$upstream_end %||% NA
    dn_s <- e$downstream_start %||% NA

    near <- function(u, v) !is.na(u) & !is.na(v) & abs(u - v) <= tol
    # Inclusion junctions touch the alternative exon; the skipping junction
    # jumps straight from the upstream exon to the downstream one.
    inc <- sel & (near(j$x1, a_s) | near(j$x0, a_e))
    skp <- sel & near(j$x0, up_e) & near(j$x1, dn_s)
    j$role[inc] <- "main"
    j$role[skp] <- "alt"
  }
  x$junctions <- j
  x
}
