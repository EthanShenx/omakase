# ---------------------------------------------------------------------------
# Transcript annotation.
#
# GTF and GFF3 disagree about almost everything except that exons have parents,
# so both are read through rtracklayer and then reduced to the handful of
# columns the model panel needs.
# ---------------------------------------------------------------------------

#' Read a GTF or GFF3 annotation
#'
#' Imports exon and CDS records and returns the tidy transcript models the
#' annotation panel draws. Format is detected from the file extension, and
#' gzip-compressed files are handled transparently.
#'
#' @param path Path to a `.gtf`, `.gff`, `.gff3` file, optionally gzipped.
#' @param region Optional region to restrict the import to. Supplying one is
#'   much faster than reading a whole-genome annotation, and requires a
#'   tabix-indexed file for the fastest path.
#' @param feature Which record types to keep. `"exon"` is enough to draw a
#'   transcript; adding `"CDS"` lets the panel draw coding regions taller than
#'   untranslated ones.
#' @param gene Restrict to one or more gene names.
#' @param transcript Restrict to one or more transcript identifiers.
#'
#' @return A data frame with columns `tx_id`, `gene_id`, `gene_name`, `chrom`,
#'   `start`, `end`, `strand` and `feature`.
#'
#' @examples
#' gtf <- system.file("extdata", "annotation.gtf", package = "omakase")
#' if (nzchar(gtf)) head(read_annotation(gtf))
#'
#' @export
read_annotation <- function(path, region = NULL, feature = c("exon", "CDS"),
                            gene = NULL, transcript = NULL) {
  need_pkg("rtracklayer", "reading GTF/GFF annotation")
  if (inherits(path, "data.frame")) return(path)
  if (!file.exists(path)) om_abort("Annotation not found: {.path {path}}.")

  args <- list(con = path)
  if (!is.null(region)) {
    r <- parse_region(region)
    if (!is.na(r$start)) args$which <- region_to_granges(new_region(
      r$chrom, r$start, r$end, "*"))
  }
  gr <- tryCatch(do.call(rtracklayer::import, args), error = function(e) {
    # `which` needs an index; without one, read the file and filter afterwards.
    if (!is.null(args$which)) {
      g <- rtracklayer::import(path)
      return(IRanges::subsetByOverlaps(g, args$which))
    }
    om_abort(c("Could not read {.path {path}}.", "x" = conditionMessage(e)))
  })
  if (!length(gr)) return(empty_models())

  md <- S4Vectors::mcols(gr)
  type <- as.character(col(md, "type") %||% NA)
  gr <- gr[type %in% feature]
  if (!length(gr)) return(empty_models())
  md <- S4Vectors::mcols(gr)

  # GTF says transcript_id, GFF3 says Parent or ID; take whichever is present.
  tx <- first_present(md, c("transcript_id", "Parent", "ID", "transcript_name"))
  gid <- first_present(md, c("gene_id", "geneID", "gene"))
  gnm <- first_present(md, c("gene_name", "Name", "gene"))

  out <- data.frame(
    tx_id = as.character(tx),
    gene_id = as.character(gid),
    gene_name = as.character(gnm),
    chrom = as.character(GenomicRanges::seqnames(gr)),
    start = GenomicRanges::start(gr),
    end = GenomicRanges::end(gr),
    strand = as.character(GenomicRanges::strand(gr)),
    feature = as.character(S4Vectors::mcols(gr)$type),
    stringsAsFactors = FALSE
  )
  if (!is.null(gene)) {
    out <- out[out$gene_name %in% gene | out$gene_id %in% gene, , drop = FALSE]
  }
  if (!is.null(transcript)) {
    out <- out[out$tx_id %in% transcript, , drop = FALSE]
  }
  out <- out[!is.na(out$tx_id), , drop = FALSE]
  rownames(out) <- NULL
  out
}

#' @noRd
empty_models <- function() {
  data.frame(tx_id = character(0), gene_id = character(0),
             gene_name = character(0), chrom = character(0),
             start = numeric(0), end = numeric(0), strand = character(0),
             feature = character(0), stringsAsFactors = FALSE)
}

# Pull the first metadata column that exists, unlisting the CharacterList that
# GFF3 Parent attributes come back as.
#' @noRd
first_present <- function(md, names) {
  for (nm in names) {
    if (nm %in% names(md)) {
      v <- md[[nm]]
      if (inherits(v, "List") || is.list(v)) {
        v <- vapply(v, function(z) if (length(z)) as.character(z)[1] else NA_character_,
                    character(1))
      }
      return(v)
    }
  }
  rep(NA_character_, nrow(md))
}

#' @noRd
resolve_annotation <- function(annotation, gene = NULL) {
  if (is.null(annotation)) return(NULL)
  if (is.data.frame(annotation)) {
    a <- as_df(annotation)
    if (!is.null(gene)) {
      a <- a[a$gene_name %in% gene | (col(a, "gene_id") %||% "") %in% gene, , drop = FALSE]
    }
    return(a)
  }
  read_annotation(annotation, gene = gene)
}

# Slice an annotation to a window and cast it into the `models` slot shape.
# Transcripts are kept whole if any exon overlaps, so a model is never drawn
# with its ends silently missing.
#' @noRd
models_in_region <- function(ann, r) {
  if (is.null(ann) || !nrow(ann)) return(empty_df(SLOT_SPEC$models))
  hit <- ann$chrom == r$chrom & ann$start <= r$end & ann$end >= r$start
  keep_tx <- unique(ann$tx_id[hit])
  a <- ann[ann$tx_id %in% keep_tx & ann$chrom == r$chrom, , drop = FALSE]
  if (!nrow(a)) return(empty_df(SLOT_SPEC$models))

  data.frame(
    tx_id = a$tx_id,
    role = NA_character_,
    feature = a$feature,
    start = a$start,
    end = a$end,
    strand = a$strand,
    stringsAsFactors = FALSE
  )
}

#' Add transcript models to a sashimi data object
#'
#' Attaches models from an annotation to loci that already carry coverage,
#' which is useful when tracks came from a source that knows nothing about
#' transcripts - a tag BED, or a table of precomputed counts.
#'
#' @param x A `sashimi_data` object.
#' @param annotation A GTF/GFF3 path or a data frame from [read_annotation()].
#' @param gene Restrict to these gene names.
#' @param max_tx Keep at most this many transcripts per locus, longest first.
#'   Annotations with thirty isoforms make an unreadable panel.
#' @return The object with its `models` slot populated.
#' @examples
#' sd <- sashimi_data(loci = data.frame(locus_id = "a", gene_name = "A",
#'   chrom = "chr1", strand = "+", win_lo = 1, win_hi = 100))
#' ann <- data.frame(tx_id = "t1", gene_id = "g1", gene_name = "A",
#'   chrom = "chr1", start = 10, end = 50, strand = "+", feature = "exon")
#' add_models(sd, ann)$models
#' @export
add_models <- function(x, annotation, gene = NULL, max_tx = NULL) {
  x <- as_sashimi_data(x)
  ann <- resolve_annotation(annotation, gene)
  if (is.null(ann) || !nrow(ann)) return(x)

  parts <- lapply(seq_len(nrow(x$loci)), function(i) {
    gi <- x$loci[i, ]
    r <- new_region(gi$chrom, gi$win_lo, gi$win_hi, gi$strand)
    m <- models_in_region(ann, r)
    if (!nrow(m)) return(NULL)
    if (!is.null(max_tx)) m <- keep_longest_tx(m, max_tx)
    m$locus_id <- gi$locus_id
    m
  })
  x$models <- normalise_slot(rbind_all(parts), "models")
  x
}

# The gene name and strand implied by the annotation records overlapping a
# window. Where several genes overlap, the one covering the most of the window
# wins - that is the gene the figure is about.
#' @noRd
annotation_hint <- function(ann, r) {
  if (is.null(ann) || !nrow(ann)) return(list(gene = NULL, strand = NULL))
  hit <- ann$chrom == r$chrom & ann$start <= r$end & ann$end >= r$start
  a <- ann[hit, , drop = FALSE]
  if (!nrow(a)) return(list(gene = NULL, strand = NULL))

  cov <- pmin(a$end, r$end) - pmax(a$start, r$start) + 1
  nm <- a$gene_name
  nm[is.na(nm)] <- a$gene_id[is.na(nm)]
  if (all(is.na(nm))) return(list(gene = NULL, strand = NULL))

  by_gene <- tapply(cov, nm, sum)
  best <- names(by_gene)[which.max(by_gene)]
  st <- unique(a$strand[nm == best])
  list(gene = best, strand = if (length(st) == 1L) st else NULL)
}

#' @noRd
keep_longest_tx <- function(m, max_tx) {
  spans <- tapply(m$end, m$tx_id, max) - tapply(m$start, m$tx_id, min)
  keep <- names(sort(spans, decreasing = TRUE))[seq_len(min(max_tx, length(spans)))]
  m[m$tx_id %in% keep, , drop = FALSE]
}
