# ---------------------------------------------------------------------------
# Open reading frames.
#
# To say what an alternative start site does to a transcript you need the
# protein each isoform encodes, which means calling an ORF on the spliced mRNA
# rather than on the genome. These are the pieces that do that: splice the
# exons, find the ORF, count upstream AUGs.
# ---------------------------------------------------------------------------

#' Splice exons into an mRNA sequence
#'
#' Concatenates exon sequences in transcription order and reverse-complements
#' them for a minus-strand transcript, giving the mature mRNA an ORF should be
#' called on.
#'
#' @param exons A data frame with `start` and `end`, one row per exon.
#' @param strand `"+"` or `"-"`.
#' @param genome A `BSgenome` object, a `DNAStringSet` of chromosomes, or a
#'   path to an indexed FASTA.
#' @param chrom The contig the exons lie on.
#'
#' @return A `DNAString`.
#'
#' @examples
#' if (requireNamespace("Biostrings", quietly = TRUE)) {
#'   gen <- Biostrings::DNAStringSet(c(chr1 = strrep("ATGC", 100)))
#'   splice_mrna(data.frame(start = c(1, 51), end = c(20, 70)), "+", gen, "chr1")
#' }
#'
#' @export
splice_mrna <- function(exons, strand, genome, chrom) {
  need_pkg("Biostrings", "calling open reading frames")
  e <- as_df(exons)
  require_cols(e, c("start", "end"), "exon table")
  e <- e[order(e$start), , drop = FALSE]

  seqs <- lapply(seq_len(nrow(e)), function(i) {
    get_seq(genome, chrom, e$start[i], e$end[i])
  })
  mrna <- Reduce(function(a, b) Biostrings::xscat(a, b), seqs)
  mrna <- Biostrings::DNAString(as.character(mrna))
  if (identical(strand, "-")) mrna <- Biostrings::reverseComplement(mrna)
  mrna
}

#' @noRd
get_seq <- function(genome, chrom, start, end) {
  if (is.character(genome) && length(genome) == 1L) {
    need_pkg("Rsamtools", "reading a FASTA genome")
    gr <- region_to_granges(new_region(chrom, start, end, "*"))
    return(Rsamtools::scanFa(genome, gr)[[1]])
  }
  if (inherits(genome, "BSgenome")) {
    need_pkg("BSgenome", "reading a BSgenome")
    return(BSgenome::getSeq(genome, chrom, start = start, end = end))
  }
  # DNAStringSet or list of DNAString
  if (is.null(genome[[chrom]])) {
    om_abort("Genome has no sequence named {.val {chrom}}.")
  }
  Biostrings::subseq(genome[[chrom]], start = start, end = end)
}

#' Call the open reading frame of an mRNA
#'
#' Finds the ORF and returns the peptide it encodes together with the length of
#' the 5' untranslated region and the number of upstream AUGs.
#'
#' @details
#' Two conventions are available. `"first"` takes the first AUG in the
#' sequence that opens a reading frame terminated by a stop codon, which is the
#' scanning model and the one used to classify start-site consequences.
#' `"longest"` takes the AUG that yields the longest peptide, which is more
#' robust when the 5' end of the transcript is not trustworthy.
#'
#' @param mrna A `DNAString`, `DNAStringSet` of length 1, or a character
#'   string.
#' @param rule `"first"` or `"longest"`.
#' @param min_aa Ignore reading frames shorter than this many amino acids.
#' @param require_stop Require an in-frame stop codon. When `FALSE`, a frame
#'   running off the 3' end still counts, which matters for transcripts whose
#'   3' end was not sequenced.
#'
#' @return A list with `protein` (character, `""` when no ORF was found),
#'   `start`, `end`, `utr5_len`, `aa_len` and `n_uATG`.
#'
#' @examples
#' find_orf("AAAATGGGGCCCTAA")
#'
#' @export
find_orf <- function(mrna, rule = c("first", "longest"), min_aa = 10,
                     require_stop = TRUE) {
  need_pkg("Biostrings", "calling open reading frames")
  rule <- match.arg(rule)
  s <- if (is.character(mrna)) mrna else as.character(mrna)
  s <- toupper(gsub("[^ACGTN]", "N", s))
  n <- nchar(s)

  empty <- list(protein = "", start = NA_integer_, end = NA_integer_,
                utr5_len = NA_integer_, aa_len = 0L, n_uATG = 0L)
  if (n < 6) return(empty)

  starts <- gregexpr("ATG", s, fixed = TRUE)[[1]]
  starts <- starts[starts > 0]
  if (!length(starts)) return(empty)

  best <- NULL
  for (st in starts) {
    # Trim to a whole number of codons so translate() has nothing to complain
    # about.
    avail <- n - st + 1
    ncod <- avail %/% 3
    if (ncod < min_aa + 1) next
    sub <- substr(s, st, st + ncod * 3 - 1)
    aa <- as.character(Biostrings::translate(
      Biostrings::DNAString(sub), if.fuzzy.codon = "solve"
    ))
    stop_at <- regexpr("*", aa, fixed = TRUE)
    has_stop <- stop_at > 0
    if (require_stop && !has_stop) next
    pep <- if (has_stop) substr(aa, 1, stop_at - 1) else aa
    if (nchar(pep) < min_aa) next

    cand <- list(
      protein = pep, start = st,
      end = st + (nchar(pep) + as.integer(has_stop)) * 3 - 1,
      utr5_len = st - 1L, aa_len = nchar(pep)
    )
    if (identical(rule, "first")) {
      best <- cand
      break
    }
    if (is.null(best) || cand$aa_len > best$aa_len) best <- cand
  }
  if (is.null(best)) return(empty)

  best$n_uATG <- sum(starts < best$start)
  best
}

#' Count upstream AUGs in a 5' UTR
#'
#' The number of AUG triplets before the main start codon. An AUG gained in a
#' lengthened 5' UTR can create an upstream open reading frame and suppress
#' translation of the main one, so this is the quantity that makes a "longer
#' 5' UTR" consequence mean something.
#'
#' @param mrna A `DNAString` or character string.
#' @param orf_start 1-based position of the main start codon.
#' @return An integer count.
#' @examples
#' count_uatg("ATGCCCAAAATGGGG", orf_start = 10)
#' @export
count_uatg <- function(mrna, orf_start) {
  s <- if (is.character(mrna)) mrna else as.character(mrna)
  if (is.na(orf_start) || orf_start <= 3) return(0L)
  utr <- substr(toupper(s), 1, orf_start - 1)
  m <- gregexpr("ATG", utr, fixed = TRUE)[[1]]
  sum(m > 0)
}

#' Call ORFs for a set of transcript models
#'
#' Convenience wrapper that splices each transcript and calls its ORF, giving
#' the per-isoform table [classify_consequence()] consumes.
#'
#' @param models A data frame with `tx_id`, `chrom`, `strand`, `start`, `end`,
#'   one row per exon.
#' @param genome A `BSgenome`, `DNAStringSet`, or path to an indexed FASTA.
#' @param rule,min_aa,require_stop Passed to [find_orf()].
#'
#' @return A data frame with one row per transcript: `tx_id`, `chrom`,
#'   `strand`, `tx_start`, `tx_end`, `first_exon_start`, `first_exon_end`,
#'   `tss_pos`, `protein`, `aa_len`, `utr5_len`, `n_uATG`.
#'
#' @examples
#' if (requireNamespace("Biostrings", quietly = TRUE)) {
#'   # A toy genome: a start codon, a short coding stretch, then a stop.
#'   genome <- Biostrings::DNAStringSet(c(
#'     chr1 = paste0(strrep("T", 60), "ATG", strrep("GCA", 30), "TAA",
#'                   strrep("T", 60))
#'   ))
#'   models <- data.frame(
#'     tx_id = "tx1", chrom = "chr1", strand = "+",
#'     start = c(1, 100), end = c(80, 214)
#'   )
#'   orf_table(models, genome)
#' }
#'
#' @export
orf_table <- function(models, genome, rule = "first", min_aa = 10,
                      require_stop = TRUE) {
  m <- as_df(models)
  require_cols(m, c("tx_id", "chrom", "strand", "start", "end"), "model table")

  parts <- lapply(split(m, m$tx_id), function(d) {
    strand <- d$strand[1]
    chrom <- d$chrom[1]
    d <- d[order(d$start), , drop = FALSE]
    # The first exon in transcription order, which is the last one by
    # coordinate on the minus strand.
    fe <- if (identical(strand, "-")) d[nrow(d), ] else d[1, ]
    tss <- if (identical(strand, "-")) fe$end else fe$start

    orf <- tryCatch({
      mrna <- splice_mrna(d, strand, genome, chrom)
      find_orf(mrna, rule = rule, min_aa = min_aa, require_stop = require_stop)
    }, error = function(e) {
      om_warn("ORF calling failed for {.val {d$tx_id[1]}}: {conditionMessage(e)}")
      list(protein = "", start = NA, end = NA, utr5_len = NA, aa_len = 0L,
           n_uATG = 0L)
    })

    data.frame(
      tx_id = d$tx_id[1], chrom = chrom, strand = strand,
      tx_start = min(d$start), tx_end = max(d$end),
      first_exon_start = fe$start, first_exon_end = fe$end,
      tss_pos = tss,
      protein = orf$protein, aa_len = orf$aa_len,
      utr5_len = orf$utr5_len, n_uATG = orf$n_uATG,
      stringsAsFactors = FALSE
    )
  })
  out <- rbind_all(parts)
  rownames(out) <- NULL
  out
}
