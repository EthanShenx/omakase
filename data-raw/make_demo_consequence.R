#!/usr/bin/env Rscript
# Build inst/extdata/demo_consequence.tsv.
#
# The switches here are synthetic. Each row is a made-up pair of isoforms -
# a peptide, a 5' UTR length, a first exon and a transcript span - which is
# then put through classify_consequence() exactly as real data would be. So
# the shipped table is genuinely the classifier's own output, and the demo
# exercises the real decision tree, while no unpublished result leaves the
# machine it was computed on.

devtools::load_all(".", quiet = TRUE)
set.seed(20260812)

n <- 180
aa <- c("A", "R", "N", "D", "C", "Q", "E", "G", "H", "I", "L", "K", "M",
        "F", "P", "S", "T", "W", "Y", "V")
pep <- function(k) paste0("M", paste(sample(aa, k - 1, replace = TRUE), collapse = ""))

# Intended mix of outcomes, chosen to resemble a real start-site switch set:
# most switches touch the N terminus, a minority are purely regulatory.
kind <- sample(
  c("utr_longer", "utr_shorter", "utr_equal", "promoter_swap",
    "alt_nterm", "extension", "truncation", "orf_loss", "distal"),
  n, replace = TRUE,
  prob = c(0.11, 0.06, 0.03, 0.12, 0.34, 0.18, 0.11, 0.02, 0.03)
)

rows <- lapply(seq_len(n), function(i) {
  k <- kind[i]
  strand <- sample(c("+", "-"), 1)
  base <- pep(sample(90:420, 1))
  m_utr <- sample(40:400, 1)

  # Defaults describe a switch that changes nothing; each branch perturbs
  # only what its consequence requires.
  a_pep <- base
  a_utr <- m_utr
  m_fe <- c(1000, 1000 + sample(120:400, 1))
  a_fe <- m_fe + sample(20:90, 1)
  m_tx <- c(1000, 1000 + sample(8000:40000, 1))
  a_tx <- c(a_fe[1], m_tx[2])
  m_tss <- if (strand == "-") m_fe[2] else m_fe[1]
  a_tss <- if (strand == "-") a_fe[2] else a_fe[1]

  if (k == "utr_longer") {
    a_utr <- m_utr + sample(30:600, 1)
  } else if (k == "utr_shorter") {
    a_utr <- max(5, m_utr - sample(20:min(300, m_utr - 5), 1))
  } else if (k == "utr_equal") {
    a_utr <- m_utr
  } else if (k == "promoter_swap") {
    # A first exon that does not touch the reference one, and lies downstream.
    a_fe <- c(m_fe[2] + sample(2000:20000, 1), 0)
    a_fe[2] <- a_fe[1] + sample(100:350, 1)
    a_tx <- c(a_fe[1], m_tx[2])
    a_tss <- if (strand == "-") a_fe[2] else a_fe[1]
    if (strand == "-") {
      # Downstream on the minus strand means a lower coordinate.
      a_fe <- c(m_fe[1] - sample(2000:20000, 1), 0)
      a_fe[2] <- a_fe[1] + sample(100:350, 1)
      a_tss <- a_fe[2]
      a_tx <- c(a_fe[1], m_tx[2])
    }
  } else if (k == "alt_nterm") {
    a_pep <- paste0(pep(sample(8:60, 1)), substring(base, sample(20:60, 1)))
  } else if (k == "extension") {
    a_pep <- paste0(pep(sample(5:45, 1)), base)
  } else if (k == "truncation") {
    a_pep <- substring(base, sample(10:70, 1))
  } else if (k == "orf_loss") {
    a_pep <- ""
  } else if (k == "distal") {
    # Two transcripts sharing a gene label but not a locus.
    a_tx <- c(m_tx[2] + sample(5000:50000, 1), 0)
    a_tx[2] <- a_tx[1] + sample(5000:30000, 1)
    a_fe <- c(a_tx[1], a_tx[1] + 200)
    a_tss <- a_fe[1]
  }

  data.frame(
    gene_id = sprintf("ENSDEMOG%08d", i),
    gene_name = sprintf("Demo%03d", i),
    strand = strand,
    main_tx = sprintf("DEMO.%d.main", i),
    alt_tx = sprintf("DEMO.%d.alt", i),
    main_protein = base, alt_protein = a_pep,
    main_utr5_len = m_utr, alt_utr5_len = a_utr,
    main_aa_len = nchar(base), alt_aa_len = nchar(a_pep),
    main_n_uATG = sample(0:2, 1),
    alt_n_uATG = sample(0:4, 1),
    main_first_exon_start = m_fe[1], main_first_exon_end = m_fe[2],
    alt_first_exon_start = a_fe[1], alt_first_exon_end = a_fe[2],
    main_tss = m_tss, alt_tss = a_tss,
    main_tx_start = m_tx[1], main_tx_end = m_tx[2],
    alt_tx_start = a_tx[1], alt_tx_end = a_tx[2],
    both_full_length = as.integer(stats::runif(1) < 0.75),
    stringsAsFactors = FALSE
  )
})

sw <- do.call(rbind, rows)
out <- classify_consequence(sw)

# Ship the classified table without the peptide columns, which are long and
# add nothing to a plotting demo.
keep <- c("gene_id", "gene_name", "strand", "category", "subtype",
          "d_utr5_bp", "n_uATG_gained", "uorf_gained", "body_overlap",
          "main_aa_len", "alt_aa_len", "aa_delta",
          "main_tx", "alt_tx", "both_full_length")
out <- out[order(out$category, out$subtype, out$gene_name), keep]

path <- "inst/extdata/demo_consequence.tsv"
utils::write.table(out, path, sep = "\t", quote = FALSE, row.names = FALSE)

cat("wrote", path, "with", nrow(out), "rows\n\n")
print(table(out$category, out$subtype))
