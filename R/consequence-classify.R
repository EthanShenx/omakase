# ---------------------------------------------------------------------------
# What an alternative start site does to the transcript.
#
# Given the reference isoform and the alternative one, compare the proteins
# they encode and the 5' UTRs that precede them, and put the switch in one of
# three classes. The rule that "a 5' UTR change wins" matters: if the protein
# is untouched, the consequence is regulatory, however different the transcript
# looks in the genome browser.
# ---------------------------------------------------------------------------

#' Consequence categories and subtypes
#'
#' @return A named list with `categories` and `subtypes`, in the order they are
#'   drawn.
#' @examples
#' consequence_levels()
#' @export
consequence_levels <- function() {
  list(
    categories = c("5'UTR change", "Promoter swap", "N-terminal/CDS",
                   "Unclassified", "Distal (excluded)"),
    subtypes = c("longer", "shorter", "equal", "alt first exon",
                 "alt N-terminus", "N-term extension", "N-term truncation",
                 "ORF loss", "ORF gain", "no ORF in either",
                 "non-overlapping isoforms")
  )
}

#' Human-readable labels for consequence subtypes
#'
#' @return A named character vector mapping subtype codes to display labels.
#' @examples
#' consequence_labels()[["longer"]]
#' @export
consequence_labels <- function() {
  c(
    "longer" = "Longer 5'UTR",
    "shorter" = "Shorter 5'UTR",
    "equal" = "Equal-length 5'UTR",
    "alt first exon" = "Alternative first exon",
    "alt N-terminus" = "Alternative N-terminus",
    "N-term extension" = "N-terminal extension",
    "N-term truncation" = "N-terminal truncation",
    "ORF loss" = "ORF loss",
    "ORF gain" = "ORF gain",
    "no ORF in either" = "No ORF in either",
    "non-overlapping isoforms" = "Non-overlapping isoforms"
  )
}

#' Classify the consequence of an alternative start site
#'
#' Compares a reference ("main") isoform with an alternative one and reports
#' what the switch between them does to the transcript.
#'
#' @details
#' Writing \eqn{P_m} and \eqn{P_a} for the peptides the two isoforms encode,
#' and \eqn{\ell_m}, \eqn{\ell_a} for their 5' UTR lengths, the decision tree
#' is:
#'
#' \enumerate{
#'   \item \eqn{P_m \ne \emptyset,\ P_a = \emptyset} gives **ORF loss**;
#'     \eqn{P_m = \emptyset,\ P_a \ne \emptyset} gives **ORF gain**; both empty
#'     is unclassified. Calling a pair with no ORF in either isoform an "ORF
#'     loss" would blame the start site for a failure that was already there.
#'   \item \eqn{P_a = P_m}, so the protein is untouched and the change is
#'     purely 5': if the two first exons are disjoint **and** the alternative
#'     site lies downstream, it is a **promoter swap**; otherwise a **5' UTR
#'     change**, `longer`, `shorter` or `equal` by the sign of
#'     \eqn{\ell_a - \ell_m}.
#'   \item \eqn{P_a \ne P_m}: if \eqn{P_m} ends with \eqn{P_a} the alternative
#'     protein is a truncation; if \eqn{P_a} ends with \eqn{P_m} it is an
#'     extension; otherwise the N-terminus is simply different.
#' }
#'
#' A guard runs first. The reciprocal overlap of the two transcript spans,
#' \deqn{o = \frac{\max(0, \min(e_1, e_2) - \max(s_1, s_2))}{\min(e_1 - s_1,\ e_2 - s_2)},}
#' is zero when the two "isoforms" do not overlap at all. That is not an
#' alternative start site of one transcription unit, it is two separate
#' transcripts sharing a gene label, and such pairs are marked
#' `"Distal (excluded)"` rather than classified.
#'
#' @param x A data frame with one row per switch. Required columns:
#'   `main_protein`, `alt_protein`, `main_utr5_len`, `alt_utr5_len`,
#'   `main_first_exon_start`, `main_first_exon_end`, `alt_first_exon_start`,
#'   `alt_first_exon_end`, `main_tss`, `alt_tss`, `strand`. Optional:
#'   `gene_id`, `gene_name`, `main_tx_start`, `main_tx_end`, `alt_tx_start`,
#'   `alt_tx_end` (needed for the overlap guard), `main_n_uATG`,
#'   `alt_n_uATG`, `main_aa_len`, `alt_aa_len`.
#' @param min_overlap Pairs whose reciprocal transcript-span overlap is at or
#'   below this are excluded as distal mispairings.
#'
#' @return The input with added columns: `category`, `subtype`, `body_overlap`,
#'   `d_utr5_bp`, `aa_delta`, `n_uATG_gained` and `uorf_gained`.
#'
#' @examples
#' sw <- data.frame(
#'   gene_name = c("A", "B", "C"),
#'   strand = "+",
#'   main_protein = c("MKVLA", "MKVLA", "MKVLA"),
#'   alt_protein  = c("MKVLA", "VLA",   ""),
#'   main_utr5_len = c(100, 100, 100), alt_utr5_len = c(200, 100, 100),
#'   main_first_exon_start = 1000, main_first_exon_end = 1200,
#'   alt_first_exon_start = c(1050, 1050, 1050),
#'   alt_first_exon_end = c(1250, 1250, 1250),
#'   main_tss = 1000, alt_tss = 1050,
#'   main_tx_start = 1000, main_tx_end = 5000,
#'   alt_tx_start = 1050, alt_tx_end = 5000
#' )
#' classify_consequence(sw)[, c("gene_name", "category", "subtype")]
#'
#' @export
classify_consequence <- function(x, min_overlap = 0) {
  d <- as_df(x)
  need <- c("main_protein", "alt_protein", "main_utr5_len", "alt_utr5_len",
            "main_first_exon_start", "main_first_exon_end",
            "alt_first_exon_start", "alt_first_exon_end",
            "main_tss", "alt_tss", "strand")
  require_cols(d, need, "switch table")

  na_to <- function(v, fill) ifelse(is.na(v), fill, v)
  Pm <- na_to(as.character(d$main_protein), "")
  Pa <- na_to(as.character(d$alt_protein), "")
  lm_ <- as.numeric(d$main_utr5_len)
  la <- as.numeric(d$alt_utr5_len)

  downstream <- ifelse(d$strand == "-", d$alt_tss < d$main_tss,
                       d$alt_tss > d$main_tss)
  disjoint <- !intervals_overlap(d$main_first_exon_start, d$main_first_exon_end,
                                 d$alt_first_exon_start, d$alt_first_exon_end)

  category <- character(nrow(d))
  subtype <- character(nrow(d))

  for (i in seq_len(nrow(d))) {
    r <- classify_one(Pm[i], Pa[i], lm_[i], la[i], disjoint[i], downstream[i])
    category[i] <- r[1]
    subtype[i] <- r[2]
  }

  # Distal guard, applied after classification so the reason for exclusion is
  # explicit rather than hidden in a filter upstream.
  ov <- rep(NA_real_, nrow(d))
  if (all(c("main_tx_start", "main_tx_end", "alt_tx_start", "alt_tx_end")
          %in% names(d))) {
    ov <- reciprocal_overlap(d$main_tx_start, d$main_tx_end,
                             d$alt_tx_start, d$alt_tx_end)
    distal <- !is.na(ov) & ov <= min_overlap
    category[distal] <- "Distal (excluded)"
    subtype[distal] <- "non-overlapping isoforms"
  }

  d$category <- category
  d$subtype <- subtype
  d$body_overlap <- ov
  d$d_utr5_bp <- la - lm_
  if (all(c("main_aa_len", "alt_aa_len") %in% names(d))) {
    d$aa_delta <- as.numeric(d$alt_aa_len) - as.numeric(d$main_aa_len)
  }
  if (all(c("main_n_uATG", "alt_n_uATG") %in% names(d))) {
    d$n_uATG_gained <- as.numeric(d$alt_n_uATG) - as.numeric(d$main_n_uATG)
    d$uorf_gained <- as.integer(d$n_uATG_gained >= 1)
  }
  d
}

# The decision tree for a single switch. Kept as a scalar function because the
# branches are easier to check against the specification this way than as a
# stack of vectorised ifelse() calls.
#' @noRd
classify_one <- function(Pm, Pa, lm_, la, disjoint, downstream) {
  if (nzchar(Pm) && !nzchar(Pa)) return(c("N-terminal/CDS", "ORF loss"))
  if (!nzchar(Pm) && nzchar(Pa)) return(c("N-terminal/CDS", "ORF gain"))
  if (!nzchar(Pm) && !nzchar(Pa)) return(c("Unclassified", "no ORF in either"))

  if (identical(Pa, Pm)) {
    if (isTRUE(disjoint) && isTRUE(downstream)) {
      return(c("Promoter swap", "alt first exon"))
    }
    if (!is.na(la) && !is.na(lm_)) {
      if (la > lm_) return(c("5'UTR change", "longer"))
      if (la < lm_) return(c("5'UTR change", "shorter"))
    }
    return(c("5'UTR change", "equal"))
  }
  if (endsWith(Pm, Pa)) return(c("N-terminal/CDS", "N-term truncation"))
  if (endsWith(Pa, Pm)) return(c("N-terminal/CDS", "N-term extension"))
  c("N-terminal/CDS", "alt N-terminus")
}

#' @noRd
intervals_overlap <- function(s1, e1, s2, e2) {
  !(e1 < s2 | e2 < s1)
}

#' @noRd
reciprocal_overlap <- function(s1, e1, s2, e2) {
  ov <- pmax(0, pmin(e1, e2) - pmax(s1, s2))
  ov / pmax(1, pmin(e1 - s1, e2 - s2))
}

#' Summarise a classified switch table
#'
#' Counts and proportions per category and subtype, in the order the donut
#' draws them.
#'
#' @param x A data frame from [classify_consequence()], or any table with
#'   `category` and `subtype`.
#' @param by `"subtype"` (the default) or `"category"`.
#' @param filter An optional expression evaluated in `x` to subset rows before
#'   counting - for example `both_full_length == 1`.
#' @param drop Categories to exclude. Defaults to the two non-biological ones.
#'
#' @return A data frame with `category`, `subtype` (when `by = "subtype"`),
#'   `n`, `proportion` and `label`.
#'
#' @examples
#' d <- data.frame(category = c("5'UTR change", "5'UTR change", "Promoter swap"),
#'                 subtype = c("longer", "shorter", "alt first exon"))
#' consequence_summary(d)
#'
#' @export
consequence_summary <- function(x, by = c("subtype", "category"),
                                filter = NULL,
                                drop = c("Unclassified", "Distal (excluded)")) {
  by <- match.arg(by)
  d <- as_df(x)
  require_cols(d, c("category", by), "consequence table")

  f <- rlang::enquo(filter)
  if (!rlang::quo_is_null(f)) {
    keep <- rlang::eval_tidy(f, data = d)
    d <- d[!is.na(keep) & keep, , drop = FALSE]
  }
  if (length(drop)) d <- d[!d$category %in% drop, , drop = FALSE]
  if (!nrow(d)) om_abort("No rows left to summarise.")

  lv <- consequence_levels()
  if (identical(by, "category")) {
    tab <- as.data.frame(table(category = d$category), stringsAsFactors = FALSE)
    names(tab)[2] <- "n"
    tab$category <- factor(tab$category,
                           levels = intersect(lv$categories, tab$category))
    tab <- tab[order(tab$category), , drop = FALSE]
    tab$label <- as.character(tab$category)
  } else {
    tab <- as.data.frame(table(category = d$category, subtype = d$subtype),
                         stringsAsFactors = FALSE)
    names(tab)[3] <- "n"
    tab <- tab[tab$n > 0, , drop = FALSE]
    # Subtype order, not count order: the donut should read the same way from
    # one dataset to the next.
    tab$subtype <- factor(tab$subtype,
                          levels = intersect(lv$subtypes, tab$subtype))
    tab$category <- factor(tab$category,
                           levels = intersect(lv$categories, tab$category))
    tab <- tab[order(tab$subtype), , drop = FALSE]
    tab$label <- unname(consequence_labels()[as.character(tab$subtype)])
    tab$label[is.na(tab$label)] <- as.character(tab$subtype)[is.na(tab$label)]
  }
  tab$proportion <- 100 * tab$n / sum(tab$n)
  rownames(tab) <- NULL
  tab
}
