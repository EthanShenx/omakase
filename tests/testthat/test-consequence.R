base_pep <- "MKVLAAGIVGTRSTQ"

one <- function(main_protein, alt_protein, main_utr5_len = 100,
                alt_utr5_len = 100, strand = "+",
                main_first_exon_start = 1000, main_first_exon_end = 1200,
                alt_first_exon_start = 1050, alt_first_exon_end = 1250,
                main_tss = 1000, alt_tss = 1050, ...) {
  d <- data.frame(
    main_protein = main_protein, alt_protein = alt_protein,
    main_utr5_len = main_utr5_len, alt_utr5_len = alt_utr5_len,
    strand = strand,
    main_first_exon_start = main_first_exon_start,
    main_first_exon_end = main_first_exon_end,
    alt_first_exon_start = alt_first_exon_start,
    alt_first_exon_end = alt_first_exon_end,
    main_tss = main_tss, alt_tss = alt_tss,
    stringsAsFactors = FALSE, ...
  )
  classify_consequence(d)
}

test_that("an unchanged protein with a longer 5' UTR is a 5'UTR change", {
  r <- one(base_pep, base_pep, main_utr5_len = 100, alt_utr5_len = 300)
  expect_equal(r$category, "5'UTR change")
  expect_equal(r$subtype, "longer")
  expect_equal(r$d_utr5_bp, 200)
})

test_that("a shorter and an equal 5' UTR are distinguished", {
  expect_equal(one(base_pep, base_pep, 300, 100)$subtype, "shorter")
  expect_equal(one(base_pep, base_pep, 100, 100)$subtype, "equal")
})

test_that("a disjoint downstream first exon is a promoter swap", {
  r <- one(base_pep, base_pep,
           main_first_exon_start = 1000, main_first_exon_end = 1200,
           alt_first_exon_start = 5000, alt_first_exon_end = 5200,
           main_tss = 1000, alt_tss = 5000)
  expect_equal(r$category, "Promoter swap")
  expect_equal(r$subtype, "alt first exon")
})

test_that("a disjoint but upstream first exon is not a promoter swap", {
  # The 5'-UTR rule wins: the protein is untouched, so the change is regulatory.
  r <- one(base_pep, base_pep, main_utr5_len = 100, alt_utr5_len = 400,
           main_first_exon_start = 5000, main_first_exon_end = 5200,
           alt_first_exon_start = 1000, alt_first_exon_end = 1200,
           main_tss = 5000, alt_tss = 1000)
  expect_equal(r$category, "5'UTR change")
})

test_that("downstream is read against the strand", {
  # On the minus strand a lower coordinate is downstream.
  r <- one(base_pep, base_pep, strand = "-",
           main_first_exon_start = 5000, main_first_exon_end = 5200,
           alt_first_exon_start = 1000, alt_first_exon_end = 1200,
           main_tss = 5200, alt_tss = 1200)
  expect_equal(r$category, "Promoter swap")
})

test_that("suffix relationships give truncation and extension", {
  expect_equal(one(base_pep, substring(base_pep, 5))$subtype,
               "N-term truncation")
  expect_equal(one(base_pep, paste0("MQQP", base_pep))$subtype,
               "N-term extension")
})

test_that("an unrelated protein is an alternative N-terminus", {
  expect_equal(one(base_pep, "MWWWWGIVGTRSTQXYZ")$subtype, "alt N-terminus")
})

test_that("ORF loss requires an ORF in the reference isoform", {
  expect_equal(one(base_pep, "")$subtype, "ORF loss")
  expect_equal(one("", base_pep)$subtype, "ORF gain")
  # Neither isoform coding is a pre-existing failure, not a consequence of the
  # switch, so it must not be reported as a loss.
  r <- one("", "")
  expect_equal(r$category, "Unclassified")
  expect_equal(r$subtype, "no ORF in either")
})

test_that("non-overlapping transcript spans are excluded as distal", {
  r <- one(base_pep, base_pep,
           main_tx_start = 1000, main_tx_end = 2000,
           alt_tx_start = 50000, alt_tx_end = 60000)
  expect_equal(r$category, "Distal (excluded)")
  expect_equal(r$subtype, "non-overlapping isoforms")
  expect_equal(r$body_overlap, 0)
})

test_that("overlapping spans survive the distal guard", {
  r <- one(base_pep, base_pep, main_utr5_len = 100, alt_utr5_len = 300,
           main_tx_start = 1000, main_tx_end = 9000,
           alt_tx_start = 1050, alt_tx_end = 9000)
  expect_equal(r$category, "5'UTR change")
  expect_gt(r$body_overlap, 0)
})

test_that("uORF gain is computed from the upstream AUG counts", {
  r <- one(base_pep, base_pep, 100, 400, main_n_uATG = 0, alt_n_uATG = 2)
  expect_equal(r$n_uATG_gained, 2)
  expect_equal(r$uorf_gained, 1L)
})

test_that("classify_consequence names every missing required column", {
  expect_error(classify_consequence(data.frame(main_protein = "M")),
               "alt_protein")
})

test_that("consequence_summary orders subtypes canonically, not by count", {
  d <- data.frame(
    category = c("N-terminal/CDS", "5'UTR change", "5'UTR change"),
    subtype = c("alt N-terminus", "shorter", "longer")
  )
  s <- consequence_summary(d)
  expect_equal(as.character(s$subtype), c("longer", "shorter", "alt N-terminus"))
  expect_equal(sum(s$proportion), 100)
})

test_that("consequence_summary drops non-biological categories by default", {
  d <- data.frame(
    category = c("5'UTR change", "Distal (excluded)", "Unclassified"),
    subtype = c("longer", "non-overlapping isoforms", "no ORF in either")
  )
  s <- consequence_summary(d)
  expect_equal(nrow(s), 1L)
  expect_equal(s$n, 1L)
})

test_that("consequence_summary honours a filter expression", {
  d <- data.frame(
    category = c("5'UTR change", "5'UTR change"),
    subtype = c("longer", "shorter"),
    keep = c(1, 0)
  )
  s <- consequence_summary(d, filter = keep == 1)
  expect_equal(nrow(s), 1L)
  expect_equal(as.character(s$subtype), "longer")
})

test_that("the shipped demo table classifies into the expected classes", {
  f <- system.file("extdata", "demo_consequence.tsv", package = "omakase")
  skip_if(!nzchar(f))
  d <- utils::read.delim(f)
  expect_true(all(d$category %in% consequence_levels()$categories))
  expect_true(all(d$subtype %in% consequence_levels()$subtypes))
  expect_gt(nrow(d), 100)
})
