demo_loci <- function(ids = "a") {
  data.frame(locus_id = ids, gene_name = toupper(ids), chrom = "chr1",
             strand = "+", win_lo = 1000, win_hi = 2000,
             stringsAsFactors = FALSE)
}

test_that("an empty object has correctly typed empty slots", {
  x <- sashimi_data()
  expect_s3_class(x, "sashimi_data")
  for (s in c("loci", "tracks", "junctions", "models", "psi", "features")) {
    expect_equal(nrow(x[[s]]), 0L, info = s)
  }
  expect_type(x$loci$locus_id, "character")
  expect_type(x$tracks$value, "double")
})

test_that("optional columns are filled with their defaults", {
  x <- sashimi_data(
    loci = demo_loci(),
    tracks = data.frame(locus_id = "a", group = "g", pos = 1:3, value = 1:3)
  )
  expect_true("sample" %in% names(x$tracks))
  expect_true(all(is.na(x$tracks$sample)))
  expect_equal(unique(x$tracks$strand), "*")
})

test_that("a reversed junction is normalised so x0 < x1", {
  x <- sashimi_data(
    loci = demo_loci(),
    junctions = data.frame(locus_id = "a", group = "g", x0 = 1800, x1 = 1200,
                           count = 5)
  )
  expect_equal(x$junctions$x0, 1200)
  expect_equal(x$junctions$x1, 1800)
})

test_that("missing required columns are all reported at once", {
  expect_error(
    sashimi_data(loci = data.frame(locus_id = "a")),
    "gene_name"
  )
})

test_that("validation catches a locus_id no locus declares", {
  x <- sashimi_data(
    loci = demo_loci("a"),
    tracks = data.frame(locus_id = "b", group = "g", pos = 1, value = 1)
  )
  expect_error(validate_sashimi_data(x), "locus_id")
  expect_warning(validate_sashimi_data(x, strict = FALSE), "locus_id")
})

test_that("validation catches an inverted window", {
  x <- sashimi_data(loci = data.frame(
    locus_id = "a", gene_name = "A", chrom = "chr1", strand = "+",
    win_lo = 2000, win_hi = 1000
  ))
  expect_error(validate_sashimi_data(x), "window")
})

test_that("subsetting keeps only the selected locus across every slot", {
  x <- sashimi_data(
    loci = demo_loci(c("a", "b")),
    tracks = data.frame(locus_id = c("a", "b"), group = "g", pos = 1,
                        value = 1),
    junctions = data.frame(locus_id = c("a", "b"), group = "g", x0 = 1,
                           x1 = 2, count = 1)
  )
  y <- x["a"]
  expect_equal(loci(y), "a")
  expect_equal(nrow(y$tracks), 1L)
  expect_equal(nrow(y$junctions), 1L)
})

test_that("subsetting also accepts a gene name", {
  x <- sashimi_data(loci = demo_loci(c("a", "b")))
  expect_equal(loci(x["A"]), "a")
})

test_that("combining objects unions the loci and drops duplicates", {
  a <- sashimi_data(loci = demo_loci("a"))
  b <- sashimi_data(loci = demo_loci("b"))
  ab <- combine_sashimi(a, b, a)
  expect_equal(sort(loci(ab)), c("a", "b"))
})

test_that("write then read round-trips through a directory", {
  d <- withr::local_tempdir()
  x <- sashimi_data(
    loci = demo_loci(),
    tracks = data.frame(locus_id = "a", group = "g", pos = 1:5,
                        value = c(1, 2, 3, 2, 1)),
    meta = list(source = "test", bin = 25)
  )
  write_sashimi_data(x, d, prefix = "demo")
  y <- read_sashimi_dir(d, prefix = "demo")
  expect_equal(loci(y), "a")
  expect_equal(y$tracks$value, x$tracks$value)
  expect_equal(y$meta$bin, "25")
})

test_that("as_sashimi_data accepts a list of slot tables", {
  x <- as_sashimi_data(list(loci = demo_loci()))
  expect_s3_class(x, "sashimi_data")
  expect_equal(loci(x), "a")
})

test_that("sashimi_from_tables renames columns and derives locus_id", {
  x <- sashimi_from_tables(
    loci = data.frame(gene_name = "DEMO", chrom = "chr1", strand = "+",
                      win_lo = 1, win_hi = 100),
    tracks = data.frame(gene_name = "DEMO", stage = "early", pos = 1:10,
                        rpm = 1:10),
    rename = list(tracks = c(group = "stage", value = "rpm")),
    locus_col = "gene_name"
  )
  expect_equal(loci(x), "DEMO")
  expect_equal(unique(x$tracks$group), "early")
  expect_equal(x$tracks$value, as.numeric(1:10))
})

test_that("a column named like a prefix of another is not mistaken for it", {
  # `$` partial-matches, so a table carrying `n_uATG_gained` but no `n` must
  # still be treated as unsummarised.
  d <- data.frame(category = "5'UTR change", subtype = "longer",
                  n_uATG_gained = 2)
  s <- consequence_summary(d)
  expect_equal(s$n, 1L)
})
