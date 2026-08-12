mk <- function(values, samples, groups = NULL) {
  n <- length(values) / length(samples)
  sashimi_data(
    loci = data.frame(locus_id = "a", gene_name = "A", chrom = "chr1",
                      strand = "+", win_lo = 1, win_hi = n),
    tracks = data.frame(
      locus_id = "a",
      group = rep(groups %||% samples, each = n),
      sample = rep(samples, each = n),
      pos = rep(seq_len(n), length(samples)),
      value = values
    )
  )
}
`%||%` <- function(x, y) if (is.null(x)) y else x

test_that("cpm makes every sample sum to a million", {
  x <- mk(c(rep(1, 10), rep(3, 10)), c("s1", "s2"))
  y <- suppressMessages(normalize_tracks(x, "cpm"))
  tot <- tapply(y$tracks$value, y$tracks$sample, sum)
  expect_equal(as.numeric(tot), c(1e6, 1e6))
})

test_that("explicit library sizes are used when supplied", {
  x <- mk(rep(1, 20), c("s1", "s2"))
  y <- normalize_tracks(x, "cpm", library_sizes = c(s1 = 1e6, s2 = 2e6))
  v <- tapply(y$tracks$value, y$tracks$sample, sum)
  expect_equal(as.numeric(v["s1"]), 10)
  expect_equal(as.numeric(v["s2"]), 5)
})

test_that("a missing library size is an error, not a silent default", {
  x <- mk(rep(1, 20), c("s1", "s2"))
  expect_error(normalize_tracks(x, "cpm", library_sizes = c(s1 = 1e6)), "s2")
})

test_that("max scaling puts every sample's peak at one", {
  x <- mk(c(1:10, seq(2, 20, by = 2)), c("s1", "s2"))
  y <- normalize_tracks(x, "max")
  expect_equal(as.numeric(tapply(y$tracks$value, y$tracks$sample, max)),
               c(1, 1))
})

test_that("DESeq2 size factors recover a known depth ratio", {
  # s2 is exactly twice s1 everywhere, so the factors must differ two-fold.
  set.seed(1)
  base <- stats::rpois(200, 40) + 1
  x <- mk(c(base, base * 2), c("s1", "s2"))
  y <- normalize_tracks(x, "size_factor")
  f <- y$meta$normalize_factors
  expect_equal(unname(f["s2"] / f["s1"]), 2, tolerance = 0.05)
})

test_that("normalisation rescales junction counts alongside the coverage", {
  x <- mk(rep(1, 20), c("s1", "s2"))
  x$junctions <- omakase:::normalise_slot(
    data.frame(locus_id = "a", group = c("s1", "s2"), sample = c("s1", "s2"),
               x0 = 2, x1 = 8, count = c(100, 100)),
    "junctions"
  )
  y <- normalize_tracks(x, "manual", factors = c(s1 = 1, s2 = 4))
  expect_equal(y$junctions$count, c(100, 25))
})

test_that("method 'none' leaves the object untouched", {
  x <- mk(rep(3, 10), "s1")
  expect_equal(normalize_tracks(x, "none")$tracks$value, x$tracks$value)
})

test_that("aggregation collapses replicates to one track per group", {
  x <- mk(c(1:5, 3:7), c("s1", "s2"), groups = c("g", "g"))
  y <- aggregate_tracks(x, "mean")
  expect_equal(nrow(y$tracks), 5L)
  expect_equal(y$tracks$value, c(2, 3, 4, 5, 6))
  expect_true(all(is.na(y$tracks$sample)))
})

test_that("junctions can be summed while coverage is averaged", {
  x <- mk(c(1:5, 1:5), c("s1", "s2"), groups = c("g", "g"))
  x$junctions <- omakase:::normalise_slot(
    data.frame(locus_id = "a", group = "g", x0 = 2, x1 = 4,
               count = c(10, 30)),
    "junctions"
  )
  y <- aggregate_tracks(x, "mean", junction_fun = "sum")
  expect_equal(y$junctions$count, 40)
  expect_equal(y$tracks$value, as.numeric(1:5))
})

test_that("PSI is the reference share of the total", {
  j <- data.frame(
    locus_id = "a", group = c("early", "early", "late", "late"),
    role = c("main", "alt", "main", "alt"), count = c(90, 10, 20, 80)
  )
  p <- compute_psi(j)
  expect_equal(p$psi[p$group == "early"], 0.9)
  expect_equal(p$psi[p$group == "late"], 0.2)
})

test_that("PSI is NA rather than zero when nothing was observed", {
  j <- data.frame(locus_id = "a", group = "g", role = c("main", "alt"),
                  count = c(0, 0))
  expect_true(is.na(compute_psi(j)$psi))
})

test_that("the rMATS form corrects for effective length", {
  j <- data.frame(locus_id = "a", group = "g", role = c("main", "alt"),
                  count = c(100, 100))
  p <- compute_psi(j, method = "rmats", len_inc = 200, len_skip = 100)
  # 100/200 vs 100/100 -> 0.5/(0.5+1) = 1/3
  expect_equal(p$psi, 1 / 3)
})

test_that("compute_psi needs a role column", {
  j <- data.frame(locus_id = "a", group = "g", count = 1)
  expect_error(compute_psi(j), "role")
})

test_that("delta_psi subtracts the second group from the first", {
  p <- data.frame(locus_id = c("a", "a"), group = c("early", "late"),
                  psi = c(0.9, 0.2))
  d <- delta_psi(p, "early", "late")
  expect_equal(d$dpsi, 0.7)
})
