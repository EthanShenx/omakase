iv <- data.frame(start = c(1200, 4000), end = c(3000, 9000))

test_that("phi and its inverse round-trip exactly", {
  m <- intron_map(iv, 1000, 10000)
  x <- c(1000, 1100, 1500, 3500, 9500, 10000, 1200, 3000)
  expect_equal(expand_coords(compress_coords(x, m), m), x)
})

test_that("phi is strictly increasing across the window", {
  m <- intron_map(iv, 1000, 10000)
  y <- compress_coords(seq(1000, 10000, by = 5), m)
  expect_true(all(diff(y) > 0))
})

test_that("compression shortens the window and leaves exons at scale 1", {
  m <- intron_map(iv, 1000, 10000)
  plotted <- diff(range(compress_coords(c(1000, 10000), m)))
  expect_lt(plotted, 9000)

  # 1000-1200 is exonic, so its drawn length is unchanged.
  seg <- diff(compress_coords(c(1000, 1200), m))
  expect_equal(seg, 200)
})

test_that("the power rule draws an intron at L^gamma", {
  m <- intron_map(data.frame(start = 1000, end = 11000), 0, 12000,
                  method = "power", gamma = 0.7)
  drawn <- diff(compress_coords(c(1000, 11000), m))
  expect_equal(drawn, 10000^0.7, tolerance = 1e-6)
})

test_that("the scale rule divides the intron length", {
  m <- intron_map(data.frame(start = 100, end = 1100), 0, 2000,
                  method = "scale", scale = 5)
  expect_equal(diff(compress_coords(c(100, 1100), m)), 200)
})

test_that("the fixed rule caps the drawn intron", {
  m <- intron_map(data.frame(start = 100, end = 5100), 0, 6000,
                  method = "fixed", cap = 250)
  expect_equal(diff(compress_coords(c(100, 5100), m)), 250)
})

test_that("method 'none' is the identity", {
  m <- intron_map(iv, 1000, 10000, method = "none")
  x <- c(1000, 5000, 10000)
  expect_equal(compress_coords(x, m), x)
  expect_equal(expand_coords(x, m), x)
})

test_that("short introns are left alone", {
  m <- intron_map(data.frame(start = 1000, end = 1050), 0, 2000,
                  min_intron = 100)
  expect_equal(compress_coords(1500, m), 1500)
})

test_that("overlapping introns are merged so the map stays monotone", {
  m <- intron_map(data.frame(start = c(1000, 1500), end = c(2000, 3000)),
                  0, 4000)
  y <- compress_coords(seq(0, 4000, by = 10), m)
  expect_true(all(diff(y) > 0))
})

test_that("coordinates beyond the window extrapolate at scale 1", {
  m <- intron_map(iv, 1000, 10000)
  lo <- compress_coords(c(900, 1000), m)
  expect_equal(diff(lo), 100)
  hi <- compress_coords(c(10000, 10100), m)
  expect_equal(diff(hi), 100)
})

test_that("an empty or NULL intron set gives the identity map", {
  expect_equal(compress_coords(500, intron_map(NULL, 0, 1000)), 500)
  expect_equal(compress_coords(500, intron_map(data.frame(), 0, 1000)), 500)
})

test_that("intron_trans reports breaks in genome coordinates", {
  m <- intron_map(iv, 1000, 10000)
  b <- intron_trans(m)$breaks(c(1000, 10000))
  expect_true(all(b >= 1000 & b <= 10000))
  expect_true(length(b) >= 2)
})

test_that("introns are derived from exon models", {
  models <- data.frame(
    locus_id = "a", tx_id = "t1",
    start = c(100, 500, 900), end = c(200, 600, 1000)
  )
  iv2 <- omakase:::introns_from_models(models)
  expect_equal(nrow(iv2), 2L)
  expect_equal(iv2$start, c(201, 601))
  expect_equal(iv2$end, c(499, 899))
})
