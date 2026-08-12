test_that("every arc shape starts and ends on the baseline", {
  for (s in arc_shapes()) {
    p <- arc_path(100, 200, h = 10, shape = s)
    expect_equal(p$y[1], 0, info = s)
    expect_equal(p$y[nrow(p)], 0, info = s)
    expect_equal(p$x[1], 100, info = s)
    expect_equal(p$x[nrow(p)], 200, info = s)
  }
})

test_that("every arc shape peaks at exactly h", {
  for (s in arc_shapes()) {
    p <- arc_path(0, 1, h = 5, shape = s)
    expect_equal(max(p$y), 5, tolerance = 1e-9, info = s)
  }
})

test_that("arcs never dip below the baseline", {
  for (s in arc_shapes()) {
    p <- arc_path(0, 1, h = 3, shape = s)
    expect_true(all(p$y >= -1e-12), info = s)
  }
})

test_that("a negative height mirrors the arc below the axis", {
  up <- arc_path(0, 10, h = 4, shape = "sine")
  dn <- arc_path(0, 10, h = -4, shape = "sine")
  expect_equal(up$y, -dn$y)
  expect_equal(min(dn$y), -4, tolerance = 1e-9)
})

test_that("the point count is forced odd so the apex is sampled", {
  # An even count straddles t = 0.5 and would miss the peak.
  p <- arc_path(0, 1, h = 1, shape = "sine", n = 100)
  expect_true(nrow(p) %% 2 == 1)
  expect_equal(max(p$y), 1, tolerance = 1e-12)
})

test_that("arc_path rejects unknown shapes", {
  expect_error(arc_path(0, 1, 1, shape = "spiral"))
})

test_that("non-finite input yields an empty path rather than an error", {
  expect_equal(nrow(arc_path(NA, 10, 5)), 0L)
  expect_equal(nrow(arc_path(0, 10, NA)), 0L)
})

test_that("arc widths grow with count and stay inside the clamp", {
  w <- arc_widths(c(1, 10, 100, 1000), rule = "log", range = c(0.2, 1.5))
  expect_true(all(diff(w) > 0))
  expect_true(all(w >= 0.2 & w <= 1.5))
})

test_that("constant width ignores the counts", {
  expect_equal(arc_widths(c(1, 999), rule = "constant", w0 = 0.4),
               c(0.4, 0.4))
})

test_that("arc heights scale monotonically with the chosen quantity", {
  h <- arc_heights(c(1, 10, 100), rule = "log", min_h = 1, max_h = 5)
  expect_equal(h[1], 1)
  expect_equal(h[3], 5)
  expect_true(all(diff(h) > 0))

  hs <- arc_heights(c(1, 1, 1), rule = "span", min_h = 2, max_h = 6,
                    span = c(10, 100, 1000))
  expect_equal(hs[1], 2)
  expect_equal(hs[3], 6)
})

test_that("identical values collapse to the mid height rather than the floor", {
  h <- arc_heights(c(7, 7, 7), rule = "log", min_h = 1, max_h = 5)
  expect_equal(unique(h), 3)
})

test_that("span scaling requires spans", {
  expect_error(arc_heights(1:3, rule = "span"), "span")
})

test_that("build_arcs draws smooth curves, not three-point triangles", {
  # Regression: a local variable once shadowed the point-count argument, so
  # every arc collapsed to its two ends and its apex.
  j <- data.frame(locus_id = "a", group = "g", x0 = 10, x1 = 90,
                  count = 100, role = NA_character_, label = NA_character_)
  b <- omakase:::build_arcs(j, ymax = 1, shape = "sine", height_rule = "auto",
                            height_frac = c(0.8, 1.2), n = 121, side = "above")
  expect_equal(nrow(b$paths), 121L)
  # A sine is well above the straight line joining its endpoints at t = 0.25.
  expect_equal(b$paths$y[31], 0.8 * sin(pi * 0.25), tolerance = 1e-9)
})

test_that("build_arcs staggers arcs that share an anchor", {
  j <- data.frame(locus_id = "a", group = "g", x0 = c(10, 30), x1 = 90,
                  count = c(5, 5), role = c("main", "alt"),
                  label = NA_character_)
  b <- omakase:::build_arcs(j, ymax = 10, shape = "sine",
                            height_rule = "constant",
                            height_frac = c(0.8, 1.2), n = 121, side = "above")
  expect_equal(sort(b$heights), c(8, 12))
})

test_that("arc labels survive the internal re-sort", {
  # build_arcs sorts its input, so labels must travel with the rows.
  j <- data.frame(locus_id = "a", group = "g", x0 = c(50, 10), x1 = c(90, 40),
                  count = c(1, 2), role = NA_character_,
                  label = c("second", "first"))
  b <- omakase:::build_arcs(j, ymax = 1, shape = "sine",
                            height_rule = "constant",
                            height_frac = 1, n = 21, side = "above")
  expect_equal(b$labels$label, c("first", "second"))
  expect_equal(b$labels$count, c(2, 1))
})
