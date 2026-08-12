ex <- function(...) system.file("extdata", ..., package = "omakase")
has_example <- nzchar(ex("samples.tsv"))

test_that("regions parse in every accepted form", {
  r <- parse_region("chr10:27,035,000-27,050,000")
  expect_equal(r$chrom, "chr10")
  expect_equal(r$start, 27035000)
  expect_equal(r$end, 27050000)
  expect_equal(r$strand, "*")

  expect_equal(parse_region("chr10:+:100:200")$strand, "+")
  expect_equal(parse_region("chr1:100..200")$end, 200)
  expect_equal(parse_region("chr1:100-200:-")$strand, "-")
  expect_true(is.na(parse_region("chrM")$start))
})

test_that("BED-style regions are shifted to one-based coordinates", {
  expect_equal(parse_region("chr1:99-200", one_based = FALSE)$start, 100)
})

test_that("an unparseable or inverted region is rejected", {
  expect_error(parse_region("not a region"))
  expect_error(parse_region("chr1:200-100"), "before")
})

test_that("adaptive formatting never prints a non-zero value as zero", {
  expect_equal(format_activity(c(0, 0.004, 0.37, 4.2, 103.01, 4218)),
               c("0", "<0.01", "0.37", "4.2", "103", "4,218"))
})

test_that("percentages drop the decimal only when exact", {
  expect_equal(format_percent(c(25, 12.5)), c("25%", "12.5%"))
})

test_that("a manifest is read with or without a header", {
  skip_if_not(has_example)
  m <- read_manifest(ex("samples.tsv"))
  expect_equal(nrow(m), 6L)
  expect_equal(sort(unique(m$group)),
               c("Endothelial", "Epithelial", "Mesenchymal"))
  expect_true(all(file.exists(m$path)))

  f <- withr::local_tempfile()
  writeLines(c("s1\ta.bam\tCtrl", "s2\tb.bam\tTreat"), f)
  m2 <- read_manifest(f)
  expect_equal(m2$sample, c("s1", "s2"))
  expect_equal(m2$group, c("Ctrl", "Treat"))
})

test_that("a bare BAM path is treated as a one-sample manifest", {
  m <- read_manifest("/tmp/whatever.bam")
  expect_equal(nrow(m), 1L)
  expect_equal(m$sample, "whatever")
})

test_that("a palette file is read in ggsashimi's format", {
  f <- withr::local_tempfile()
  writeLines(c("orange", "cornflowerblue", "#008000"), f)
  expect_equal(read_palette(f), c("orange", "cornflowerblue", "#008000"))
})

test_that("palettes interpolate rather than recycle when short", {
  p <- omakase_palette("omakase", 20)
  expect_length(p, 20)
  expect_equal(anyDuplicated(p), 0L)
})

test_that("STAR junction files are read onto the right coordinates", {
  f <- withr::local_tempfile()
  writeLines(c("chr1\t1001\t1999\t1\t1\t1\t42\t3\t30",
               "chr1\t3001\t3999\t2\t2\t0\t2\t0\t20"), f)
  j <- read_junctions(f)
  expect_equal(nrow(j), 2L)
  expect_equal(j$x0, c(1000, 3000))
  expect_equal(j$x1, c(2000, 4000))
  expect_equal(j$count, c(42, 2))
  expect_equal(j$strand, c("+", "-"))

  expect_equal(nrow(read_junctions(f, min_count = 10)), 1L)
  expect_equal(nrow(read_junctions(f, annotated_only = TRUE)), 1L)
  expect_equal(read_junctions(f, include_multimappers = TRUE)$count[1], 45)
})

test_that("tag BEDs are read as one-based positions", {
  f <- withr::local_tempfile()
  writeLines(c("chr1\t999\t1000\ttag\t3\t+", "chr1\t1999\t2000\ttag\t1\t-"), f)
  t <- read_tag_bed(f)
  expect_equal(t$pos, c(1000, 2000))
  expect_equal(t$count, c(3, 1))
  expect_equal(nrow(read_tag_bed(f, strand = "+")), 1L)
})

test_that("tag tracks bin, scale and extend toward the 3' end", {
  tags <- data.frame(chrom = "chr1", pos = 1000, count = 10, strand = "+",
                     sample = "s1", group = "g")
  x <- sashimi_from_tags(tags, "chr1:900-1400", bin = 50, footprint = 200,
                         normalize = "none")
  tr <- x$tracks
  # Signal starts at the tag's bin and runs 200 bp downstream, not upstream.
  expect_equal(tr$value[tr$pos == 950], 0)
  expect_gt(tr$value[tr$pos == 1000], 0)
  expect_gt(tr$value[tr$pos == 1150], 0)
  expect_equal(tr$value[tr$pos == 1250], 0)
})

test_that("minus-strand tags extend toward lower coordinates", {
  tags <- data.frame(chrom = "chr1", pos = 1400, count = 10, strand = "-",
                     sample = "s1", group = "g")
  x <- sashimi_from_tags(tags, "chr1:1000-1500", bin = 50, footprint = 200,
                         normalize = "none", strand = "-")
  tr <- x$tracks
  expect_gt(tr$value[tr$pos == 1250], 0)
  expect_equal(tr$value[tr$pos == 1050], 0)
})

test_that("rMATS layouts cover every event type", {
  for (et in event_types()) {
    expect_true(!is.null(omakase:::RMATS_LAYOUT[[et]]), info = et)
  }
})

test_that("plotting returns objects rather than writing files", {
  x <- sashimi_data(
    loci = data.frame(locus_id = "a", gene_name = "DEMO", chrom = "chr1",
                      strand = "+", win_lo = 1000, win_hi = 2000),
    tracks = data.frame(locus_id = "a", group = c(rep("x", 51), rep("y", 51)),
                        pos = rep(seq(1000, 2000, 20), 2),
                        value = c(1:51, 51:1)),
    junctions = data.frame(locus_id = "a", group = c("x", "y"), x0 = 1200,
                           x1 = 1800, count = c(40, 5)),
    models = data.frame(locus_id = "a", tx_id = "t1", role = "main",
                        start = c(1100, 1700), end = c(1300, 1900))
  )
  p <- plot_sashimi(x)
  expect_s3_class(p, "patchwork")
  expect_equal(attr(p, "omakase_panels"), 2L)
  expect_s3_class(sashimi_track(x), "ggplot")
  expect_s3_class(sashimi_annotation(x), "ggplot")
})

test_that("every arc shape and shrink method draws without error", {
  x <- sashimi_data(
    loci = data.frame(locus_id = "a", gene_name = "D", chrom = "chr1",
                      strand = "+", win_lo = 1000, win_hi = 9000),
    tracks = data.frame(locus_id = "a", group = "g",
                        pos = seq(1000, 9000, 100), value = 1),
    junctions = data.frame(locus_id = "a", group = "g", x0 = 1500, x1 = 8500,
                           count = 10),
    models = data.frame(locus_id = "a", tx_id = "t1", role = "main",
                        start = c(1100, 8000), end = c(1500, 8800))
  )
  for (s in arc_shapes()) {
    expect_s3_class(plot_sashimi(x, arc_shape = s), "patchwork")
  }
  for (m in shrink_methods()) {
    expect_s3_class(plot_sashimi(x, shrink = TRUE, shrink_method = m),
                    "patchwork")
  }
})

test_that("every consequence style draws without error", {
  f <- ex("demo_consequence.tsv")
  skip_if(!nzchar(f))
  d <- utils::read.delim(f)
  for (s in c("donut", "bar", "lollipop", "stacked")) {
    expect_s3_class(plot_consequence(d, style = s), "ggplot")
  }
})

test_that("reading a real BAM gives coverage and junctions", {
  skip_if_not(has_example)
  skip_if_not_installed("GenomicAlignments")
  x <- sashimi_from_bam(ex("samples.tsv"), "chr10:27040584-27048100",
                        annotation = ex("annotation.gtf"), min_count = 10)
  expect_gt(nrow(x$tracks), 0)
  expect_gt(nrow(x$junctions), 0)
  expect_gt(nrow(x$models), 0)
  # The gene name and strand come from the annotation, not the region string.
  expect_equal(x$loci$gene_name, "ABI1")
  expect_equal(x$loci$strand, "-")
  # Junctions are contained in the window by default.
  expect_true(all(x$junctions$x0 >= 27040584))
  expect_true(all(x$junctions$x1 <= 27048100))
  expect_s3_class(plot_sashimi(x, aggregate = "mean"), "patchwork")
})

test_that("saving a figure writes a file", {
  skip_on_cran()
  x <- sashimi_data(
    loci = data.frame(locus_id = "a", gene_name = "A", chrom = "chr1",
                      strand = "+", win_lo = 1, win_hi = 10),
    tracks = data.frame(locus_id = "a", group = "g", pos = 1:10, value = 1:10)
  )
  f <- withr::local_tempfile(fileext = ".pdf")
  save_sashimi(plot_sashimi(x), f)
  expect_true(file.exists(f))
})

test_that("overlay groups panels and collapses their shared junctions", {
  sd <- sashimi_data(
    loci = data.frame(locus_id = "a", gene_name = "A", chrom = "chr1",
                      strand = "+", win_lo = 1, win_hi = 100),
    tracks = data.frame(locus_id = "a", group = rep(c("x", "y", "z"), each = 10),
                        pos = rep(1:10 * 10, 3), value = 1),
    junctions = data.frame(locus_id = "a", group = c("x", "y", "z"),
                           x0 = 20, x1 = 80, count = c(10, 20, 60))
  )
  # three panels by default, two once y and z share one
  expect_equal(attr(plot_sashimi(sd), "omakase_panels"), 3L)
  p <- plot_sashimi(sd, overlay = list(x = "x", yz = c("y", "z")))
  expect_equal(attr(p, "omakase_panels"), 2L)

  # the shared junction becomes a single arc at the combined count
  g <- ggplot2::ggplot_build(sashimi_track(sd, group = c("y", "z")))
  lab <- g$data[[which(vapply(g$data, function(d) "label" %in% names(d),
                              logical(1)))[1]]]
  expect_equal(nrow(lab), 1L)
  expect_equal(lab$label[1], "40")           # mean of 20 and 60
  expect_equal(
    ggplot2::ggplot_build(
      sashimi_track(sd, group = c("y", "z"), overlay_junction_fun = "sum")
    )$data[[which(vapply(
      ggplot2::ggplot_build(sashimi_track(sd, group = c("y", "z"),
                                          overlay_junction_fun = "sum"))$data,
      function(d) "label" %in% names(d), logical(1)))[1]]]$label[1],
    "80"
  )
})

test_that("overlay by column name reads the tracks table", {
  sd <- sashimi_data(
    loci = data.frame(locus_id = "a", gene_name = "A", chrom = "chr1",
                      strand = "+", win_lo = 1, win_hi = 10),
    tracks = data.frame(locus_id = "a", group = c("s1", "s2"),
                        pos = 1, value = 1, batch = c("b1", "b1"))
  )
  expect_equal(attr(plot_sashimi(sd, overlay = "batch"), "omakase_panels"), 1L)
  expect_error(plot_sashimi(sd, overlay = "nope"), "nope")
  expect_error(plot_sashimi(sd, overlay = list(a = "missing")), "missing")
})

test_that("exon_scale shrinks exons and stays invertible", {
  m <- intron_map(data.frame(start = 2000, end = 8000), 1000, 10000,
                  exon_scale = 2)
  x <- c(1000, 1500, 9000, 10000)
  expect_equal(expand_coords(compress_coords(x, m), m), x)
  # the 1000 bp exonic run before the intron is drawn at half length
  expect_equal(diff(compress_coords(c(1000, 2000), m)), 500)

  m2 <- intron_map(NULL, 0, 1000, exon_scale = 4)
  expect_equal(compress_coords(1000, m2), 250)
  expect_equal(expand_coords(compress_coords(750, m2), m2), 750)
})

test_that("presets are bundles that explicit arguments override", {
  expect_true(all(c("default", "tss", "junction", "minimal", "igv") %in%
                    presets()))
  sd <- sashimi_data(
    loci = data.frame(locus_id = "a", gene_name = "A", chrom = "1",
                      strand = "+", win_lo = 1, win_hi = 100),
    tracks = data.frame(locus_id = "a", group = "g", pos = 1:10 * 10, value = 1)
  )
  for (pr in presets()) expect_s3_class(plot_sashimi(sd, preset = pr), "patchwork")
  # the tss preset prints a UCSC-style contig name
  ann <- ggplot2::ggplot_build(
    sashimi_annotation(sd, chrom_style = "ucsc")
  )
  expect_true(any(vapply(ann$data, function(d)
    any(grepl("^chr1 ", as.character(d$label %||% ""))), logical(1))))
})

test_that("named arc heights are keyed by role, not by position", {
  j <- data.frame(locus_id = "a", group = "g", x0 = c(10, 30), x1 = 90,
                  count = c(5, 5), role = c("main", "ATSS"),
                  label = NA_character_)
  b <- omakase:::build_arcs(j, ymax = 10, shape = "sine",
                            height_rule = "constant",
                            height_frac = c(main = 0.8, ATSS = 1.2),
                            n = 21, side = "above")
  expect_equal(b$junctions$.h[b$junctions$role == "main"], 8)
  expect_equal(b$junctions$.h[b$junctions$role == "ATSS"], 12)
})

test_that("BED12 blocks expand into CDS and UTR exons", {
  f <- withr::local_tempfile()
  writeLines(paste("chr1", 999, 5000, "tx1", 0, "+", 1199, 4000,
                   "31,120,180", 2, "500,800,", "0,3200,", sep = "\t"), f)
  d <- read_bed(f)
  expect_true(all(c("CDS", "UTR") %in% d$feature))
  # first block 1000-1499: 1000-1199 is UTR, 1200-1499 is CDS
  expect_equal(d$start[1], 1000)
  expect_equal(d$end[1], 1199)
  expect_equal(d$feature[1], "UTR")
  expect_equal(d$color[1], "#1F78B4")
})

test_that("junctions round-trip through a BED file", {
  sd <- sashimi_data(
    loci = data.frame(locus_id = "a", gene_name = "A", chrom = "chr1",
                      strand = "+", win_lo = 1, win_hi = 5000),
    junctions = data.frame(locus_id = "a", group = "g", x0 = 1000, x1 = 2000,
                           count = 42)
  )
  f <- withr::local_tempfile(fileext = ".bed")
  write_junctions(sd, f)
  back <- read_junctions(f, format = "bed")
  expect_equal(back$x0, 1000)
  expect_equal(back$x1, 2000)
  expect_equal(back$count, 42)
})

test_that("genome tracks build from a sashimi object and from parts", {
  sd <- sashimi_data(
    loci = data.frame(locus_id = "a", gene_name = "A", chrom = "chr1",
                      strand = "-", win_lo = 1000, win_hi = 5000),
    tracks = data.frame(locus_id = "a", group = "s1",
                        pos = seq(1000, 5000, 100), value = 1),
    models = data.frame(locus_id = "a", tx_id = "t1", feature = "CDS",
                        start = c(1200, 3000), end = c(1600, 4200))
  )
  expect_s3_class(plot_tracks(sd), "patchwork")
  expect_s3_class(
    plot_tracks(tracks = list(track_models(sd), track_axis()),
                region = "chr1:1000-5000"),
    "patchwork"
  )
  expect_s3_class(track_models(sd), "omakase_track")
  expect_error(plot_tracks(region = "chr1:1-10"), "tracks")
})
