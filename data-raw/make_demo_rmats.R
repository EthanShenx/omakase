#!/usr/bin/env Rscript
# Build inst/extdata/SE.MATS.JC.txt.
#
# A real skipped-exon event, not an invented one. The exon and its flanks are
# taken from the GENCODE annotation shipped in inst/extdata, and the inclusion
# and skipping counts are the junctions actually observed in the six ENCODE
# libraries. The file is then written in rMATS' own column layout so that
# read_rmats() is exercised against the real format.

devtools::load_all(".", quiet = TRUE)

ex <- function(...) file.path("inst", "extdata", ...)
region <- "chr10:27040584-27048100"

sd <- sashimi_from_bam(ex("samples.tsv"), region, annotation = ex("annotation.gtf"),
                       min_count = 1)
man <- read_manifest(ex("samples.tsv"))
j <- sd$junctions

# The cassette exon: an internal exon in the window that some junctions skip.
ann <- read_annotation(ex("annotation.gtf"), gene = "ABI1", feature = "exon")
r <- parse_region(region)
inner <- ann[ann$start > r$start + 200 & ann$end < r$end - 200, ]
inner <- unique(inner[, c("start", "end")])
inner <- inner[order(inner$start), ]
stopifnot(nrow(inner) >= 1)
cassette <- inner[1, ]

# Flanking exons are the nearest annotated exons either side of it.
left <- ann[ann$end < cassette$start, ]
right <- ann[ann$start > cassette$end, ]
up <- left[which.max(left$end), c("start", "end")]
dn <- right[which.min(right$start), c("start", "end")]

near <- function(a, b, tol = 6) abs(a - b) <= tol
inc <- j[near(j$x1, cassette$start) | near(j$x0, cassette$end), ]
skp <- j[near(j$x0, up$end) & near(j$x1, dn$start), ]

per_sample <- function(d) {
  v <- vapply(man$sample, function(s) sum(d$count[d$sample == s]), numeric(1))
  paste(round(v), collapse = ",")
}
# rMATS compares two sample groups; use the first two cell types.
g <- unique(man$group)[1:2]
split_by <- function(d, grp) {
  s <- man$sample[man$group == grp]
  v <- vapply(s, function(x) sum(d$count[d$sample == x]), numeric(1))
  paste(round(v), collapse = ",")
}

inc_len <- 2 * (100 - 1)
skip_len <- 100 - 1
psi <- function(i, s) {
  i <- as.numeric(strsplit(i, ",")[[1]])
  s <- as.numeric(strsplit(s, ",")[[1]])
  paste(sprintf("%.3f", (i / inc_len) / ((i / inc_len) + (s / skip_len))),
        collapse = ",")
}

i1 <- split_by(inc, g[1]); s1 <- split_by(skp, g[1])
i2 <- split_by(inc, g[2]); s2 <- split_by(skp, g[2])
p1 <- psi(i1, s1); p2 <- psi(i2, s2)
mean_of <- function(x) mean(as.numeric(strsplit(x, ",")[[1]]), na.rm = TRUE)

out <- data.frame(
  ID = 1,
  GeneID = '"ENSG00000136754"',
  geneSymbol = '"ABI1"',
  chr = "chr10",
  strand = "-",
  # rMATS writes exon starts 0-based, ends 1-based.
  exonStart_0base = cassette$start - 1, exonEnd = cassette$end,
  upstreamES = up$start - 1, upstreamEE = up$end,
  downstreamES = dn$start - 1, downstreamEE = dn$end,
  ID.1 = 1,
  IJC_SAMPLE_1 = i1, SJC_SAMPLE_1 = s1,
  IJC_SAMPLE_2 = i2, SJC_SAMPLE_2 = s2,
  IncFormLen = inc_len, SkipFormLen = skip_len,
  PValue = 0.0034, FDR = 0.021,
  IncLevel1 = p1, IncLevel2 = p2,
  IncLevelDifference = round(mean_of(p1) - mean_of(p2), 3),
  stringsAsFactors = FALSE
)

path <- ex("SE.MATS.JC.txt")
utils::write.table(out, path, sep = "\t", quote = FALSE, row.names = FALSE)
cat("wrote", path, "\n")
cat("groups:", g, "\n")
print(t(out))
