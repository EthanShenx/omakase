# Read coverage and junctions from alignment files

Reads one or more BAM/CRAM files over a region and returns a
[`sashimi_data()`](https://EthanShenx.github.io/omakase/reference/sashimi_data.md)
object holding binned coverage and spliced-junction counts. This is the
main entry point when starting from alignments.

## Usage

``` r
sashimi_from_bam(
  bam,
  region,
  annotation = NULL,
  bin = NULL,
  n_bins = 800,
  group_col = NULL,
  label_col = NULL,
  strand = "none",
  keep_strand = "both",
  min_count = 1,
  min_mapq = 0,
  flags = NULL,
  per_sample = TRUE,
  gene = NULL,
  flank = 0,
  junction_overlap = c("within", "any")
)
```

## Arguments

- bam:

  A path to a BAM/CRAM/SAM file, a vector of paths, or a manifest read
  by
  [`read_manifest()`](https://EthanShenx.github.io/omakase/reference/read_manifest.md)
  (a path to a TSV works directly). A `.sam` file is converted to a
  sorted, indexed BAM in the session's temporary directory the first
  time it is read, which is what `rmats2sashimiplot` does for its
  `--s1`/`--s2` inputs.

- region:

  A region string such as `"chr10:27035000-27050000"`, a
  [`parse_region()`](https://EthanShenx.github.io/omakase/reference/parse_region.md)
  object, or a data frame of regions with columns `chrom`, `start`,
  `end` and optionally `name`.

- annotation:

  Optional GTF/GFF3 path, or the result of
  [`read_annotation()`](https://EthanShenx.github.io/omakase/reference/read_annotation.md),
  used to draw transcript models under the tracks.

- bin:

  Bin width in base pairs. `1` gives per-base coverage, which is exact
  but slow to draw over a wide window; the default of `NULL` picks a bin
  that gives roughly `n_bins` points.

- n_bins:

  Target number of bins when `bin` is `NULL`.

- group_col, label_col:

  Manifest columns used for grouping and panel labels; see
  [`read_manifest()`](https://EthanShenx.github.io/omakase/reference/read_manifest.md).

- strand:

  Strand specificity, see
  [`strand_modes()`](https://EthanShenx.github.io/omakase/reference/strand_modes.md).

- keep_strand:

  Which strand to keep once reads are assigned: `"both"`, `"+"` or
  `"-"`. Only meaningful when `strand` is not `"none"`.

- min_count:

  Drop junctions with fewer supporting reads.

- min_mapq:

  Ignore alignments below this mapping quality.

- flags:

  A
  [`Rsamtools::scanBamFlag()`](https://rdrr.io/pkg/Rsamtools/man/ScanBamParam-class.html)
  value. The default drops unmapped reads, secondary alignments,
  duplicates and QC failures.

- per_sample:

  Keep one track per sample. When `FALSE`, samples sharing a group are
  summed as they are read, which is cheaper in memory.

- gene:

  Restrict `annotation` to this gene name.

- flank:

  Widen the region by this many base pairs on each side.

- junction_overlap:

  Which junctions to keep: `"within"` (the default) keeps only junctions
  whose donor and acceptor both fall inside the region, `"any"` also
  keeps those with one end beyond it, whose arcs then run off the edge
  of the panel.

## Value

A `sashimi_data` object.

## Details

Coverage is the depth of reference-consuming CIGAR operations (`M`, `D`,
`=`, `X`), binned to `bin` base pairs by taking the mean depth within
each bin. Junctions are `N` operations, counted once per read and keyed
by their donor and acceptor coordinates; a junction supported by fewer
than `min_count` reads is dropped.

For a stranded library, set `strand` so that reads are assigned to the
correct strand, and `keep_strand` to draw only one of them. For an
unstranded library leave `strand = "none"` and everything is pooled.

## Examples

``` r
bams <- system.file("extdata", "samples.tsv", package = "omakase")
gtf <- system.file("extdata", "annotation.gtf", package = "omakase")
if (nzchar(bams)) {
  sd <- sashimi_from_bam(bams, "chr10:27040584-27048100", annotation = gtf,
                         min_count = 10)
  plot_sashimi(sd, aggregate = "mean")
}

```
