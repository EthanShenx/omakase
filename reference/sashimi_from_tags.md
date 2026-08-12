# Build a sashimi data object from 5'-tag data

Turns single-base 5'-end tags (CAGE, STRT, CamoTSS) into the binned
tracks a sashimi panel draws, and optionally into arcs that run from
each start site to a shared downstream anchor.

## Usage

``` r
sashimi_from_tags(
  tags,
  region,
  annotation = NULL,
  bin = 25,
  footprint = 250,
  normalize = c("tpm", "none"),
  library_sizes = NULL,
  aggregate = c("mean", "median", "sum", "none"),
  strand = "gene",
  activity = NULL,
  group_col = NULL,
  label_col = NULL,
  gene = NULL
)
```

## Arguments

- tags:

  A manifest read by
  [`read_manifest()`](https://EthanShenx.github.io/omakase/reference/read_manifest.md),
  a path to one, a vector of BED paths, or a data frame of tags with
  `chrom`, `pos`, `count`, `strand`.

- region:

  A region string, or a data frame of regions with `chrom`, `start`,
  `end` and optionally `name` and `strand`.

- annotation:

  Optional GTF/GFF3 for transcript models.

- bin:

  Bin width in base pairs.

- footprint:

  Distance each tag is extended toward the 3' end.

- normalize:

  `"tpm"` scales each library to tags per million; `"none"` leaves raw
  counts.

- library_sizes:

  Named vector of total tag counts per sample. Computed from the file
  when absent, which counts only tags inside the window unless the file
  is region-restricted already.

- aggregate:

  How to combine the samples of a group: `"mean"`, `"median"`, `"sum"`,
  or `"none"` to keep one track per sample.

- strand:

  Keep only tags on this strand. `"gene"` uses each region's own strand,
  which is usually what you want.

- activity:

  A data frame of start-site activities to draw as arcs, with columns
  `locus_id` (or `name`), `group`, `x0` (the start site), `x1` (the
  anchor), `count`, and optionally `role`.

- group_col, label_col:

  Manifest columns for grouping and labels.

- gene:

  Restrict the annotation to this gene.

## Value

A `sashimi_data` object.

## Details

For each sample the tags on the requested strand inside the window are
binned at `bin` base pairs and scaled to tags per million, \\t_i \times
10^6 / N_j\\, using that library's total tag count. Each tag is then
extended `footprint` base pairs toward the 3' end, which gives the track
the body of a coverage profile without moving the start site. Finally
the samples of a group are averaged.

Arcs in this mode are not junctions - a 5'-tag protocol produces no
spliced reads to count. An arc here is a pointer from a start site to
the body of the transcript it starts, and its label is that site's
activity. Supply those numbers through `activity`; the resulting figure
is the one the `omakase` sashimi style was designed for.

## Examples

``` r
tags <- data.frame(chrom = "chr1", pos = c(1200, 1200, 1800),
                   count = 1, strand = "+", sample = "s1", group = "early")
sashimi_from_tags(tags, "chr1:1000-2000", bin = 25, footprint = 100)
#> <sashimi_data>: 1 locus, 1 group
#> • tracks: 41 rows
#> loci: chr1:1,000-2,000
```
