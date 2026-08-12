# Build a sashimi data object from junction files

For projects that keep junction tables rather than alignments. Produces
arcs but no coverage, which draws a perfectly readable figure - the arcs
carry the splicing signal, and the filled area is context.

## Usage

``` r
sashimi_from_junctions(
  junctions,
  region,
  annotation = NULL,
  group_col = NULL,
  label_col = NULL,
  min_count = 1,
  format = "auto",
  gene = NULL,
  ...
)
```

## Arguments

- junctions:

  A path, a vector of paths, a manifest read by
  [`read_manifest()`](https://EthanShenx.github.io/omakase/reference/read_manifest.md),
  or a data frame of junctions.

- region:

  A region string, or a data frame of regions.

- annotation:

  Optional GTF/GFF3 for transcript models.

- group_col, label_col:

  Manifest columns for grouping and labels.

- min_count:

  Drop junctions with fewer supporting reads.

- format:

  Junction file format, see
  [`read_junctions()`](https://EthanShenx.github.io/omakase/reference/read_junctions.md).

- gene:

  Restrict the annotation to this gene.

- ...:

  Passed to
  [`read_junctions()`](https://EthanShenx.github.io/omakase/reference/read_junctions.md).

## Value

A `sashimi_data` object with an empty `tracks` slot.

## Examples

``` r
f <- tempfile()
writeLines(c("chr1\t1200\t1800\t1\t1\t1\t42\t0\t30"), f)
sashimi_from_junctions(f, "chr1:1000-2000")
#> <sashimi_data>: 1 locus, 0 groups
#> • junctions: 1 row
#> loci: chr1:1,000-2,000
```
