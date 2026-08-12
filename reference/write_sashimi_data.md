# Write a sashimi data object to disk

Writes one file per non-empty slot, plus a `*_methods.tsv` recording
`meta`. The layout mirrors what
[`sashimi_from_tables()`](https://EthanShenx.github.io/omakase/reference/sashimi_from_tables.md)
reads back, so a slow counting step can be run once and re-plotted
cheaply.

## Usage

``` r
write_sashimi_data(x, dir, prefix = "sashimi", parquet = FALSE)
```

## Arguments

- x:

  A `sashimi_data` object.

- dir:

  Output directory, created if absent.

- prefix:

  File name prefix.

- parquet:

  If `TRUE` and the arrow package is available, write the (largest)
  `tracks` slot as parquet instead of TSV.

## Value

The character vector of written paths, invisibly.

## Examples

``` r
sd <- sashimi_data(loci = data.frame(locus_id = "a", gene_name = "A",
  chrom = "chr1", strand = "+", win_lo = 1, win_hi = 10))
write_sashimi_data(sd, tempdir(), prefix = "demo")
```
