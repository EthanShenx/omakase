# Coerce an object to sashimi data

Used by the plotting functions so they accept a `sashimi_data` object, a
list of slot tables, or a directory path interchangeably.

## Usage

``` r
as_sashimi_data(x, ...)
```

## Arguments

- x:

  A `sashimi_data` object, a named list of slot data frames, or a path
  to a directory of tables.

- ...:

  Passed to
  [`sashimi_from_tables()`](https://EthanShenx.github.io/omakase/reference/sashimi_from_tables.md).

## Value

A `sashimi_data` object.

## Examples

``` r
as_sashimi_data(list(loci = data.frame(locus_id = "a", gene_name = "A",
  chrom = "chr1", strand = "+", win_lo = 1, win_hi = 10)))
#> <sashimi_data>: 1 locus, 0 groups
#> loci: A
```
