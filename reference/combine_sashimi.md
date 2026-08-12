# Combine sashimi data objects

Row-binds matching slots, so several regions read separately can be
plotted or written as one.

## Usage

``` r
combine_sashimi(...)
```

## Arguments

- ...:

  `sashimi_data` objects.

## Value

A single `sashimi_data` object.

## Examples

``` r
a <- sashimi_data(loci = data.frame(locus_id = "a", gene_name = "A",
  chrom = "chr1", strand = "+", win_lo = 1, win_hi = 10))
b <- sashimi_data(loci = data.frame(locus_id = "b", gene_name = "B",
  chrom = "chr2", strand = "-", win_lo = 5, win_hi = 20))
c_ab <- combine_sashimi(a, b)
loci(c_ab)
#> [1] "a" "b"
```
