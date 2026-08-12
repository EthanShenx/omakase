# Subset a sashimi data object by locus

Subset a sashimi data object by locus

## Usage

``` r
# S3 method for class 'sashimi_data'
x[i, ...]
```

## Arguments

- x:

  A `sashimi_data` object.

- i:

  Locus identifiers, gene names, or a numeric/logical index into `loci`.

- ...:

  Unused.

## Value

A `sashimi_data` object holding only the selected loci.

## Examples

``` r
sd <- sashimi_data(loci = data.frame(
  locus_id = c("a", "b"), gene_name = c("A", "B"), chrom = "chr1",
  strand = "+", win_lo = 1, win_hi = 10
))
sd["a"]
#> <sashimi_data>: 1 locus, 0 groups
#> loci: A
```
