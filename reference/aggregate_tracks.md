# Aggregate replicate tracks into one track per group

Replicates are usually plotted one panel each, but with many samples
that is unreadable. Aggregation collapses each group to a single track,
the way `ggsashimi`'s `--aggr` does, with the option to keep junction
counts summed while the coverage is averaged.

## Usage

``` r
aggregate_tracks(
  x,
  fun = c("mean", "median", "sum", "max", "none"),
  junction_fun = NULL
)
```

## Arguments

- x:

  A `sashimi_data` object.

- fun:

  One of `"mean"`, `"median"`, `"sum"`, `"max"`, or `"none"` to leave
  the object untouched.

- junction_fun:

  Aggregation applied to junction counts; defaults to the same as `fun`.
  `"sum"` is the honest choice when arcs should reflect total support
  across replicates.

## Value

The aggregated `sashimi_data` object.

## Examples

``` r
sd <- sashimi_data(
  loci = data.frame(locus_id = "a", gene_name = "A", chrom = "chr1",
                    strand = "+", win_lo = 1, win_hi = 3),
  tracks = data.frame(locus_id = "a", group = "g",
                      sample = rep(c("s1", "s2"), each = 3),
                      pos = rep(1:3, 2), value = c(1, 2, 3, 3, 4, 5))
)
aggregate_tracks(sd, "mean")$tracks
#>   locus_id group pos strand value sample bin
#> 1        a     g   1      *     2   <NA>  NA
#> 2        a     g   2      *     3   <NA>  NA
#> 3        a     g   3      *     4   <NA>  NA
```
