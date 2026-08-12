# Summarise a classified switch table

Counts and proportions per category and subtype, in the order the donut
draws them.

## Usage

``` r
consequence_summary(
  x,
  by = c("subtype", "category"),
  filter = NULL,
  drop = c("Unclassified", "Distal (excluded)")
)
```

## Arguments

- x:

  A data frame from
  [`classify_consequence()`](https://EthanShenx.github.io/omakase/reference/classify_consequence.md),
  or any table with `category` and `subtype`.

- by:

  `"subtype"` (the default) or `"category"`.

- filter:

  An optional expression evaluated in `x` to subset rows before
  counting - for example `both_full_length == 1`.

- drop:

  Categories to exclude. Defaults to the two non-biological ones.

## Value

A data frame with `category`, `subtype` (when `by = "subtype"`), `n`,
`proportion` and `label`.

## Examples

``` r
d <- data.frame(category = c("5'UTR change", "5'UTR change", "Promoter swap"),
                subtype = c("longer", "shorter", "alt first exon"))
consequence_summary(d)
#>        category        subtype n                  label proportion
#> 1  5'UTR change         longer 1           Longer 5'UTR   33.33333
#> 2  5'UTR change        shorter 1          Shorter 5'UTR   33.33333
#> 3 Promoter swap alt first exon 1 Alternative first exon   33.33333
```
