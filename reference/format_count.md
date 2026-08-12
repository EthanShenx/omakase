# Format an integer count

Junction read counts are whole numbers, so they get thousands separators
and nothing else. Used as the default `arc_label_format` when values
come from a BAM rather than from a normalised activity table.

## Usage

``` r
format_count(v, big_mark = ",")
```

## Arguments

- v:

  Numeric vector.

- big_mark:

  Thousands separator.

## Value

A character vector.

## Examples

``` r
format_count(c(3, 45, 1200))
#> [1] "3"     "45"    "1,200"
```
