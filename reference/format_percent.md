# Format a percentage

Whole numbers when the value is exact, one decimal otherwise, so a donut
of round percentages is not littered with trailing zeros.

## Usage

``` r
format_percent(v)
```

## Arguments

- v:

  Numeric vector of percentages.

## Value

A character vector.

## Examples

``` r
format_percent(c(25, 12.5, 33.333))
#> [1] "25%"   "12.5%" "33.3%"
```
