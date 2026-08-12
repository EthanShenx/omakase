# Format a genomic coordinate

Format a genomic coordinate

## Usage

``` r
format_coord(v, big_mark = ",")
```

## Arguments

- v:

  Numeric vector of coordinates.

- big_mark:

  Thousands separator.

## Value

A character vector.

## Examples

``` r
format_coord(27035000)
#> [1] "27,035,000"
```
