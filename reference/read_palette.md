# Read a palette file

Reads a one-colour-per-line file, the same format `ggsashimi` accepts,
so an existing palette file works unchanged. R colour names and hex
values are both valid, and only the first column is read.

## Usage

``` r
read_palette(path)
```

## Arguments

- path:

  Path to the file, or `NULL` to return `NULL`.

## Value

A character vector of colours, or `NULL`.

## Examples

``` r
p <- tempfile()
writeLines(c("orange", "cornflowerblue", "#008000"), p)
read_palette(p)
#> [1] "orange"         "cornflowerblue" "#008000"       
```
