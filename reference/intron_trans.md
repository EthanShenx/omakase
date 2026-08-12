# A ggplot2 axis transform for compressed genomic coordinates

Wraps an
[`intron_map()`](https://EthanShenx.github.io/omakase/reference/intron_map.md)
as a `scales` transform so that `scale_x_continuous(transform = ...)`
places tick marks at round *genomic* coordinates while drawing the data
in compressed space. This is what keeps the axis honest when introns are
shrunk.

## Usage

``` r
intron_trans(map, n = 5)
```

## Arguments

- map:

  An `omakase_intron_map`.

- n:

  Target number of axis breaks.

## Value

A transform object from
[`scales::new_transform()`](https://scales.r-lib.org/reference/new_transform.html).

## Examples

``` r
m <- intron_map(data.frame(start = 200, end = 8000), 100, 10000)
tr <- intron_trans(m)
tr$breaks(c(100, 10000))
#> [1] 8269.768 9269.768
```
