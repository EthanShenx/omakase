# Apply or invert an intron compression map

`compress_coords()` is \\\phi\\: genome coordinates to plot coordinates.
`expand_coords()` is \\\phi^{-1}\\. Values outside the mapped window are
extrapolated at scale 1, so an arc anchored just beyond the window still
lands somewhere sensible.

## Usage

``` r
compress_coords(x, map)

expand_coords(x, map)
```

## Arguments

- x:

  Numeric vector of coordinates.

- map:

  An `omakase_intron_map` from
  [`intron_map()`](https://EthanShenx.github.io/omakase/reference/intron_map.md).

## Value

A numeric vector the same length as `x`.

## Examples

``` r
m <- intron_map(data.frame(start = 200, end = 800), 100, 1000)
compress_coords(c(100, 500, 1000), m)
#> [1] 100.0000 244.0226 488.0452
expand_coords(compress_coords(500, m), m)
#> [1] 500
```
