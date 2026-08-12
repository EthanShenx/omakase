# Colour palettes

Colour palettes

## Usage

``` r
omakase_palette(name = "omakase", n = NULL)
```

## Arguments

- name:

  Palette name: `"omakase"` (the default, ColorBrewer *Paired*),
  `"paired_dark"` (its saturated members only), `"set2"` (ColorBrewer
  *Set2*), `"ocean"`, `"ember"`, `"okabe"` (colour-vision-safe), or
  `"mono"`.

- n:

  Number of colours needed. When `n` exceeds the palette length the
  palette is interpolated rather than recycled, so groups stay
  distinguishable.

## Value

A character vector of hex colours.

## Examples

``` r
omakase_palette("omakase", 3)
#> [1] "#A6CEE3" "#1F78B4" "#B2DF8A"
omakase_palette("okabe", 12)
#>  [1] "#E69F00" "#8AAC94" "#3EAEC8" "#07A07D" "#82C458" "#C4CF56" "#2B869D"
#>  [8] "#606861" "#D4600F" "#CE7179" "#814C6A" "#000000"
```
