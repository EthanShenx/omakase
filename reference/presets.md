# Plotting presets

Named bundles of arguments for
[`plot_sashimi()`](https://EthanShenx.github.io/omakase/reference/plot_sashimi.md).
A preset sets the house style for a kind of figure; any argument you
pass explicitly overrides it.

## Usage

``` r
presets()
```

## Value

A character vector of preset names.

## Details

- `default`:

  ColorBrewer *Paired*, arcs above the coverage, counts on the curve.

- `tss`:

  Start-site activity figures: *Set2* tracks, a slate/coral pair for the
  two model rows, staggered arc heights at 0.80 and 1.20 of the panel,
  adaptive-precision activity labels and a PSI gutter.

- `junction`:

  Junction counts from alignments: integer labels, width tracking read
  support, span-scaled heights, no PSI gutter.

- `minimal`:

  One ink, unboxed labels, no features.

- `igv`:

  X-spline arcs on a tinted panel, in the manner of IGV.

## Examples

``` r
presets()
#> [1] "default"  "tss"      "junction" "minimal"  "igv"     
```
