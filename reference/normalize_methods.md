# Coverage normalisation methods

Coverage normalisation methods

## Usage

``` r
normalize_methods()
```

## Value

A character vector of method names accepted by
[`normalize_tracks()`](https://EthanShenx.github.io/omakase/reference/normalize_tracks.md)
and the `normalize` argument of
[`plot_sashimi()`](https://EthanShenx.github.io/omakase/reference/plot_sashimi.md).

## Examples

``` r
normalize_methods()
#> [1] "none"        "cpm"         "rpm"         "rpkm"        "size_factor"
#> [6] "manual"      "max"         "sum"        
```
