# Build a coverage track

A filled profile with its value range printed as a bracket on the left,
the way pyGenomeTracks labels a bedGraph.

## Usage

``` r
track_coverage(
  coverage,
  title = NULL,
  color = "#435469",
  group = NULL,
  ymax = NULL,
  show_range = TRUE,
  log_y = FALSE,
  height = 0.7
)
```

## Arguments

- coverage:

  A data frame with `pos` and `value`, a `sashimi_data` object, or a
  path to a bedGraph/bigWig-style file.

- title:

  Track title.

- color:

  Fill colour.

- group:

  When `coverage` is a `sashimi_data`, which group to draw.

- ymax:

  Upper limit; `NULL` uses the track's own maximum.

- show_range:

  Print the `0`-to-`ymax` bracket at the left.

- log_y:

  Draw on a `log10(1 + value)` axis.

- height:

  Relative height of this track.

## Value

An `omakase_track` object.

## Examples

``` r
cv <- data.frame(pos = 1:100, value = abs(sin(seq(0, pi, length.out = 100))))
track_coverage(cv, title = "RPM")
#> <omakase track> coverage - RPM (100 rows)
```
