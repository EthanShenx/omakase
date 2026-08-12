# Build a feature track

A row of blocks - peaks, repeats, motifs, start sites - optionally
labelled. With `shape = "marker"` each feature is drawn as a small
downward pointer rather than a block, which is what a single-base
position wants.

## Usage

``` r
track_features(
  features,
  title = NULL,
  color = OM_FEATURE_FILL,
  labels = TRUE,
  shape = c("block", "marker"),
  collapse = TRUE,
  height = 0.4
)
```

## Arguments

- features:

  A data frame with `start`, `end` and optionally `name`, a
  `sashimi_data` object, or a path to a BED file.

- title:

  Track title.

- color:

  Fill; `"bed_rgb"` uses each record's own colour.

- labels:

  Print feature names.

- shape:

  `"block"` or `"marker"`.

- collapse:

  Draw every feature on one row rather than spreading overlapping ones
  across rows.

- height:

  Relative height of this track.

## Value

An `omakase_track` object.

## Examples

``` r
f <- data.frame(start = c(100, 500), end = c(200, 560),
                name = c("peak1", "peak2"))
track_features(f, title = "peaks")
#> <omakase track> features - peaks (2 rows)
```
