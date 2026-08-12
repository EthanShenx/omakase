# Build a coordinate axis track

Build a coordinate axis track

## Usage

``` r
track_axis(
  unit = c("auto", "bp", "Kb", "Mb"),
  n = 6,
  show_chrom = TRUE,
  height = 0.42
)
```

## Arguments

- unit:

  `"auto"`, `"bp"`, `"Kb"` or `"Mb"`.

- n:

  Target number of ticks.

- show_chrom:

  Print the contig name under the axis.

- height:

  Relative height of this track.

## Value

An `omakase_track` object.

## Examples

``` r
track_axis()
#> <omakase track> axis
```
