# Build a transcript model track

Draws transcript models the way a genome browser does: exons as blocks,
introns as a line, and the direction of transcription marked by chevrons
along that line. With `style = "UCSC"` the coding part of each exon is
drawn taller than the untranslated part.

## Usage

``` r
track_models(
  models,
  title = NULL,
  color = "#1F78B4",
  style = c("UCSC", "flat"),
  labels = FALSE,
  chevrons = 18,
  height = 1,
  row_by = c("tx_id", "role")
)
```

## Arguments

- models:

  A data frame of exons with `start`, `end` and `tx_id`, such as
  [`read_bed12()`](https://EthanShenx.github.io/omakase/reference/read_bed.md)
  or
  [`read_annotation()`](https://EthanShenx.github.io/omakase/reference/read_annotation.md)
  returns; a `sashimi_data` object, whose `models` slot is used; or a
  path to a BED/GTF/GFF3 file.

- title:

  Track title, printed to the right of the panel.

- color:

  Fill for the blocks. `"bed_rgb"` uses each record's own colour when
  the source carried one; a named vector keyed by `role` colours by
  role.

- style:

  `"UCSC"` draws CDS taller than UTR; `"flat"` draws every block the
  same height.

- labels:

  Print each transcript's name beside its row.

- chevrons:

  Number of direction marks per intron; `0` for none.

- height:

  Relative height of this track.

- row_by:

  Which column separates rows: `"tx_id"` (one row per transcript, the
  default) or `"role"` (one row per role).

## Value

A `omakase_track` object, for
[`plot_tracks()`](https://EthanShenx.github.io/omakase/reference/plot_tracks.md).

## Examples

``` r
m <- data.frame(tx_id = "tx1", start = c(1000, 3000), end = c(1500, 4000),
                strand = "+", feature = "CDS")
track_models(m, title = "models")
#> <omakase track> models - models (2 rows)
```
