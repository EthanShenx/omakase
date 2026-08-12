# Save a genome-track figure

Save a genome-track figure

## Usage

``` r
save_tracks(plot, file, width = 8, height = NULL, dpi = 300, ...)
```

## Arguments

- plot:

  A plot from
  [`plot_tracks()`](https://EthanShenx.github.io/omakase/reference/plot_tracks.md).

- file:

  Output path; the extension selects the device.

- width:

  Width in inches.

- height:

  Height in inches, or `NULL` to size it from the track count.

- dpi:

  Resolution for raster devices.

- ...:

  Passed to
  [`ggplot2::ggsave()`](https://ggplot2.tidyverse.org/reference/ggsave.html).

## Value

`file`, invisibly.

## Examples

``` r
sd <- sashimi_data(
  loci = data.frame(locus_id = "a", gene_name = "A", chrom = "chr1",
                    strand = "+", win_lo = 1, win_hi = 100),
  models = data.frame(locus_id = "a", tx_id = "t", start = 10, end = 90)
)
save_tracks(plot_tracks(sd), file.path(tempdir(), "tracks.pdf"))
```
