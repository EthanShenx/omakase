# Save a sashimi figure

A thin wrapper over
[`ggplot2::ggsave()`](https://ggplot2.tidyverse.org/reference/ggsave.html)
that defaults the height to the one the panel stack wants - `0.95`
inches per coverage panel plus `2.0` for the annotation - and picks a
device from the file extension.

## Usage

``` r
save_sashimi(plot, file, width = 5.9, height = NULL, dpi = 300, ...)
```

## Arguments

- plot:

  A plot from
  [`plot_sashimi()`](https://EthanShenx.github.io/omakase/reference/plot_sashimi.md).

- file:

  Output path; the extension selects the device.

- width:

  Width in inches.

- height:

  Height in inches, or `NULL` to compute it.

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
                    strand = "+", win_lo = 1, win_hi = 10),
  tracks = data.frame(locus_id = "a", group = "g", pos = 1:10, value = 1:10)
)
f <- file.path(tempdir(), "demo.pdf")
save_sashimi(plot_sashimi(sd), f)
```
