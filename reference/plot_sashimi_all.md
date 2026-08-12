# Draw a sashimi figure for every locus

Loops
[`plot_sashimi()`](https://EthanShenx.github.io/omakase/reference/plot_sashimi.md)
over the loci in an object, optionally writing one file per locus.

## Usage

``` r
plot_sashimi_all(
  x,
  dir = NULL,
  device = "pdf",
  width = 5.9,
  height = NULL,
  dpi = 300,
  quiet = FALSE,
  ...
)
```

## Arguments

- x:

  A `sashimi_data` object.

- dir:

  Directory to write into, or `NULL` to return the plots without
  writing.

- device:

  File extension: `"pdf"`, `"png"`, `"svg"`, `"tiff"`, `"jpeg"`.

- width:

  Figure width in inches.

- height:

  Figure height in inches, or `NULL` to size it from the number of
  panels.

- dpi:

  Resolution for raster devices.

- quiet:

  Suppress the per-file progress message.

- ...:

  Passed to
  [`plot_sashimi()`](https://EthanShenx.github.io/omakase/reference/plot_sashimi.md).

## Value

A named list of plot objects, invisibly when writing to disk.

## Examples

``` r
sd <- sashimi_data(
  loci = data.frame(locus_id = c("a", "b"), gene_name = c("A", "B"),
                    chrom = "chr1", strand = "+",
                    win_lo = 1, win_hi = 100),
  tracks = data.frame(locus_id = rep(c("a", "b"), each = 10), group = "g",
                      pos = rep(seq(1, 100, length.out = 10), 2),
                      value = runif(20))
)
plots <- plot_sashimi_all(sd)
names(plots)
#> [1] "a" "b"
```
