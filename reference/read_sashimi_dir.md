# Read a directory of sashimi tables

The counterpart to
[`write_sashimi_data()`](https://EthanShenx.github.io/omakase/reference/write_sashimi_data.md):
reads back the per-slot files it wrote. Files are matched by the slot
name appearing in the file name, so both `sashimi_tracks.tsv` and a
prefixed `mystudy_tracks.parquet` are found.

## Usage

``` r
read_sashimi_dir(dir, prefix = NULL)
```

## Arguments

- dir:

  Directory to read.

- prefix:

  Optional file name prefix to require.

## Value

A `sashimi_data` object.

## Examples

``` r
d <- tempfile(); dir.create(d)
sd <- sashimi_data(loci = data.frame(locus_id = "a", gene_name = "A",
  chrom = "chr1", strand = "+", win_lo = 1, win_hi = 10))
write_sashimi_data(sd, d)
read_sashimi_dir(d)
#> <sashimi_data>: 1 locus, 0 groups
#> loci: A
```
