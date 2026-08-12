# Normalise coverage tracks

Rescales the `value` column of a `sashimi_data` object's `tracks` slot
(and, where it makes sense, the matching junction counts) so panels from
libraries of different depth are comparable.

## Usage

``` r
normalize_tracks(
  x,
  method = "none",
  library_sizes = NULL,
  factors = NULL,
  bin_width = NULL,
  scale_junctions = TRUE
)
```

## Arguments

- x:

  A `sashimi_data` object.

- method:

  One of
  [`normalize_methods()`](https://EthanShenx.github.io/omakase/reference/normalize_methods.md).

- library_sizes:

  Named numeric vector of \\N_j\\, named by sample. If absent for
  `cpm`/`rpm`/`rpkm`, the sum of the plotted signal is used instead, and
  a note is emitted - a window-local total is not a library size, and
  only makes panels comparable in a rough sense.

- factors:

  Named numeric vector of divisors for `method = "manual"`.

- bin_width:

  Bin width in base pairs, needed for `rpkm`. Taken from the track
  spacing when `NULL`.

- scale_junctions:

  If `TRUE`, junction counts are divided by the same per-sample factor
  as the coverage, so arc labels stay on the same scale as the track
  beneath them.

## Value

The `sashimi_data` object with rescaled values and a record of the
scaling in `meta`.

## Details

With \\c_i\\ the raw signal in bin \\i\\ of sample \\j\\, \\N_j\\ the
library size and \\L\\ the bin width in kilobases:

- `none`:

  \\c_i\\

- `cpm`, `rpm`:

  \\c_i \times 10^6 / N_j\\

- `rpkm`:

  \\c_i \times 10^9 / (N_j L)\\

- `size_factor`:

  \\c_i / s_j\\, with DESeq2's median-of-ratios factor \\s_j =
  \mathrm{median}\_i \left( k\_{ij} / (\prod_v k\_{iv})^{1/m} \right)\\
  over the \\m\\ samples, taken across bins with non-zero signal in
  every sample

- `manual`:

  \\c_i / f_j\\ for user-supplied factors `factors`

- `max`:

  each sample scaled to a maximum of 1

- `sum`:

  each sample scaled to sum to 1

## Examples

``` r
sd <- sashimi_data(
  loci = data.frame(locus_id = "a", gene_name = "A", chrom = "chr1",
                    strand = "+", win_lo = 1, win_hi = 100),
  tracks = data.frame(locus_id = "a", group = "g", sample = "s1",
                      pos = 1:100, value = 1)
)
normalize_tracks(sd, "cpm")$tracks$value[1]
#> ! No `library_sizes` given; using the total signal in the plotted window.
#> ℹ Pass real library sizes for a normalisation that is comparable between genes.
#> [1] 10000
```
