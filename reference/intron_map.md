# Build a genome-to-plot coordinate map that compresses introns

Given a window and a set of intronic intervals, constructs a monotone
piecewise-linear map \\\phi\\ that leaves exonic sequence at scale 1 and
draws each intron of length \\L = b - a\\ at a reduced length \\g(L)\\:

## Usage

``` r
intron_map(
  introns,
  win_lo,
  win_hi,
  method = "power",
  gamma = 0.7,
  cap = 200,
  scale = 5,
  coef = 20,
  min_intron = 100,
  exon_scale = 1
)
```

## Arguments

- introns:

  A two-column data frame or matrix of intron `start`/`end` coordinates,
  or `NULL` for no compression. Overlapping intervals are merged;
  intervals are sorted internally.

- win_lo, win_hi:

  The plotted window. Introns are clipped to it.

- method:

  One of
  [`shrink_methods()`](https://EthanShenx.github.io/omakase/reference/shrink_methods.md).

- gamma:

  Exponent for `method = "power"`.

- cap:

  Maximum drawn intron length for `method = "fixed"`.

- scale:

  Divisor for `method = "scale"`.

- coef:

  Multiplier \\c\\ for `method = "log"`.

- min_intron:

  Introns shorter than this are left uncompressed; shrinking a 60 bp
  intron buys no width and distorts short-intron genes.

- exon_scale:

  Divisor applied to exonic stretches, the counterpart of
  `rmats2sashimiplot`'s `--exon_s`. The default of `1` leaves exons at
  true scale; larger values shrink them too, which is occasionally
  wanted when a single huge exon dominates the view.

## Value

An object of class `omakase_intron_map` with elements `genome` and
`plot` (the knot coordinates), plus the settings used. Apply it with
[`compress_coords()`](https://EthanShenx.github.io/omakase/reference/compress_coords.md)
and invert it with
[`expand_coords()`](https://EthanShenx.github.io/omakase/reference/compress_coords.md).

## Details

- `none`:

  \\g(L) = L\\

- `power`:

  \\g(L) = L^{\gamma}\\, default \\\gamma = 0.7\\ (the `ggsashimi` rule)

- `log`:

  \\g(L) = c \log(1 + L)\\

- `fixed`:

  \\g(L) = \min(L, k)\\

- `scale`:

  \\g(L) = L / s\\\\ (the `rmats2sashimiplot --intron_s` rule)

Writing the cumulative shift after the \\j\\th intron as \\\Delta_j =
\sum\_{i \le j} (L_i - g(L_i))\\, the map is \\\phi(x) = x -
\Delta\_{j(x)}\\ on exonic stretches, and interpolates linearly across
each compressed intron. Because \\\phi\\ is strictly increasing it has
an exact inverse, which is what
[`plot_sashimi()`](https://EthanShenx.github.io/omakase/reference/plot_sashimi.md)
uses to label the axis with true genomic coordinates.

## Examples

``` r
m <- intron_map(data.frame(start = c(120, 400), end = c(300, 900)),
                win_lo = 100, win_hi = 1000)
compress_coords(c(100, 300, 1000), m)
#> [1] 100.0000 157.9045 435.4004
# phi is exactly invertible
expand_coords(compress_coords(c(150, 950), m), m)
#> [1] 150 950
```
