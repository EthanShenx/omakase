# Construct a sashimi data object

The container every `sashimi_from_*()` reader returns and every plotting
function accepts. It is a plain list of tidy data frames, so it can be
built by hand when counts come from a source omakase does not read
directly.

## Usage

``` r
sashimi_data(
  loci = NULL,
  tracks = NULL,
  junctions = NULL,
  models = NULL,
  psi = NULL,
  features = NULL,
  meta = list()
)
```

## Arguments

- loci, tracks, junctions, models, psi, features:

  Data frames as described above. Any may be `NULL`.

- meta:

  A named list of provenance/parameter values.

## Value

An object of class `sashimi_data`.

## Slots

- `loci`:

  One row per plotted locus. Required: `locus_id`, `gene_name`, `chrom`,
  `strand`, `win_lo`, `win_hi`. Optional: `anchor` (a shared downstream
  point arcs may point at), `main_apex`, `alt_apex` (start sites to mark
  with a triangle).

- `tracks`:

  Binned coverage. Required: `locus_id`, `group`, `pos`, `value`.
  Optional: `sample`, `strand`, `bin`.

- `junctions`:

  Arcs. Required: `locus_id`, `group`, `x0`, `x1`, `count`. Optional:
  `role`, `label`, `sample`, `strand`. `x0`/`x1` are genomic
  coordinates; `count` drives the arc label, and optionally its width
  and height.

- `models`:

  Transcript models drawn beneath the tracks. Required: `locus_id`,
  `tx_id`, `start`, `end`. Optional: `role` (used to colour and to order
  the rows), `feature` (`"exon"`, `"CDS"`, `"UTR"`), `strand`.

- `psi`:

  Per-group inclusion values printed in the right-hand gutter. Required:
  `locus_id`, `group`, `psi`. Optional: `numerator`, `denominator`.

- `features`:

  An extra annotation row under each model - repeats, motifs, peaks,
  anything. Required: `locus_id`, `start`, `end`, `name`. Optional:
  `role`, `class`, `color`.

- `meta`:

  A named list recording how the object was built. Printed by
  [`summary()`](https://rdrr.io/r/base/summary.html) and written to the
  methods table by
  [`write_sashimi_data()`](https://EthanShenx.github.io/omakase/reference/write_sashimi_data.md).

## Examples

``` r
loci <- data.frame(
  locus_id = "demo", gene_name = "DEMO", chrom = "chr1",
  strand = "+", win_lo = 1000, win_hi = 2000
)
tracks <- data.frame(
  locus_id = "demo", group = "A",
  pos = seq(1000, 2000, by = 10),
  value = abs(sin(seq(0, pi, length.out = 101))) * 10
)
sashimi_data(loci = loci, tracks = tracks)
#> <sashimi_data>: 1 locus, 1 group
#> • tracks: 101 rows
#> loci: DEMO
```
