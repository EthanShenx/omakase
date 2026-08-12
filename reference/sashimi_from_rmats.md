# Build a sashimi data object from rMATS events and alignments

Reads an rMATS event file, turns each event into a plotting window with
the event's exons as transcript models, and reads coverage and junctions
from the supplied BAM files over those windows. The result is one locus
per event, ready for
[`plot_sashimi_all()`](https://EthanShenx.github.io/omakase/reference/plot_sashimi_all.md).

## Usage

``` r
sashimi_from_rmats(
  events,
  bam,
  event_type = NULL,
  fdr = NULL,
  min_dpsi = NULL,
  top = NULL,
  genes = NULL,
  chrom_style = "keep",
  flank = 300,
  group_col = NULL,
  label_col = NULL,
  psi_from = c("rmats", "junctions", "none"),
  ...
)
```

## Arguments

- events:

  Path to an rMATS output file, or a data frame from
  [`read_rmats()`](https://EthanShenx.github.io/omakase/reference/read_rmats.md).

- bam:

  A manifest or vector of BAM paths, as for
  [`sashimi_from_bam()`](https://EthanShenx.github.io/omakase/reference/sashimi_from_bam.md).

- event_type:

  One of
  [`event_types()`](https://EthanShenx.github.io/omakase/reference/event_types.md);
  inferred when `NULL`.

- fdr:

  Keep only events with `FDR` at or below this.

- min_dpsi:

  Keep only events whose absolute `IncLevelDifference` is at least this.

- top:

  Keep only the `top` most significant events.

- genes:

  Restrict to these gene names or IDs.

- chrom_style:

  Rewrite contig names: `"keep"`, `"ucsc"` (ensure a `chr` prefix) or
  `"ensembl"` (strip it), for when the event file and the BAM header
  disagree.

- flank:

  Padding added around the event when computing its window.

- group_col, label_col:

  Manifest columns for grouping and labels.

- psi_from:

  `"rmats"` takes the inclusion levels rMATS reported; `"junctions"`
  recomputes PSI from the junctions read out of the BAMs; `"none"`
  leaves the slot empty.

- ...:

  Passed to
  [`sashimi_from_bam()`](https://EthanShenx.github.io/omakase/reference/sashimi_from_bam.md)
  (`bin`, `strand`, `min_count`, `min_mapq`, and so on).

## Value

A `sashimi_data` object with one locus per event.

## Details

The two isoforms of an event are drawn as two model rows: the inclusion
form (role `"main"`) carries the alternative exon, the skipping form
(role `"alt"`) does not. That is the same reading as
`rmats2sashimiplot`'s two-row layout, but as data rather than as pixels,
so the rows can be recoloured or relabelled afterwards.

## Examples

``` r
ev <- system.file("extdata", "SE.MATS.JC.txt", package = "omakase")
bams <- system.file("extdata", "samples.tsv", package = "omakase")
if (nzchar(ev) && nzchar(bams)) {
  sd <- sashimi_from_rmats(ev, bams, min_count = 5)
  plot_sashimi(sd)
}

```
