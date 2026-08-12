# Read an rMATS event table

Reads a `*.MATS.JC.txt` or `*.MATS.JCEC.txt` file and returns one row
per event with the coordinates needed to draw it, plus rMATS's own
inclusion levels and significance.

## Usage

``` r
read_rmats(
  path,
  event_type = NULL,
  fdr = NULL,
  min_dpsi = NULL,
  top = NULL,
  genes = NULL,
  chrom_style = "keep",
  flank = 300
)
```

## Arguments

- path:

  Path to the rMATS output file.

- event_type:

  One of
  [`event_types()`](https://EthanShenx.github.io/omakase/reference/event_types.md).
  Inferred from the file name when `NULL`.

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

## Value

A data frame with one row per event: `event_id`, `event_type`,
`gene_id`, `gene_name`, `chrom`, `strand`, the exon coordinates, the
window (`win_lo`, `win_hi`), and rMATS's `psi_1`, `psi_2`, `dpsi`,
`pvalue`, `fdr`.

## Examples

``` r
f <- system.file("extdata", "SE.MATS.JC.txt", package = "omakase")
if (nzchar(f)) read_rmats(f, "SE", fdr = 0.05, top = 20)
#>   event_id event_type           gene_id gene_name chrom strand upstream_start
#> 1     SE.1         SE "ENSG00000136754"      ABI1 chr10      -       27040527
#>   upstream_end alt_start  alt_end downstream_start downstream_end   win_lo
#> 1     27040712  27044584 27044670         27047991       27048164 27040227
#>     win_hi  psi_1  psi_2   dpsi pvalue   fdr inc_counts_1 skip_counts_1
#> 1 27048464 0.4465 0.6645 -0.218 0.0034 0.021      342,440       237,244
#>   inc_counts_2 skip_counts_2 inc_len skip_len
#> 1      412,443        77,147     198       99
```
