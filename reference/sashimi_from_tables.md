# Build a sashimi data object from tidy tables

Accepts data frames or paths to TSV/CSV/parquet files, one per slot, and
maps their columns onto the
[`sashimi_data()`](https://EthanShenx.github.io/omakase/reference/sashimi_data.md)
contract. Column names that differ from the contract can be remapped
through `rename`, so tables written by an existing pipeline usually need
no editing.

## Usage

``` r
sashimi_from_tables(
  loci = NULL,
  tracks = NULL,
  junctions = NULL,
  models = NULL,
  psi = NULL,
  features = NULL,
  rename = list(),
  locus_col = NULL,
  meta = list()
)
```

## Arguments

- loci, tracks, junctions, models, psi, features:

  Data frames, or paths to files. Any may be `NULL`.

- rename:

  A named list of named character vectors, one per slot, mapping
  contract column names to the names used in your tables. For example
  `rename = list(tracks = c(group = "stage", value = "rpm"))` reads a
  table whose columns are `stage` and `rpm`.

- locus_col:

  The column identifying the locus, if it is not `locus_id`. A common
  case is a gene name column shared by every table.

- meta:

  A named list recorded on the object.

## Value

A `sashimi_data` object.

## Examples

``` r
loci <- data.frame(gene_name = "DEMO", chrom = "chr1", strand = "+",
                   win_lo = 1000, win_hi = 2000)
tracks <- data.frame(gene_name = "DEMO", stage = "early",
                     pos = seq(1000, 2000, 100), rpm = runif(11))
sashimi_from_tables(
  loci = loci, tracks = tracks,
  rename = list(tracks = c(group = "stage", value = "rpm")),
  locus_col = "gene_name"
)
#> <sashimi_data>: 1 locus, 1 group
#> • tracks: 11 rows
#> loci: DEMO
```
