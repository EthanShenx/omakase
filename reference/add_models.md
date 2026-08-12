# Add transcript models to a sashimi data object

Attaches models from an annotation to loci that already carry coverage,
which is useful when tracks came from a source that knows nothing about
transcripts - a tag BED, or a table of precomputed counts.

## Usage

``` r
add_models(x, annotation, gene = NULL, max_tx = NULL)
```

## Arguments

- x:

  A `sashimi_data` object.

- annotation:

  A GTF/GFF3 path or a data frame from
  [`read_annotation()`](https://EthanShenx.github.io/omakase/reference/read_annotation.md).

- gene:

  Restrict to these gene names.

- max_tx:

  Keep at most this many transcripts per locus, longest first.
  Annotations with thirty isoforms make an unreadable panel.

## Value

The object with its `models` slot populated.

## Examples

``` r
sd <- sashimi_data(loci = data.frame(locus_id = "a", gene_name = "A",
  chrom = "chr1", strand = "+", win_lo = 1, win_hi = 100))
ann <- data.frame(tx_id = "t1", gene_id = "g1", gene_name = "A",
  chrom = "chr1", start = 10, end = 50, strand = "+", feature = "exon")
add_models(sd, ann)$models
#>   tx_id role feature start end strand locus_id
#> 1    t1 <NA>    exon    10  50      +        a
```
