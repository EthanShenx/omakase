# Read a GTF or GFF3 annotation

Imports exon and CDS records and returns the tidy transcript models the
annotation panel draws. Format is detected from the file extension, and
gzip-compressed files are handled transparently.

## Usage

``` r
read_annotation(
  path,
  region = NULL,
  feature = c("exon", "CDS"),
  gene = NULL,
  transcript = NULL
)
```

## Arguments

- path:

  Path to a `.gtf`, `.gff`, `.gff3` file, optionally gzipped.

- region:

  Optional region to restrict the import to. Supplying one is much
  faster than reading a whole-genome annotation, and requires a
  tabix-indexed file for the fastest path.

- feature:

  Which record types to keep. `"exon"` is enough to draw a transcript;
  adding `"CDS"` lets the panel draw coding regions taller than
  untranslated ones.

- gene:

  Restrict to one or more gene names.

- transcript:

  Restrict to one or more transcript identifiers.

## Value

A data frame with columns `tx_id`, `gene_id`, `gene_name`, `chrom`,
`start`, `end`, `strand` and `feature`.

## Examples

``` r
gtf <- system.file("extdata", "annotation.gtf", package = "omakase")
if (nzchar(gtf)) head(read_annotation(gtf))
#> Warning: replacing previous import ‘S4Arrays::makeNindexFromArrayViewport’ by ‘DelayedArray::makeNindexFromArrayViewport’ when loading ‘SummarizedExperiment’
#>               tx_id            gene_id gene_name chrom    start      end strand
#> 1 ENST00000376170.4 ENSG00000136754.12      ABI1 chr10 27149676 27149821      -
#> 2 ENST00000376170.4 ENSG00000136754.12      ABI1 chr10 27149676 27149792      -
#> 3 ENST00000376170.4 ENSG00000136754.12      ABI1 chr10 27112067 27112234      -
#> 4 ENST00000376170.4 ENSG00000136754.12      ABI1 chr10 27112067 27112234      -
#> 5 ENST00000376170.4 ENSG00000136754.12      ABI1 chr10 27065994 27066170      -
#> 6 ENST00000376170.4 ENSG00000136754.12      ABI1 chr10 27065994 27066170      -
#>   feature
#> 1    exon
#> 2     CDS
#> 3    exon
#> 4     CDS
#> 5    exon
#> 6     CDS
```
