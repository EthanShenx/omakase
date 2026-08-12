# Splice exons into an mRNA sequence

Concatenates exon sequences in transcription order and
reverse-complements them for a minus-strand transcript, giving the
mature mRNA an ORF should be called on.

## Usage

``` r
splice_mrna(exons, strand, genome, chrom)
```

## Arguments

- exons:

  A data frame with `start` and `end`, one row per exon.

- strand:

  `"+"` or `"-"`.

- genome:

  A `BSgenome` object, a `DNAStringSet` of chromosomes, or a path to an
  indexed FASTA.

- chrom:

  The contig the exons lie on.

## Value

A `DNAString`.

## Examples

``` r
if (requireNamespace("Biostrings", quietly = TRUE)) {
  gen <- Biostrings::DNAStringSet(c(chr1 = strrep("ATGC", 100)))
  splice_mrna(data.frame(start = c(1, 51), end = c(20, 70)), "+", gen, "chr1")
}
#> 40-letter DNAString object
#> seq: ATGCATGCATGCATGCATGCGCATGCATGCATGCATGCAT
```
