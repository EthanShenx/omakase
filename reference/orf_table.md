# Call ORFs for a set of transcript models

Convenience wrapper that splices each transcript and calls its ORF,
giving the per-isoform table
[`classify_consequence()`](https://EthanShenx.github.io/omakase/reference/classify_consequence.md)
consumes.

## Usage

``` r
orf_table(models, genome, rule = "first", min_aa = 10, require_stop = TRUE)
```

## Arguments

- models:

  A data frame with `tx_id`, `chrom`, `strand`, `start`, `end`, one row
  per exon.

- genome:

  A `BSgenome`, `DNAStringSet`, or path to an indexed FASTA.

- rule, min_aa, require_stop:

  Passed to
  [`find_orf()`](https://EthanShenx.github.io/omakase/reference/find_orf.md).

## Value

A data frame with one row per transcript: `tx_id`, `chrom`, `strand`,
`tx_start`, `tx_end`, `first_exon_start`, `first_exon_end`, `tss_pos`,
`protein`, `aa_len`, `utr5_len`, `n_uATG`.

## Examples

``` r
if (requireNamespace("Biostrings", quietly = TRUE)) {
  # A toy genome: a start codon, a short coding stretch, then a stop.
  genome <- Biostrings::DNAStringSet(c(
    chr1 = paste0(strrep("T", 60), "ATG", strrep("GCA", 30), "TAA",
                  strrep("T", 60))
  ))
  models <- data.frame(
    tx_id = "tx1", chrom = "chr1", strand = "+",
    start = c(1, 100), end = c(80, 214)
  )
  orf_table(models, genome)
}
#>   tx_id chrom strand tx_start tx_end first_exon_start first_exon_end tss_pos
#> 1   tx1  chr1      +        1    214                1             80       1
#>   protein aa_len utr5_len n_uATG
#> 1              0       NA      0
```
