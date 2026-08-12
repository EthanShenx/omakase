# Call the open reading frame of an mRNA

Finds the ORF and returns the peptide it encodes together with the
length of the 5' untranslated region and the number of upstream AUGs.

## Usage

``` r
find_orf(mrna, rule = c("first", "longest"), min_aa = 10, require_stop = TRUE)
```

## Arguments

- mrna:

  A `DNAString`, `DNAStringSet` of length 1, or a character string.

- rule:

  `"first"` or `"longest"`.

- min_aa:

  Ignore reading frames shorter than this many amino acids.

- require_stop:

  Require an in-frame stop codon. When `FALSE`, a frame running off the
  3' end still counts, which matters for transcripts whose 3' end was
  not sequenced.

## Value

A list with `protein` (character, `""` when no ORF was found), `start`,
`end`, `utr5_len`, `aa_len` and `n_uATG`.

## Details

Two conventions are available. `"first"` takes the first AUG in the
sequence that opens a reading frame terminated by a stop codon, which is
the scanning model and the one used to classify start-site consequences.
`"longest"` takes the AUG that yields the longest peptide, which is more
robust when the 5' end of the transcript is not trustworthy.

## Examples

``` r
find_orf("AAAATGGGGCCCTAA")
#> $protein
#> [1] ""
#> 
#> $start
#> [1] NA
#> 
#> $end
#> [1] NA
#> 
#> $utr5_len
#> [1] NA
#> 
#> $aa_len
#> [1] 0
#> 
#> $n_uATG
#> [1] 0
#> 
```
