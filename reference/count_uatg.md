# Count upstream AUGs in a 5' UTR

The number of AUG triplets before the main start codon. An AUG gained in
a lengthened 5' UTR can create an upstream open reading frame and
suppress translation of the main one, so this is the quantity that makes
a "longer 5' UTR" consequence mean something.

## Usage

``` r
count_uatg(mrna, orf_start)
```

## Arguments

- mrna:

  A `DNAString` or character string.

- orf_start:

  1-based position of the main start codon.

## Value

An integer count.

## Examples

``` r
count_uatg("ATGCCCAAAATGGGG", orf_start = 10)
#> [1] 1
```
