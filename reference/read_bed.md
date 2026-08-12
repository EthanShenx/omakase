# Read a BED file

Reads BED3 through BED12. For BED12 the blocks are expanded into exons
and each exon is split against the thick range, so the result carries a
`feature` column of `"CDS"` and `"UTR"` that a genome track can draw at
different heights.

## Usage

``` r
read_bed(path, region = NULL, expand = TRUE, name_col = 4)

read_bed12(path, region = NULL, expand = TRUE, name_col = 4)
```

## Arguments

- path:

  Path to the BED file, optionally gzipped.

- region:

  Optional region to restrict to.

- expand:

  For BED12, expand blocks into one row per exon (the default). `FALSE`
  keeps one row per transcript.

- name_col:

  Column holding the transcript name; defaults to column 4.

## Value

A data frame with `chrom`, `start`, `end`, `tx_id`, `strand`, `feature`,
and `color` when the file carries an `itemRgb` column. Coordinates are
1-based inclusive.

## Examples

``` r
f <- tempfile()
writeLines(paste("chr1", 1000, 5000, "tx1", 0, "+", 1200, 4000,
                 "31,120,180", 2, "500,800,", "0,3200,", sep = "\t"), f)
read_bed(f)
#>   chrom start  end tx_id strand   color feature
#> 1  chr1  1001 1200   tx1      + #1F78B4     UTR
#> 2  chr1  1201 1500   tx1      + #1F78B4     CDS
#> 3  chr1  4201 5000   tx1      + #1F78B4     UTR
```
