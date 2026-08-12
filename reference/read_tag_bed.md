# Read a 5'-tag BED file

Reads single-base transcription start site tags. Any BED-like file
works: columns 1-3 are contig, start, end, column 5 is taken as a count
if present, and column 6 as strand.

## Usage

``` r
read_tag_bed(path, region = NULL, strand = NULL, count_col = 5)
```

## Arguments

- path:

  Path to the BED file, optionally gzipped.

- region:

  Optional region to restrict to.

- strand:

  Keep only tags on this strand (`"+"`, `"-"`, or `NULL` for both).

- count_col:

  Column holding a per-tag count. `NULL` treats every row as one tag.

## Value

A data frame with `chrom`, `pos`, `count`, `strand`.

## Examples

``` r
f <- tempfile()
writeLines(c("chr1\t999\t1000\ttag1\t3\t+"), f)
read_tag_bed(f)
#>   chrom  pos count strand
#> 1  chr1 1000     3      +
```
