# Read a splice junction file

Reads STAR `SJ.out.tab`, a regtools/TopHat junction BED, or a plain
BED-like table of junctions. The format is detected from the columns.

## Usage

``` r
read_junctions(
  path,
  format = c("auto", "star", "bed", "plain"),
  region = NULL,
  min_count = 1,
  include_multimappers = FALSE,
  annotated_only = FALSE
)
```

## Arguments

- path:

  Path to the junction file.

- format:

  `"auto"`, `"star"`, `"bed"`, or `"plain"`.

- region:

  Optional region to restrict to.

- min_count:

  Drop junctions with fewer supporting reads.

- include_multimappers:

  For STAR input, add the multi-mapping read count to the unique count.

- annotated_only:

  For STAR input, keep only junctions the annotation already contains.

## Value

A data frame with `chrom`, `x0`, `x1`, `count`, `strand`.

## Details

STAR's `SJ.out.tab` has no header and nine columns: contig, first intron
base, last intron base, strand code (`0` undefined, `1` `+`, `2` `-`),
intron motif, annotated flag, uniquely-mapping reads, multi-mapping
reads, maximum overhang. The uniquely-mapping count is used unless
`include_multimappers` is set.

A regtools junction BED has twelve columns, where the thick start/end
and block sizes describe the anchors flanking the intron; the intron
itself is recovered from them.

## Examples

``` r
f <- tempfile()
writeLines("chr1\t100\t200\t1\t1\t1\t42\t3\t30", f)
read_junctions(f)
#>   chrom x0  x1 count strand
#> 1  chr1 99 201    42      +
```
