# Write junctions to a BED file

Writes the `junctions` slot as a BED12 file in the layout regtools and
TopHat use: two anchor blocks either side of the intron, the score
column carrying the supporting read count. This is the counterpart of
`ggsashimi`'s `--junctions-bed`, and lets counted junctions be loaded
into a genome browser or fed back in through
[`read_junctions()`](https://EthanShenx.github.io/omakase/reference/read_junctions.md).

## Usage

``` r
write_junctions(x, file, anchor = 25, name_prefix = "JUNC")
```

## Arguments

- x:

  A `sashimi_data` object, or a junction data frame with `x0`, `x1` and
  `count`.

- file:

  Output path.

- anchor:

  Width of the flanking blocks, in base pairs.

- name_prefix:

  Prefix for the junction names.

## Value

`file`, invisibly.

## Examples

``` r
j <- data.frame(locus_id = "a", group = "g", x0 = 1000, x1 = 2000,
                count = 42, chrom = "chr1")
f <- tempfile(fileext = ".bed")
write_junctions(j, f)
readLines(f)
#> [1] "chr1\t975\t2024\tJUNC00001\t42\t.\t975\t2024\t255,0,0\t2\t25,25,\t0,1024,"
```
