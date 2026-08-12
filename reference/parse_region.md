# Parse a genomic region string

Accepts the usual browser syntax and the colon-separated form used by
`rmats2sashimiplot`, and is forgiving about thousands separators and
whitespace.

## Usage

``` r
parse_region(x, one_based = TRUE)
```

## Arguments

- x:

  A region string, or a `list`/`data.frame` already carrying `chrom`,
  `start`, `end` (returned unchanged after validation).

- one_based:

  Logical. If `TRUE` (the default) the returned `start` is treated as
  1-based inclusive, matching the way coordinates are written in genome
  browsers and GTF files. Set to `FALSE` to interpret the input as
  0-based half-open BED coordinates, which are then converted
  internally.

## Value

A list with elements `chrom`, `start`, `end`, `strand` and `label`.
`strand` is `"*"` unless the input specified one.

## Details

Recognised forms:

- `"chr10:27035000-27050000"`

- `"chr10:27,035,000-27,050,000"`

- `"chr10:27035000..27050000"`

- `"chr10:+:27035000:27050000"` (rmats2sashimiplot `-c` style, with
  strand)

- `"chr10"` (whole contig; `start`/`end` become `NA`)

## Examples

``` r
parse_region("chr10:27,035,000-27,050,000")
#> <omakase region> chr10:27,035,000-27,050,000
parse_region("chr10:+:27035000:27050000")
#> <omakase region> chr10:27,035,000-27,050,000 (+)
```
