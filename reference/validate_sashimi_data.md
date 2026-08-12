# Validate a sashimi data object

Checks slot classes and cross-slot referential integrity: every
`locus_id` appearing in `tracks`, `junctions`, `models`, `psi` or
`features` must exist in `loci`.

## Usage

``` r
validate_sashimi_data(x, strict = TRUE)
```

## Arguments

- x:

  A `sashimi_data` object.

- strict:

  If `TRUE`, an unknown `locus_id` is an error; otherwise a warning and
  the offending rows are reported but kept.

## Value

`x`, invisibly.

## Examples

``` r
validate_sashimi_data(sashimi_data())
```
