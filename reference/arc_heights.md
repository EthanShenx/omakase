# Compute arc apex heights from junction counts

Compute arc apex heights from junction counts

## Usage

``` r
arc_heights(
  count,
  rule = "constant",
  min_h = 0,
  max_h = 1,
  base = 10,
  span = NULL
)
```

## Arguments

- count:

  Numeric vector of junction counts or activities.

- rule:

  One of
  [`arc_height_rules()`](https://EthanShenx.github.io/omakase/reference/arc_height_rules.md).
  `"constant"` gives every arc the same height, which keeps the picture
  legible when counts span orders of magnitude; `"span"` scales height
  with the width of the junction, so nested junctions nest visually
  instead of crossing; the rest scale height with count.

- min_h, max_h:

  The height range to map onto, in coverage-axis units.

- base:

  Logarithm base used by `rule = "log"`.

- span:

  Junction widths, required for `rule = "span"`.

## Value

A numeric vector of heights, the same length as `count`.

## Examples

``` r
arc_heights(c(1, 10, 100), rule = "log", min_h = 1, max_h = 5)
#> [1] 1.000000 2.738664 5.000000
```
