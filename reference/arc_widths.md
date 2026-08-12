# Compute arc line widths from junction counts

The classic sashimi convention is that a junction supported by more
reads is drawn with a thicker curve. The default here is \$\$w = w_0
\left(\log_b (n + 1)\right)^{\alpha}\$\$ clamped to `range`. Passing
`rule = "miso"` instead reproduces MISO's \\w = \log(b)\\(n+1)^{0.33}
\times 0.1\\, and `rule = "constant"` disables the effect entirely.

## Usage

``` r
arc_widths(
  count,
  rule = c("log", "miso", "linear", "constant"),
  w0 = 0.5,
  base = 10,
  alpha = 1,
  range = c(0.15, 2)
)
```

## Arguments

- count:

  Numeric vector of junction counts.

- rule:

  `"log"` (the default), `"miso"`, `"linear"`, or `"constant"`.

- w0:

  Scale factor \\w_0\\.

- base:

  Logarithm base \\b\\.

- alpha:

  Exponent \\\alpha\\ applied to the log term.

- range:

  Numeric length-2 clamp applied to the result.

## Value

A numeric vector of line widths.

## Examples

``` r
arc_widths(c(1, 5, 50, 500))
#> [1] 0.1505150 0.3890756 0.8537851 1.3499189
arc_widths(c(1, 5, 50, 500), rule = "constant", w0 = 0.5)
#> [1] 0.5 0.5 0.5 0.5
```
