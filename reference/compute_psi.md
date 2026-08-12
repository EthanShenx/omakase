# Compute PSI (percent spliced in) per locus and group

Compute PSI (percent spliced in) per locus and group

## Usage

``` r
compute_psi(
  x,
  main = "main",
  alt = c("alt", "ATSS", "alternative"),
  method = c("ratio", "rmats"),
  len_inc = NULL,
  len_skip = NULL
)
```

## Arguments

- x:

  A `sashimi_data` object whose `junctions` slot carries a `role`
  column, or a data frame with columns `locus_id`, `group`, `role` and
  `count`.

- main, alt:

  The `role` values treated as the reference and alternative features.

- method:

  `"ratio"` (the default) or `"rmats"`.

- len_inc, len_skip:

  Effective lengths \\\ell_I\\, \\\ell_S\\ used by `method = "rmats"`.

## Value

For a `sashimi_data` input, the object with its `psi` slot filled in.
For a data frame, a data frame of `locus_id`, `group`, `psi`,
`numerator`, `denominator`.

## Details

The `"ratio"` method is the share of activity carried by the reference
feature, \$\$\psi = \frac{A\_{\mathrm{main}}}{A\_{\mathrm{main}} +
A\_{\mathrm{alt}}},\$\$ which is well defined for TSS activities,
junction counts, or anything else additive. Where both terms are zero,
\\\psi\\ is `NA` rather than 0 - no signal is not the same as no
inclusion.

The `"rmats"` method is the length-corrected inclusion level \$\$\psi =
\frac{I / \ell_I}{I / \ell_I + S / \ell_S},\$\$ where \\I\\ and \\S\\
are inclusion and skipping counts and \\\ell_I\\, \\\ell_S\\ the
effective lengths, defaulting to the values rMATS itself uses when they
are not supplied.

## Examples

``` r
j <- data.frame(
  locus_id = "a", group = c("early", "early", "late", "late"),
  role = c("main", "alt", "main", "alt"), count = c(90, 10, 20, 80)
)
compute_psi(j)
#>   locus_id group psi numerator denominator
#> 1        a early 0.9        90         100
#> 2        a  late 0.2        20         100
```
