# Difference in PSI between two groups

The quantity rMATS reports as `IncLevelDifference`: \\\Delta\psi =
\psi\_{g_1} - \psi\_{g_2}\\, per locus.

## Usage

``` r
delta_psi(x, group1, group2)
```

## Arguments

- x:

  A `sashimi_data` object with a populated `psi` slot, or a data frame
  with `locus_id`, `group`, `psi`.

- group1, group2:

  The two group labels to compare.

## Value

A data frame with `locus_id`, `psi_1`, `psi_2` and `dpsi`.

## Examples

``` r
p <- data.frame(locus_id = c("a", "a"), group = c("early", "late"),
                psi = c(0.9, 0.2))
delta_psi(p, "early", "late")
#>   locus_id psi_1 psi_2 dpsi
#> 1        a   0.9   0.2  0.7
```
