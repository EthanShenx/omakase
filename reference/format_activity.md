# Format junction or activity values with adaptive precision

A single number of decimal places never suits a column of junction
counts or TSS activities: the values span several orders of magnitude,
so rounding to whole numbers prints small non-zero values as `"0"`,
while a fixed two decimals makes large counts unreadable.
`format_activity()` picks the precision from the magnitude of each
value, and prints an explicit `"<0.01"` so that a small non-zero value
is never shown as absent.

## Usage

``` r
format_activity(v, big_mark = ",", zero = "0")
```

## Arguments

- v:

  Numeric vector.

- big_mark:

  Thousands separator used for values of 10 or more.

- zero:

  String used for exact zeros.

## Value

A character vector the same length as `v`.

## Details

The rule is:

\$\$f(v) = \begin{cases} \texttt{"0"} & v = 0 \\ \texttt{"\<0.01"} & 0
\< v \< 0.01 \\ \mathrm{fixed}(v, 2) & 0.01 \le v \< 1 \\
\mathrm{fixed}(v, 1) & 1 \le v \< 10 \\
\mathrm{comma}(\mathrm{round}(v)) & v \ge 10 \end{cases}\$\$

## Examples

``` r
format_activity(c(0, 0.004, 0.37, 4.2, 103.01, 4218))
#> [1] "0"     "<0.01" "0.37"  "4.2"   "103"   "4,218"
```
