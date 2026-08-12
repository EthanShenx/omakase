# The omakase theme

A blank canvas for genome-track panels: no axes, no grid, no legend,
tight margins, and every piece of text black, unbolded and at the base
size. Track panels sit directly on top of one another, so any horizontal
rule or panel border between them reads as a feature of the data.

## Usage

``` r
theme_omakase(
  base_size = 9,
  base_family = "",
  margin = c(1, 4, 1, 4),
  legend = "none"
)
```

## Arguments

- base_size:

  Base font size in points. Everything else is expressed relative to it.

- base_family:

  Font family. The empty string uses the device default.

- margin:

  Plot margin in points, given as a length-4 numeric in the order top,
  right, bottom, left.

- legend:

  Legend position; `"none"` by default.

## Value

A
[`ggplot2::theme()`](https://ggplot2.tidyverse.org/reference/theme.html)
object.

## Examples

``` r
library(ggplot2)
ggplot(mtcars, aes(wt, mpg)) + geom_point() + theme_omakase()
```
