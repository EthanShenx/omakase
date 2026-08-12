# A themed variant that keeps the axes

Same type and weight as
[`theme_omakase()`](https://EthanShenx.github.io/omakase/reference/theme_omakase.md),
but with a visible x axis and y axis. Use it when a panel is being read
quantitatively - a coverage track whose depth matters, or a standalone
plot outside a sashimi stack.

## Usage

``` r
theme_omakase_axes(
  base_size = 9,
  base_family = "",
  hairline = 0.3,
  margin = c(1, 4, 1, 4),
  legend = "none",
  grid = FALSE
)
```

## Arguments

- base_size:

  Base font size in points. Everything else is expressed relative to it.

- base_family:

  Font family. The empty string uses the device default.

- hairline:

  Line width for the axis rules.

- margin:

  Plot margin in points, given as a length-4 numeric in the order top,
  right, bottom, left.

- legend:

  Legend position; `"none"` by default.

- grid:

  Whether to draw a light horizontal grid.

## Value

A
[`ggplot2::theme()`](https://ggplot2.tidyverse.org/reference/theme.html)
object.

## Examples

``` r
library(ggplot2)
ggplot(mtcars, aes(wt, mpg)) + geom_point() + theme_omakase_axes()
```
