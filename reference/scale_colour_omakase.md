# Discrete colour and fill scales using the omakase palettes

Discrete colour and fill scales using the omakase palettes

## Usage

``` r
scale_colour_omakase(palette = "omakase", ...)

scale_color_omakase(palette = "omakase", ...)

scale_fill_omakase(palette = "omakase", ...)
```

## Arguments

- palette:

  Palette name, see
  [`omakase_palettes()`](https://EthanShenx.github.io/omakase/reference/omakase_palettes.md).

- ...:

  Passed to
  [`ggplot2::discrete_scale()`](https://ggplot2.tidyverse.org/reference/discrete_scale.html).

## Value

A ggplot2 scale.

## Examples

``` r
library(ggplot2)
ggplot(mtcars, aes(wt, mpg, colour = factor(cyl))) +
  geom_point() +
  scale_colour_omakase()
```
