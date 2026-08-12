# Build the path of a junction arc

Returns the polyline for a single arc running from `x0` to `x1` and
peaking at height `h`. All shapes start and end at `y = 0` and reach
exactly `h` at their apex, so arcs of different shapes are directly
comparable.

## Usage

``` r
arc_path(x0, x1, h, shape = "sine", n = 120, y0 = 0)
```

## Arguments

- x0, x1:

  Genomic start and end of the junction.

- h:

  Apex height, in the units of the coverage axis.

- shape:

  One of
  [`arc_shapes()`](https://EthanShenx.github.io/omakase/reference/arc_shapes.md).

- n:

  Number of points in the returned path. Higher is smoother; 120 is
  plenty at print size, and `elbow` ignores it.

- y0:

  Baseline the arc springs from. Non-zero when arcs are stacked above a
  coverage track that does not start at zero.

## Value

A data frame with columns `x` and `y`.

## Details

The shapes, parameterised by \\t \in \[0, 1\]\\ with \\x(t) = x_0 +
(x_1 - x_0)t\\:

- `sine`:

  \\y = h \sin(\pi t)\\. The default: symmetric, with a gentle take-off
  from the baseline.

- `parabola`:

  \\y = 4 h\\ t (1 - t)\\. Slightly fuller shoulders than the sine.

- `bezier`:

  Cubic Bezier \\B(t) = (1-t)^3 P_0 + 3t(1-t)^2 P_1 + 3t^2(1-t) P_2 +
  t^3 P_3\\ with control points lifted to \\4h/3\\ so the apex lands on
  \\h\\; this is the curve MISO's `sashimi_plot` draws.

- `xspline`:

  An approximation of the X-spline used by `ggsashimi`, which rises
  steeply from the donor and flattens across the middle.

- `elbow`:

  Straight risers joined by a flat top - useful when many arcs overlap
  and curvature becomes noise.

- `arch`:

  A semi-ellipse, \\y = h\sqrt{1 - (2t - 1)^2}\\, the roundest of the
  set.

## Examples

``` r
head(arc_path(100, 200, h = 10))
#>          x         y
#> 1 100.0000 0.0000000
#> 2 100.8333 0.2617695
#> 3 101.6667 0.5233596
#> 4 102.5000 0.7845910
#> 5 103.3333 1.0452846
#> 6 104.1667 1.3052619
# every shape peaks at h
vapply(arc_shapes(), function(s) max(arc_path(0, 1, 5, shape = s)$y), numeric(1))
#>     sine parabola   bezier  xspline    elbow     arch 
#>        5        5        5        5        5        5 
```
