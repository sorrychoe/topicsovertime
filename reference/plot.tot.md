# Plot topic trends or topic time densities

Plot topic trends or topic time densities

## Usage

``` r
# S3 method for class 'tot'
plot(
  x,
  type = c("trend", "density"),
  topics = NULL,
  n = 200L,
  legend = TRUE,
  ...
)
```

## Arguments

- x:

  A fitted `"tot"` model.

- type:

  `"trend"` plots
  [`topic_time()`](https://sorrychoe.github.io/topicovertime/reference/topic_time.md)
  (topic proportions over time); `"density"` plots each topic's Beta
  density over time.

- topics:

  Topics to plot (indices or names). Default: all.

- n:

  Number of time points.

- legend:

  Draw a legend.

- ...:

  Passed to
  [`graphics::matplot()`](https://rdrr.io/r/graphics/matplot.html).

## Value

The plotted matrix, invisibly.

## Examples

``` r
set.seed(1)
sim <- simulate_tot(n_docs = 100, n_topics = 3)
fit <- tot(sim$x, sim$time, k = 3, n_iter = 50)
plot(fit)

plot(fit, type = "density")
```
