# Topic proportions over time

Computes \\P(z \mid t) \propto p(t \mid \psi_z) P(z)\\, the expected
topic mixture at each time point, with \\P(z)\\ estimated from the
fitted topic sizes. Times outside the training range are clamped to its
ends.

## Usage

``` r
topic_time(object, times = NULL)
```

## Arguments

- object:

  A fitted `"tot"` model.

- times:

  Time points on the scale of the training `time`. Default: 200 evenly
  spaced points over the training range.

## Value

A matrix with one row per time point and one column per topic; rows sum
to 1. The time points are stored in attribute `"times"`.

## Examples

``` r
set.seed(1)
sim <- simulate_tot(n_docs = 100, n_topics = 3)
fit <- tot(sim$x, sim$time, k = 3, n_iter = 50)
round(topic_time(fit, c(0.1, 0.5, 0.9)), 3)
#>      Topic1 Topic2 Topic3
#> [1,]  0.000  0.001  0.999
#> [2,]  0.019  0.979  0.002
#> [3,]  1.000  0.000  0.000
#> attr(,"times")
#> [1] 0.1 0.5 0.9
```
