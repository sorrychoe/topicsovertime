# topicsovertime (development version)

* `tot()` gains `burn_in` and `thin` to average estimates over several Gibbs
  samples instead of using only the last one (the default is unchanged).
* `tot()` records the joint log-likelihood every `loglik_every` sweeps
  (default 10); check convergence with `plot(fit, type = "loglik")`.
* `tot()` gains `min_shape`, the lower bound on the Beta shape parameters.
  Values below 1 allow U-shaped topic time profiles.
* `terms()` gains `type = "frex"` and `type = "lift"` for ranking words.
* New `top_docs()` lists the documents with the highest share of each topic.
* New `tidy()` and `augment()` methods (generics, broom, tidytext) and an
  `autoplot()` method (ggplot2). Both packages are optional.
* `citation("topicsovertime")` now cites Wang and McCallum (2006) as well as
  the package.

# topicsovertime 0.1.0

* Initial release: `tot()`, `predict()`, `topic_time()`, `terms()`, `plot()`,
  `simulate_tot()`.
