# Fit a Topics over Time model

Fits the Topics over Time (TOT) model of Wang and McCallum (2006) by
collapsed Gibbs sampling. Each topic has a multinomial distribution over
words and a Beta distribution over (normalized) time, so topics are
shaped jointly by word co-occurrence and by when documents were written.

## Usage

``` r
tot(
  x,
  time,
  k,
  alpha = 50/k,
  beta = 0.1,
  n_iter = 500L,
  time_weight = 1,
  verbose = FALSE
)
```

## Arguments

- x:

  A document-term matrix of counts with documents in rows: a base
  `matrix`, a `Matrix` sparse matrix, a quanteda `dfm`, or a tm
  `DocumentTermMatrix`.

- time:

  Document timestamps, one per row of `x`: numeric, `Date` or `POSIXct`.

- k:

  Number of topics.

- alpha:

  Symmetric Dirichlet prior on document-topic proportions. Default
  `50 / k`, as in the paper.

- beta:

  Symmetric Dirichlet prior on topic-word distributions.

- n_iter:

  Number of Gibbs sweeps.

- time_weight:

  Exponent on the time likelihood, balancing time against words. A
  number (`1` gives the sampler exactly as in the paper), or
  `"inverse_length"` for `1 / N_d`, the alternative the paper describes,
  which treats each document's timestamp as one observation. The paper
  reports that results are sensitive to this setting.

- verbose:

  Print progress every 50 iterations.

## Value

An object of class `"tot"`, a list with

- phi:

  `k` x vocabulary matrix of topic-word probabilities.

- theta:

  documents x `k` matrix of document-topic proportions.

- psi:

  `k` x 2 matrix of Beta shape parameters on the normalized time scale.

- topic_size:

  Number of tokens assigned to each topic.

plus the settings used and the information needed to map times back to
the original scale.

## Details

Timestamps are mapped linearly onto the open interval (0, 1) using the
range of `time`. The Beta parameters are re-estimated by the method of
moments after every sweep, as in the paper, with both shapes kept at or
above 1 (same mean, higher concentration) so densities stay bounded and
topics cannot lock onto the earliest or latest documents.

## References

Wang, X. and McCallum, A. (2006). Topics over time: a non-Markov
continuous-time model of topical trends. *Proceedings of the 12th ACM
SIGKDD International Conference on Knowledge Discovery and Data Mining*,
424–433.
[doi:10.1145/1150402.1150450](https://doi.org/10.1145/1150402.1150450)

## See also

[`predict.tot()`](https://sorrychoe.github.io/topicovertime/reference/predict.tot.md),
[`topic_time()`](https://sorrychoe.github.io/topicovertime/reference/topic_time.md),
[`terms.tot()`](https://sorrychoe.github.io/topicovertime/reference/terms.tot.md),
[`plot.tot()`](https://sorrychoe.github.io/topicovertime/reference/plot.tot.md)

## Examples

``` r
set.seed(1)
sim <- simulate_tot(n_docs = 100, n_topics = 3)
fit <- tot(sim$x, sim$time, k = 3, n_iter = 50)
fit
#> Topics over Time model
#>   topics: 3, documents: 100, vocabulary: 100
#>   iterations: 50, alpha: 16.66667, beta: 0.1, time_weight: 1
#>   time range: 0.0355636 to 0.9930231
terms(fit, 5)
#>      Topic1 Topic2 Topic3
#> [1,] "w100" "w6"   "w6"  
#> [2,] "w26"  "w26"  "w75" 
#> [3,] "w12"  "w97"  "w30" 
#> [4,] "w87"  "w33"  "w26" 
#> [5,] "w84"  "w51"  "w86" 
```
