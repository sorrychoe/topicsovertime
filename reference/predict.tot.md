# Predict topics or timestamps for new documents

Runs Gibbs sampling on new documents with the fitted topic-word and
topic-time distributions held fixed ("fold-in").

## Usage

``` r
# S3 method for class 'tot'
predict(
  object,
  newdata,
  type = c("topics", "time"),
  time = NULL,
  n_iter = 100L,
  n_grid = 512L,
  ...
)
```

## Arguments

- object:

  A fitted `"tot"` model.

- newdata:

  A document-term matrix (see
  [`tot()`](https://sorrychoe.github.io/topicsovertime/reference/tot.md)).
  Columns are matched to the training vocabulary by name when both have
  column names; unknown words are dropped. Without names, columns must
  be in training order.

- type:

  `"topics"` for document-topic proportions, or `"time"` for predicted
  timestamps.

- time:

  Optional timestamps of `newdata`, used by `type = "topics"` so that
  time also informs the topic assignments.

- n_iter:

  Number of fold-in sweeps.

- n_grid:

  Number of grid points for `type = "time"`.

- ...:

  Unused.

## Value

For `type = "topics"`, a documents x topics matrix whose rows sum to 1.
For `type = "time"`, a vector of timestamps on the scale of the training
`time` (`NA` for documents without known words).

## Details

With `type = "time"`, topics are sampled from the words alone and the
timestamp is chosen to maximize \\\sum_i \log p(t \mid \psi\_{z_i})\\
over a grid spanning the training time range, as in Section 2 of Wang
and McCallum (2006).

## Examples

``` r
set.seed(1)
sim <- simulate_tot(n_docs = 100, n_topics = 3)
fit <- tot(sim$x, sim$time, k = 3, n_iter = 50)
head(predict(fit, sim$x[1:5, ], type = "topics"))
#>         Topic1    Topic2    Topic3
#> doc1 0.1766667 0.3266667 0.4966667
#> doc2 0.2466667 0.5166667 0.2366667
#> doc3 0.6366667 0.1866667 0.1766667
#> doc4 0.1766667 0.2766667 0.5466667
#> doc5 0.6166667 0.2066667 0.1766667
predict(fit, sim$x[1:5, ], type = "time")
#>      doc1      doc2      doc3      doc4      doc5 
#> 0.2622810 0.4177979 0.8037796 0.1835857 0.8300114 
```
