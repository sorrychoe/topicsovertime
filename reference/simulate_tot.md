# Simulate documents from a Topics over Time model

Generates a corpus with one timestamp per document, drawn from the
mixture of topic Beta distributions weighted by the document's topic
proportions (the generative view of Figure 1(b) in Wang and McCallum,
2006). Useful for examples, tests and checking parameter recovery.

## Usage

``` r
simulate_tot(
  n_docs = 200L,
  n_topics = 3L,
  n_words = 100L,
  doc_length = 50L,
  psi = NULL,
  alpha = 0.1,
  beta = 0.1
)
```

## Arguments

- n_docs:

  Number of documents.

- n_topics:

  Number of topics.

- n_words:

  Vocabulary size.

- doc_length:

  Tokens per document.

- psi:

  `n_topics` x 2 matrix of Beta shape parameters. Default: topics
  centered at evenly spaced times with concentration 20.

- alpha:

  Dirichlet prior on document-topic proportions.

- beta:

  Dirichlet prior on topic-word distributions.

## Value

A list with `x` (a sparse documents x words count matrix), `time`
(timestamps in (0, 1)), and the true `phi`, `theta` and `psi`.

## Examples

``` r
set.seed(1)
sim <- simulate_tot(n_docs = 50, n_topics = 2)
dim(sim$x)
#> [1]  50 100
summary(sim$time)
#>    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#> 0.01992 0.24836 0.77532 0.58809 0.86850 0.96212 
```
