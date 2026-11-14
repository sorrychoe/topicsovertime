# topicsovertime

An R implementation of the **Topics over Time (TOT)** model of
[Wang and McCallum (2006)](https://doi.org/10.1145/1150402.1150450).

TOT is an LDA-style topic model in which every topic has, besides its
distribution over words, a **Beta distribution over continuous time**. Topics
are found jointly from word co-occurrence and from *when* documents were
written:

- Time is neither discretized nor modelled with Markov dynamics, so there is
  no time-slice size to choose.
- The meaning (word distribution) of a topic stays constant. What changes over
  time is how often each topic occurs.
- A fitted model can predict when an undated document was written.

Inference uses collapsed Gibbs sampling written in C++ (via Rcpp).

Documentation: <https://sorrychoe.github.io/topicsovertime/>

## Installation

```r
# install.packages("remotes")
remotes::install_github("sorrychoe/topicsovertime")
```

## Usage

`tot()` takes a document-term matrix of counts (base `matrix`, a `Matrix`
sparse matrix, a quanteda `dfm`, or a tm `DocumentTermMatrix`) and one
timestamp per document (numeric, `Date` or `POSIXct`).

```r
library(topicsovertime)
library(quanteda)

toks <- tokens(data_corpus_inaugural, remove_punct = TRUE, remove_numbers = TRUE)
toks <- tokens_remove(tokens_tolower(toks), stopwords("en"))
x <- dfm_trim(dfm(toks), min_termfreq = 5)
year <- docvars(data_corpus_inaugural, "Year")

set.seed(2026)
fit <- tot(x, year, k = 8, n_iter = 300)

fit                      # summary of the model
terms(fit, 10)           # top words per topic
fit$psi                  # Beta parameters of each topic (normalized time)
plot(fit)                # topic proportions over time, P(topic | time)
plot(fit, type = "density")  # Beta density of each topic over time

topic_time(fit, c(1800, 1900, 2000))      # P(topic | time) at given times
predict(fit, x, type = "topics")          # topic proportions of documents
predict(fit, x, type = "time")            # predicted year of each document
```

To try the model without any text data, `simulate_tot()` generates a
synthetic corpus from known topics and time distributions.

## Main functions

| Function | Description |
|---|---|
| `tot()` | Fit the model by collapsed Gibbs sampling |
| `terms()` | Most probable words per topic |
| `topic_time()` | Topic proportions at given times, `P(z \| t) ∝ p(t \| ψ_z) P(z)` |
| `predict(type = "topics")` | Topic proportions of new documents (fold-in) |
| `predict(type = "time")` | Predicted timestamps of new documents |
| `plot()` | Topic trends or topic time densities |
| `simulate_tot()` | Simulate a corpus from the TOT generative process |

## Choosing `time_weight`

Each word token carries its document's timestamp, so a long document's time
counts many times in the likelihood. `time_weight` rescales the time
likelihood against the words. The paper notes that this setting is sensitive
and has to be tuned.

| `time_weight` | Effect |
|---|---|
| `1` (default) | The Gibbs sampler exactly as in the paper. With long documents, time can dominate: topics become narrow time periods with mixed words. |
| `"inverse_length"` | `1 / N_d`, the alternative suggested in the paper: each document's timestamp counts once. Topics are driven mostly by words. |
| a value in between | A compromise. On ~8,000 Korean news articles (about 200 tokens each), `0.3` gave coherent topics while still improving date prediction. |

Inspect `terms()` and `plot()` under a few values, and, if timestamps matter
for your question, check `predict(type = "time")` on held-out documents.

## Implementation notes

The sampler follows the paper. Three details differ from a literal reading of
it:

- **Beta density orientation.** The paper's conditional is printed as
  `(1 - t)^(ψ1 - 1) t^(ψ2 - 1)`, which does not match the method-of-moments
  estimates it is paired with. This package uses the standard
  `t^(ψ1 - 1) (1 - t)^(ψ2 - 1)`.
- **Time normalization.** Timestamps are mapped onto `[1e-4, 1 - 1e-4]`
  rather than exactly `[0, 1]`, so log densities stay finite for the earliest
  and latest documents.
- **Bounded Beta densities.** After each method-of-moments update, both Beta
  shapes are kept at or above 1 (same mean, higher concentration). Without
  this, a topic could collapse onto the first and last documents with a
  U-shaped `Beta(≈0, ≈0)` density and stay locked there. The cost is that
  U-shaped topics (common early *and* late, rare in between) cannot be
  represented.

The hyper-parameters `alpha` (default `50 / k`) and `beta` (default `0.1`) are
fixed, as in the paper. The number of topics `k` must be chosen by the user.

## Citation

If you use this package, please cite the original paper:

> Wang, X. and McCallum, A. (2006). Topics over time: a non-Markov
> continuous-time model of topical trends. *Proceedings of the 12th ACM SIGKDD
> International Conference on Knowledge Discovery and Data Mining*, 424–433.
> <https://doi.org/10.1145/1150402.1150450>

## License

[MIT](LICENSE.md) © Jaeseong Choe

## Have an issue?

If you run into a problem, please open an issue on
[GitHub](https://github.com/sorrychoe/topicsovertime/issues).
