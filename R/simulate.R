#' Simulate documents from a Topics over Time model
#'
#' Generates a corpus with one timestamp per document, drawn from the mixture
#' of topic Beta distributions weighted by the document's topic proportions
#' (the generative view of Figure 1(b) in Wang and McCallum, 2006). Useful for
#' examples, tests and checking parameter recovery.
#'
#' @param n_docs Number of documents.
#' @param n_topics Number of topics.
#' @param n_words Vocabulary size.
#' @param doc_length Tokens per document.
#' @param psi `n_topics` x 2 matrix of Beta shape parameters. Default: topics
#'   centered at evenly spaced times with concentration 20.
#' @param alpha Dirichlet prior on document-topic proportions.
#' @param beta Dirichlet prior on topic-word distributions.
#'
#' @return A list with `x` (a sparse documents x words count matrix),
#'   `time` (timestamps in (0, 1)), and the true `phi`, `theta` and `psi`.
#' @examples
#' set.seed(1)
#' sim <- simulate_tot(n_docs = 50, n_topics = 2)
#' dim(sim$x)
#' summary(sim$time)
#' @export
simulate_tot <- function(n_docs = 200L, n_topics = 3L, n_words = 100L,
                         doc_length = 50L, psi = NULL, alpha = 0.1,
                         beta = 0.1) {
  if (is.null(psi)) {
    mu <- seq(0.15, 0.85, length.out = n_topics)
    psi <- cbind(mu * 20, (1 - mu) * 20)
  }
  stopifnot(is.matrix(psi), nrow(psi) == n_topics, ncol(psi) == 2)
  phi <- rdirichlet(n_topics, rep(beta, n_words))
  theta <- rdirichlet(n_docs, rep(alpha, n_topics))
  z <- apply(theta, 1, function(p) sample.int(n_topics, 1, prob = p))
  time <- stats::rbeta(n_docs, psi[z, 1], psi[z, 2])
  counts <- vapply(seq_len(n_docs), function(d) {
    stats::rmultinom(1, doc_length, theta[d, ] %*% phi)[, 1]
  }, numeric(n_words))
  x <- methods::as(t(counts), "CsparseMatrix")
  dimnames(x) <- list(paste0("doc", seq_len(n_docs)),
                      paste0("w", seq_len(n_words)))
  list(x = x, time = time, phi = phi, theta = theta, psi = psi)
}

rdirichlet <- function(n, a) {
  g <- matrix(stats::rgamma(n * length(a), a), n, byrow = TRUE)
  g / rowSums(g)
}
