#' Fit a Topics over Time model
#'
#' Fits the Topics over Time (TOT) model of Wang and McCallum (2006) by
#' collapsed Gibbs sampling. Each topic has a multinomial distribution over
#' words and a Beta distribution over (normalized) time, so topics are shaped
#' jointly by word co-occurrence and by when documents were written.
#'
#' Timestamps are mapped linearly onto the open interval (0, 1) using the
#' range of `time`. The Beta parameters are re-estimated by the method of
#' moments after every sweep, as in the paper, with both shapes kept at or
#' above 1 (same mean, higher concentration) so densities stay bounded and
#' topics cannot lock onto the earliest or latest documents.
#'
#' @param x A document-term matrix of counts with documents in rows: a base
#'   `matrix`, a `Matrix` sparse matrix, a quanteda `dfm`, or a tm
#'   `DocumentTermMatrix`.
#' @param time Document timestamps, one per row of `x`: numeric, `Date` or
#'   `POSIXct`.
#' @param k Number of topics.
#' @param alpha Symmetric Dirichlet prior on document-topic proportions.
#'   Default `50 / k`, as in the paper.
#' @param beta Symmetric Dirichlet prior on topic-word distributions.
#' @param n_iter Number of Gibbs sweeps.
#' @param time_weight Exponent on the time likelihood, balancing time against
#'   words. A number (`1` gives the sampler exactly as in the paper), or
#'   `"inverse_length"` for `1 / N_d`, the alternative the paper describes,
#'   which treats each document's timestamp as one observation. The paper
#'   reports that results are sensitive to this setting.
#' @param verbose Print progress every 50 iterations.
#'
#' @return An object of class `"tot"`, a list with
#'   \describe{
#'     \item{phi}{`k` x vocabulary matrix of topic-word probabilities.}
#'     \item{theta}{documents x `k` matrix of document-topic proportions.}
#'     \item{psi}{`k` x 2 matrix of Beta shape parameters on the normalized
#'       time scale.}
#'     \item{topic_size}{Number of tokens assigned to each topic.}
#'   }
#'   plus the settings used and the information needed to map times back to
#'   the original scale.
#'
#' @references Wang, X. and McCallum, A. (2006). Topics over time: a
#'   non-Markov continuous-time model of topical trends. *Proceedings of the
#'   12th ACM SIGKDD International Conference on Knowledge Discovery and Data
#'   Mining*, 424--433. \doi{10.1145/1150402.1150450}
#'
#' @seealso [predict.tot()], [topic_time()], [terms.tot()], [plot.tot()]
#' @examples
#' set.seed(1)
#' sim <- simulate_tot(n_docs = 100, n_topics = 3)
#' fit <- tot(sim$x, sim$time, k = 3, n_iter = 50)
#' fit
#' terms(fit, 5)
#' @export
tot <- function(x, time, k, alpha = 50 / k, beta = 0.1, n_iter = 500L,
                time_weight = 1, verbose = FALSE) {
  dtm <- as_triplet(x)
  check_count(k, "k", min = 2)
  check_count(n_iter, "n_iter", min = 1)
  check_positive(alpha, "alpha")
  check_positive(beta, "beta")
  if (length(time) != dtm$D) {
    stop("`time` must have one value per document (row of `x`).", call. = FALSE)
  }
  tn <- time_scale(time)
  tok <- tokens_by_doc(dtm)
  weight <- doc_weight(time_weight, tok$n)

  res <- tot_gibbs_cpp(tok$words, tok$doc_ptr, to_unit(tn, tn$value), weight,
                       as.integer(k), dtm$V, alpha, beta, as.integer(n_iter),
                       isTRUE(verbose))

  topics <- paste0("Topic", seq_len(k))
  ndk <- t(res$ndk)
  theta <- (ndk + alpha) / (tok$n + k * alpha)
  phi <- (res$nkw + beta) / (res$nk + dtm$V * beta)
  psi <- res$psi
  dimnames(theta) <- list(dtm$docs, topics)
  dimnames(phi) <- list(topics, dtm$vocab)
  dimnames(psi) <- list(topics, c("shape1", "shape2"))

  structure(
    list(phi = phi, theta = theta, psi = psi,
         topic_size = stats::setNames(res$nk, topics),
         alpha = alpha, beta = beta, k = as.integer(k),
         n_iter = as.integer(n_iter), time_weight = time_weight,
         time = tn[c("range", "attr")], call = match.call()),
    class = "tot"
  )
}

#' @export
print.tot <- function(x, ...) {
  rng <- from_unit(x, c(0, 1) * (1 - 2 * .eps) + .eps)
  cat("Topics over Time model\n",
      "  topics: ", x$k, ", documents: ", nrow(x$theta),
      ", vocabulary: ", ncol(x$phi), "\n",
      "  iterations: ", x$n_iter, ", alpha: ", format(x$alpha),
      ", beta: ", format(x$beta), ", time_weight: ", format(x$time_weight), "\n",
      "  time range: ", format(rng[1]), " to ", format(rng[2]), "\n", sep = "")
  invisible(x)
}

# ---- input handling -------------------------------------------------------

# Convert any supported document-term matrix to 0-based triplets.
as_triplet <- function(x) {
  if (inherits(x, "simple_triplet_matrix")) {  # tm / slam, no dependency needed
    if (inherits(x, "TermDocumentMatrix")) {
      x <- list(i = x$j, j = x$i, v = x$v, nrow = x$ncol, ncol = x$nrow,
                dimnames = rev(x$dimnames))
    }
    out <- list(i = x$i - 1L, j = x$j - 1L, v = x$v, D = x$nrow, V = x$ncol,
                docs = x$dimnames[[1]], vocab = x$dimnames[[2]])
  } else {
    if (!is.matrix(x) && !methods::is(x, "Matrix")) {
      stop("`x` must be a document-term matrix (matrix, Matrix, dfm or ",
           "DocumentTermMatrix).", call. = FALSE)
    }
    m <- methods::as(methods::as(methods::as(x, "dMatrix"), "generalMatrix"),
                     "TsparseMatrix")
    out <- list(i = m@i, j = m@j, v = m@x, D = nrow(m), V = ncol(m),
                docs = rownames(x), vocab = colnames(x))
  }
  if (anyNA(out$v) || any(out$v < 0) || any(out$v != round(out$v))) {
    stop("`x` must contain non-negative integer counts.", call. = FALSE)
  }
  if (sum(out$v) == 0) stop("`x` contains no tokens.", call. = FALSE)
  out
}

# Expand counts into one entry per token, grouped by document.
tokens_by_doc <- function(dtm) {
  keep <- dtm$v > 0
  i <- dtm$i[keep]; j <- dtm$j[keep]; v <- as.integer(dtm$v[keep])
  o <- order(i)
  n <- integer(dtm$D)
  s <- rowsum(v, i)
  n[as.integer(rownames(s)) + 1L] <- as.integer(s)
  list(words = rep.int(as.integer(j[o]), v[o]),
       doc_ptr = as.integer(c(0, cumsum(n))), n = n)
}

doc_weight <- function(time_weight, n) {
  if (identical(time_weight, "inverse_length")) return(1 / pmax(n, 1))
  if (!is.numeric(time_weight) || length(time_weight) != 1 ||
      !is.finite(time_weight) || time_weight < 0) {
    stop("`time_weight` must be a non-negative number or \"inverse_length\".",
         call. = FALSE)
  }
  rep(time_weight, length(n))
}

check_count <- function(x, name, min) {
  if (!is.numeric(x) || length(x) != 1 || !is.finite(x) || x < min ||
      x != round(x)) {
    stop("`", name, "` must be a whole number >= ", min, ".", call. = FALSE)
  }
}

check_positive <- function(x, name) {
  if (!is.numeric(x) || length(x) != 1 || !is.finite(x) || x <= 0) {
    stop("`", name, "` must be a positive number.", call. = FALSE)
  }
}

# ---- time scale -----------------------------------------------------------

# Times are mapped onto [eps, 1 - eps] so Beta log-densities stay finite at
# the earliest and latest documents.
# eps is fixed; expose it if users need control over endpoint behaviour.
.eps <- 1e-4

time_scale <- function(time) {
  if (is.factor(time) || is.character(time) || is.logical(time)) {
    stop("`time` must be numeric, Date or POSIXct.", call. = FALSE)
  }
  value <- as.numeric(time)
  if (!all(is.finite(value))) stop("`time` must not contain missing or ",
                                   "infinite values.", call. = FALSE)
  rng <- range(value)
  if (rng[1] == rng[2]) stop("`time` must not be constant.", call. = FALSE)
  a <- attributes(time)
  list(value = value, range = rng, attr = a[setdiff(names(a), "names")])
}

to_unit <- function(scale, value) {
  u <- (as.numeric(value) - scale$range[1]) / diff(scale$range)
  .eps + (1 - 2 * .eps) * pmin(pmax(u, 0), 1)
}

from_unit <- function(object, u) {
  rng <- object$time$range
  v <- rng[1] + (u - .eps) / (1 - 2 * .eps) * diff(rng)
  attributes(v) <- object$time$attr
  v
}
