#' Predict topics or timestamps for new documents
#'
#' Runs Gibbs sampling on new documents with the fitted topic-word and
#' topic-time distributions held fixed ("fold-in").
#'
#' With `type = "time"`, topics are sampled from the words alone and the
#' timestamp is chosen to maximize \eqn{\sum_i \log p(t \mid \psi_{z_i})}{sum_i
#' log p(t | psi_z_i)} over a grid spanning the training time range, as in
#' Section 2 of Wang and McCallum (2006).
#'
#' @param object A fitted `"tot"` model.
#' @param newdata A document-term matrix (see [tot()]). Columns are matched to
#'   the training vocabulary by name when both have column names; unknown
#'   words are dropped. Without names, columns must be in training order.
#' @param type `"topics"` for document-topic proportions, or `"time"` for
#'   predicted timestamps.
#' @param time Optional timestamps of `newdata`, used by `type = "topics"` so
#'   that time also informs the topic assignments.
#' @param n_iter Number of fold-in sweeps.
#' @param n_grid Number of grid points for `type = "time"`.
#' @param ... Unused.
#'
#' @return For `type = "topics"`, a documents x topics matrix whose rows sum
#'   to 1. For `type = "time"`, a vector of timestamps on the scale of the
#'   training `time` (`NA` for documents without known words).
#' @examples
#' set.seed(1)
#' sim <- simulate_tot(n_docs = 100, n_topics = 3)
#' fit <- tot(sim$x, sim$time, k = 3, n_iter = 50)
#' head(predict(fit, sim$x[1:5, ], type = "topics"))
#' predict(fit, sim$x[1:5, ], type = "time")
#' @export
predict.tot <- function(object, newdata, type = c("topics", "time"),
                        time = NULL, n_iter = 100L, n_grid = 512L, ...) {
  type <- match.arg(type)
  check_count(n_iter, "n_iter", min = 1)
  dtm <- align_vocab(as_triplet(newdata), object)
  tok <- tokens_by_doc(dtm)

  if (type == "topics" && !is.null(time)) {
    if (length(time) != dtm$D) {
      stop("`time` must have one value per document of `newdata`.", call. = FALSE)
    }
    tu <- to_unit(object$time, time_scale_values(time))
  } else {
    tu <- rep(NA_real_, dtm$D)
  }
  ndk <- t(tot_foldin_cpp(tok$words, tok$doc_ptr, tu,
                          doc_weight(object$time_weight, tok$n), object$phi,
                          object$psi, object$alpha, as.integer(n_iter)))

  if (type == "topics") {
    theta <- (ndk + object$alpha) / (tok$n + object$k * object$alpha)
    dimnames(theta) <- list(dtm$docs, rownames(object$phi))
    return(theta)
  }
  # Uses the last fold-in sample; average several samples if predictions
  # look noisy on short documents.
  grid <- seq(.eps, 1 - .eps, length.out = n_grid)
  loglik <- ndk %*% t(log_beta_grid(object, grid))
  out <- from_unit(object, grid[max.col(loglik, ties.method = "first")])
  out[tok$n == 0] <- NA
  names(out) <- dtm$docs
  out
}

#' Topic proportions over time
#'
#' Computes \eqn{P(z \mid t) \propto p(t \mid \psi_z) P(z)}{P(z | t) ~ p(t |
#' psi_z) P(z)}, the expected topic mixture at each time point, with
#' \eqn{P(z)} estimated from the fitted topic sizes. Times outside the
#' training range are clamped to its ends.
#'
#' @param object A fitted `"tot"` model.
#' @param times Time points on the scale of the training `time`. Default: 200
#'   evenly spaced points over the training range.
#'
#' @return A matrix with one row per time point and one column per topic;
#'   rows sum to 1. The time points are stored in attribute `"times"`.
#' @examples
#' set.seed(1)
#' sim <- simulate_tot(n_docs = 100, n_topics = 3)
#' fit <- tot(sim$x, sim$time, k = 3, n_iter = 50)
#' round(topic_time(fit, c(0.1, 0.5, 0.9)), 3)
#' @export
topic_time <- function(object, times = NULL) {
  if (is.null(times)) times <- time_grid(object)
  lp <- log_beta_grid(object, to_unit(object$time, time_scale_values(times)))
  lp <- sweep(lp, 2, log(object$topic_size / sum(object$topic_size)), "+")
  p <- exp(lp - apply(lp, 1, max))
  p <- p / rowSums(p)
  colnames(p) <- rownames(object$phi)
  attr(p, "times") <- times
  p
}

#' Top words per topic
#'
#' @param x A fitted `"tot"` model.
#' @param n Number of words per topic.
#' @param type How to rank words within a topic:
#'   \describe{
#'     \item{`"prob"`}{Topic-word probability \eqn{\phi_{kv}}{phi_kv}.}
#'     \item{`"frex"`}{FREX: the weighted harmonic mean of a word's rank in
#'       frequency (\eqn{\phi_{kv}}{phi_kv}) and in exclusivity
#'       (\eqn{\phi_{kv} / \sum_j \phi_{jv}}{phi_kv / sum_j phi_jv}), as in
#'       Bischof and Airoldi (2012). Favors words that are both common in the
#'       topic and rare in others.}
#'     \item{`"lift"`}{\eqn{\phi_{kv} / p(v)}{phi_kv / p(v)}, where
#'       \eqn{p(v)} is the word's overall probability under the model. Favors
#'       rare, topic-specific words.}
#'   }
#' @param frex_weight Weight on exclusivity for `type = "frex"`, in `[0, 1]`.
#' @param ... Unused.
#' @return A character matrix with `n` rows and one column per topic.
#' @references Bischof, J. and Airoldi, E. (2012). Summarizing topical
#'   content with word frequency and exclusivity. *Proceedings of the 29th
#'   International Conference on Machine Learning*, 201--208.
#' @examples
#' set.seed(1)
#' sim <- simulate_tot(n_docs = 100, n_topics = 3)
#' fit <- tot(sim$x, sim$time, k = 3, n_iter = 50)
#' terms(fit, 5)
#' terms(fit, 5, type = "frex")
#' @importFrom stats terms
#' @export
terms.tot <- function(x, n = 10L, type = c("prob", "frex", "lift"),
                      frex_weight = 0.5, ...) {
  type <- match.arg(type)
  check_count(n, "n", min = 1)
  vocab <- colnames(x$phi)
  if (is.null(vocab)) vocab <- as.character(seq_len(ncol(x$phi)))
  score <- word_score(x, type, frex_weight)
  n <- min(n, length(vocab))
  out <- apply(score, 1, function(s) vocab[order(s, decreasing = TRUE)[seq_len(n)]])
  if (n == 1) out <- matrix(out, nrow = 1, dimnames = list(NULL, rownames(x$phi)))
  out
}

#' Most representative documents per topic
#'
#' @param object A fitted `"tot"` model.
#' @param n Number of documents per topic.
#' @return A matrix with `n` rows and one column per topic, holding the names
#'   (or row numbers, if the training matrix had no row names) of the
#'   documents with the highest proportion of each topic.
#' @examples
#' set.seed(1)
#' sim <- simulate_tot(n_docs = 100, n_topics = 3)
#' fit <- tot(sim$x, sim$time, k = 3, n_iter = 50)
#' top_docs(fit, 3)
#' @export
top_docs <- function(object, n = 5L) {
  check_count(n, "n", min = 1)
  docs <- rownames(object$theta)
  if (is.null(docs)) docs <- as.character(seq_len(nrow(object$theta)))
  n <- min(n, length(docs))
  out <- apply(object$theta, 2, function(p) docs[order(p, decreasing = TRUE)[seq_len(n)]])
  if (n == 1) out <- matrix(out, nrow = 1, dimnames = list(NULL, colnames(object$theta)))
  out
}

#' Plot topic trends, topic time densities or the log-likelihood
#'
#' @param x A fitted `"tot"` model.
#' @param type `"trend"` plots [topic_time()] (topic proportions over time);
#'   `"density"` plots each topic's Beta density over time; `"loglik"` plots
#'   the log-likelihood recorded during fitting, to check convergence.
#' @param topics Topics to plot (indices or names). Default: all.
#' @param n Number of time points.
#' @param legend Draw a legend.
#' @param ... Passed to [graphics::matplot()] or [graphics::plot()].
#' @return The plotted data, invisibly.
#' @seealso [autoplot.tot()] for a ggplot2 version.
#' @examples
#' set.seed(1)
#' sim <- simulate_tot(n_docs = 100, n_topics = 3)
#' fit <- tot(sim$x, sim$time, k = 3, n_iter = 50)
#' plot(fit)
#' plot(fit, type = "density")
#' plot(fit, type = "loglik")
#' @export
plot.tot <- function(x, type = c("trend", "density", "loglik"), topics = NULL,
                     n = 200L, legend = TRUE, ...) {
  type <- match.arg(type)
  if (type == "loglik") {
    ll <- loglik_data(x)
    graphics::plot(ll$iteration, ll$loglik, type = "l", xlab = "Iteration",
                   ylab = "Log-likelihood", ...)
    return(invisible(ll))
  }
  times <- time_grid(x, n)
  y <- time_curves(x, type, times, topics)
  ylab <- if (type == "trend") "P(topic | time)" else "Density (normalized time)"
  graphics::matplot(as.numeric(times), y, type = "l", lty = 1, xaxt = "n",
                    xlab = "Time", ylab = ylab, col = seq_len(ncol(y)), ...)
  if (inherits(times, "Date")) {
    graphics::axis.Date(1, times)
  } else if (inherits(times, "POSIXct")) {
    graphics::axis.POSIXct(1, times)
  } else {
    graphics::axis(1)
  }
  if (legend) {
    graphics::legend("topright", colnames(y), col = seq_len(ncol(y)), lty = 1,
                     bty = "n", cex = 0.8)
  }
  invisible(y)
}

# ---- helpers --------------------------------------------------------------

# topics x words matrix of scores used to rank words
word_score <- function(x, type, frex_weight) {
  phi <- x$phi
  if (type == "prob") return(phi)
  if (type == "lift") {
    pv <- colSums(phi * (x$topic_size / sum(x$topic_size)))
    return(sweep(phi, 2, pv, "/"))
  }
  if (!is.numeric(frex_weight) || length(frex_weight) != 1 ||
      !is.finite(frex_weight) || frex_weight < 0 || frex_weight > 1) {
    stop("`frex_weight` must be a number in [0, 1].", call. = FALSE)
  }
  excl <- sweep(phi, 2, colSums(phi), "/")
  ecdf_rank <- function(m) t(apply(m, 1, function(v) rank(v, ties.method = "max") / length(v)))
  1 / (frex_weight / ecdf_rank(excl) + (1 - frex_weight) / ecdf_rank(phi))
}

# grid points x topics matrix for type = "trend" or "density"
time_curves <- function(x, type, times, topics = NULL) {
  y <- if (type == "trend") {
    topic_time(x, times)
  } else {
    exp(log_beta_grid(x, to_unit(x$time, as.numeric(times))))
  }
  colnames(y) <- rownames(x$phi)
  if (!is.null(topics)) y <- y[, topics, drop = FALSE]
  y
}

loglik_data <- function(x) {
  if (is.null(x$loglik)) {
    stop("No log-likelihood was recorded; refit with `loglik_every` > 0.",
         call. = FALSE)
  }
  x$loglik
}

# grid points x topics matrix of Beta log-densities on the unit scale
log_beta_grid <- function(object, u) {
  vapply(seq_len(object$k), function(k) {
    stats::dbeta(u, object$psi[k, 1], object$psi[k, 2], log = TRUE)
  }, numeric(length(u)))
}

time_grid <- function(object, n = 200L) {
  from_unit(object, seq(.eps, 1 - .eps, length.out = n))
}

time_scale_values <- function(time) {
  if (is.factor(time) || is.character(time) || is.logical(time)) {
    stop("`time` must be numeric, Date or POSIXct.", call. = FALSE)
  }
  as.numeric(time)
}

align_vocab <- function(dtm, object) {
  vocab <- colnames(object$phi)
  if (!is.null(dtm$vocab) && !is.null(vocab)) {
    j <- match(dtm$vocab, vocab)[dtm$j + 1L]
    keep <- !is.na(j)
    dtm$i <- dtm$i[keep]; dtm$v <- dtm$v[keep]; dtm$j <- j[keep] - 1L
  } else if (dtm$V != ncol(object$phi)) {
    stop("`newdata` has ", dtm$V, " columns but the model has ",
         ncol(object$phi), " words; supply column names to match them.",
         call. = FALSE)
  }
  dtm$V <- ncol(object$phi)
  dtm
}
