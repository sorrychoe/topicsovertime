#' Tidy a Topics over Time model
#'
#' Methods for the [generics::tidy()] and [generics::augment()] generics,
#' which are also exported by broom and tidytext. They are registered when
#' the generics package is loaded.
#'
#' @param x A fitted `"tot"` model.
#' @param matrix Which estimates to return: `"phi"` (topic-word
#'   probabilities), `"theta"` (document-topic proportions) or `"psi"`
#'   (topic time distributions).
#' @param data For `augment()`, the document-term matrix the model was fitted
#'   on (any type accepted by [tot()]).
#' @param ... Unused.
#'
#' @return A data frame.
#'   * `tidy(matrix = "phi")`: one row per topic and word, with columns
#'     `topic`, `term` and `phi`.
#'   * `tidy(matrix = "theta")`: one row per document and topic, with columns
#'     `document`, `topic` and `theta`.
#'   * `tidy(matrix = "psi")`: one row per topic, with columns `topic`,
#'     `shape1`, `shape2` (on the normalized time scale) and `mean_time`, the
#'     mean of the Beta distribution on the scale of the training `time`.
#'   * `augment()`: one row per non-zero document-word cell of `data`, with
#'     columns `document`, `term`, `count` and `.topic`, the topic most likely
#'     to have generated that word in that document
#'     (\eqn{\arg\max_k \theta_{dk} \phi_{kv}}{argmax_k theta_dk phi_kv}).
#'     Words not in the model's vocabulary are dropped.
#' @examples
#' set.seed(1)
#' sim <- simulate_tot(n_docs = 100, n_topics = 3)
#' fit <- tot(sim$x, sim$time, k = 3, n_iter = 50)
#' if (requireNamespace("generics", quietly = TRUE)) {
#'   head(generics::tidy(fit))
#'   generics::tidy(fit, matrix = "psi")
#'   head(generics::augment(fit, sim$x))
#' }
#' @name tidy.tot
#' @exportS3Method generics::tidy
tidy.tot <- function(x, matrix = c("phi", "theta", "psi"), ...) {
  matrix <- match.arg(matrix)
  topics <- rownames(x$phi)
  if (matrix == "phi") {
    vocab <- colnames(x$phi)
    if (is.null(vocab)) vocab <- as.character(seq_len(ncol(x$phi)))
    return(data.frame(topic = rep(topics, ncol(x$phi)),
                      term = rep(vocab, each = nrow(x$phi)),
                      phi = as.vector(x$phi)))
  }
  if (matrix == "theta") {
    docs <- rownames(x$theta)
    if (is.null(docs)) docs <- as.character(seq_len(nrow(x$theta)))
    return(data.frame(document = rep(docs, ncol(x$theta)),
                      topic = rep(topics, each = nrow(x$theta)),
                      theta = as.vector(x$theta)))
  }
  data.frame(topic = topics, shape1 = x$psi[, 1], shape2 = x$psi[, 2],
             mean_time = from_unit(x, x$psi[, 1] / rowSums(x$psi)),
             row.names = NULL)
}

#' @rdname tidy.tot
#' @exportS3Method generics::augment
augment.tot <- function(x, data, ...) {
  dtm <- align_vocab(as_triplet(data), x)
  if (dtm$D != nrow(x$theta)) {
    stop("`data` must be the document-term matrix the model was fitted on.",
         call. = FALSE)
  }
  keep <- dtm$v > 0
  i <- dtm$i[keep] + 1L
  j <- dtm$j[keep] + 1L
  topic <- vapply(seq_along(i), function(r) {
    which.max(x$theta[i[r], ] * x$phi[, j[r]])
  }, integer(1))
  docs <- rownames(x$theta)
  if (is.null(docs)) docs <- as.character(seq_len(nrow(x$theta)))
  vocab <- colnames(x$phi)
  if (is.null(vocab)) vocab <- as.character(seq_len(ncol(x$phi)))
  out <- data.frame(document = docs[i], term = vocab[j], count = dtm$v[keep],
                    .topic = rownames(x$phi)[topic])
  out <- out[order(i, j), , drop = FALSE]
  rownames(out) <- NULL
  out
}

#' Plot a Topics over Time model with ggplot2
#'
#' A ggplot2 version of [plot.tot()]. Registered as a method for
#' [ggplot2::autoplot()] when ggplot2 is loaded.
#'
#' @inheritParams plot.tot
#' @param object A fitted `"tot"` model.
#' @return A ggplot object.
#' @examples
#' set.seed(1)
#' sim <- simulate_tot(n_docs = 100, n_topics = 3)
#' fit <- tot(sim$x, sim$time, k = 3, n_iter = 50)
#' if (requireNamespace("ggplot2", quietly = TRUE)) {
#'   ggplot2::autoplot(fit)
#'   ggplot2::autoplot(fit, type = "loglik")
#' }
#' @exportS3Method ggplot2::autoplot
autoplot.tot <- function(object, type = c("trend", "density", "loglik"),
                         topics = NULL, n = 200L, ...) {
  type <- match.arg(type)
  if (type == "loglik") {
    return(ggplot2::ggplot(loglik_data(object),
                           ggplot2::aes(.data$iteration, .data$loglik)) +
             ggplot2::geom_line() +
             ggplot2::labs(x = "Iteration", y = "Log-likelihood"))
  }
  times <- time_grid(object, n)
  y <- time_curves(object, type, times, topics)
  df <- data.frame(time = rep(times, ncol(y)),
                   topic = factor(rep(colnames(y), each = nrow(y)),
                                  levels = colnames(y)),
                   value = as.vector(y))
  ylab <- if (type == "trend") "P(topic | time)" else "Density (normalized time)"
  ggplot2::ggplot(df, ggplot2::aes(.data$time, .data$value,
                                   colour = .data$topic)) +
    ggplot2::geom_line() +
    ggplot2::labs(x = "Time", y = ylab, colour = "Topic")
}
