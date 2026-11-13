# Match fitted topics to true topics by word-distribution correlation.
match_topics <- function(fit, sim) {
  apply(sim$phi, 1, function(p) which.max(cor(p, t(fit$phi))))
}

beta_mean <- function(psi) psi[, 1] / (psi[, 1] + psi[, 2])

test_that("recovers topic words and time locations on simulated data", {
  set.seed(42)
  sim <- simulate_tot(n_docs = 300, n_topics = 3, n_words = 150,
                      doc_length = 60)
  fit <- tot(sim$x, sim$time, k = 3, n_iter = 150)
  m <- match_topics(fit, sim)
  expect_setequal(m, 1:3)
  # psi lives on the normalized scale; map the true means onto it
  rng <- range(sim$time)
  true_mu <- (beta_mean(sim$psi) - rng[1]) / diff(rng)
  expect_equal(unname(beta_mean(fit$psi)[m]), true_mu, tolerance = 0.1)
})

test_that("Beta orientation: early topic gets an early density", {
  set.seed(1)
  sim <- simulate_tot(n_docs = 200, n_topics = 2,
                      psi = rbind(c(2, 18), c(18, 2)))
  fit <- tot(sim$x, sim$time, k = 2, n_iter = 100)
  m <- match_topics(fit, sim)
  expect_lt(beta_mean(fit$psi)[m[1]], 0.5)
  expect_gt(beta_mean(fit$psi)[m[2]], 0.5)
  tt <- topic_time(fit, range(sim$time))
  expect_gt(tt[1, m[1]], tt[1, m[2]])
  expect_gt(tt[2, m[2]], tt[2, m[1]])
})

test_that("output structure and probabilities are valid", {
  set.seed(2)
  sim <- simulate_tot(n_docs = 50, n_topics = 2)
  fit <- tot(sim$x, sim$time, k = 2, n_iter = 20)
  expect_s3_class(fit, "tot")
  expect_equal(dim(fit$phi), c(2, 100))
  expect_equal(dim(fit$theta), c(50, 2))
  expect_equal(unname(rowSums(fit$phi)), c(1, 1))
  expect_equal(unname(rowSums(fit$theta)), rep(1, 50))
  expect_equal(sum(fit$topic_size), sum(sim$x))
  expect_true(all(fit$psi > 0))
  expect_equal(dim(terms(fit, 4)), c(4, 2))
  expect_output(print(fit), "Topics over Time")
})

test_that("set.seed makes fits reproducible", {
  sim <- simulate_tot(n_docs = 40, n_topics = 2)
  set.seed(7); a <- tot(sim$x, sim$time, k = 2, n_iter = 10)
  set.seed(7); b <- tot(sim$x, sim$time, k = 2, n_iter = 10)
  expect_identical(a$phi, b$phi)
  expect_identical(a$psi, b$psi)
})

test_that("matrix, sparse, dfm and DocumentTermMatrix inputs agree", {
  set.seed(3)
  sim <- simulate_tot(n_docs = 40, n_topics = 2)
  fit_with <- function(x) {
    set.seed(9)
    tot(x, sim$time, k = 2, n_iter = 10)$phi
  }
  ref <- fit_with(sim$x)
  expect_equal(fit_with(as.matrix(sim$x)), ref)
  tm_like <- structure(
    list(i = Matrix::summary(sim$x)$i, j = Matrix::summary(sim$x)$j,
         v = Matrix::summary(sim$x)$x, nrow = nrow(sim$x),
         ncol = ncol(sim$x), dimnames = dimnames(sim$x)),
    class = c("DocumentTermMatrix", "simple_triplet_matrix")
  )
  expect_equal(fit_with(tm_like), ref)
  skip_if_not_installed("quanteda")
  expect_equal(fit_with(quanteda::as.dfm(sim$x)), ref, ignore_attr = TRUE)
})

test_that("Date timestamps round-trip through predictions", {
  set.seed(4)
  sim <- simulate_tot(n_docs = 60, n_topics = 2)
  dates <- as.Date("2000-01-01") + round(sim$time * 3650)
  fit <- tot(sim$x, dates, k = 2, n_iter = 20)
  pt <- predict(fit, sim$x[1:3, ], type = "time", n_iter = 10)
  expect_s3_class(pt, "Date")
  expect_true(all(pt >= min(dates) & pt <= max(dates)))
  expect_s3_class(attr(topic_time(fit), "times"), "Date")
  expect_equal(colnames(topic_time(fit)), rownames(fit$phi))
})

test_that("predicted times track true times", {
  set.seed(5)
  sim <- simulate_tot(n_docs = 300, n_topics = 3, doc_length = 80)
  fit <- tot(sim$x, sim$time, k = 3, n_iter = 100)
  pt <- predict(fit, sim$x, type = "time", n_iter = 30)
  expect_gt(cor(pt, sim$time), 0.6)
})

test_that("predict(type = 'topics') aligns vocabulary by name", {
  set.seed(6)
  sim <- simulate_tot(n_docs = 50, n_topics = 2)
  fit <- tot(sim$x, sim$time, k = 2, n_iter = 20)
  shuffled <- sim$x[1:5, rev(seq_len(ncol(sim$x)))]
  set.seed(1); a <- predict(fit, sim$x[1:5, ], n_iter = 20)
  set.seed(1); b <- predict(fit, shuffled, n_iter = 20)
  expect_equal(dim(a), c(5, 2))
  expect_equal(unname(rowSums(a)), rep(1, 5))
  # same tokens, different column order: same distribution of results
  expect_equal(colMeans(a), colMeans(b), tolerance = 0.3)
  expect_silent(predict(fit, sim$x[1:2, ], time = sim$time[1:2]))
})

test_that("invalid input is rejected", {
  sim <- simulate_tot(n_docs = 20, n_topics = 2)
  expect_error(tot(sim$x, sim$time[-1], k = 2), "one value per document")
  expect_error(tot(sim$x, rep(1, 20), k = 2), "constant")
  expect_error(tot(sim$x, sim$time, k = 1), "`k`")
  expect_error(tot(sim$x * 0.5, sim$time, k = 2), "integer counts")
  expect_error(tot(sim$x, as.character(sim$time), k = 2), "numeric, Date")
  expect_error(tot(sim$x, sim$time, k = 2, time_weight = "x"), "time_weight")
  expect_error(tot(list(1), sim$time, k = 2), "document-term matrix")
})
