fit_sim <- function(seed = 1, ...) {
  set.seed(seed)
  sim <- simulate_tot(n_docs = 80, n_topics = 3)
  list(sim = sim, fit = tot(sim$x, sim$time, k = 3, ...))
}

test_that("log-likelihood is recorded without changing the sampled values", {
  sim <- simulate_tot(n_docs = 40, n_topics = 2)
  set.seed(3); a <- tot(sim$x, sim$time, k = 2, n_iter = 25)
  set.seed(3); b <- tot(sim$x, sim$time, k = 2, n_iter = 25, loglik_every = 0)
  expect_identical(a$phi, b$phi)
  expect_identical(a$psi, b$psi)
  expect_equal(a$loglik$iteration, c(10, 20, 25))
  expect_true(all(is.finite(a$loglik$loglik)))
  expect_null(b$loglik)
  expect_error(plot(b, type = "loglik"), "loglik_every")
})

test_that("averaging over samples keeps valid probabilities", {
  r <- fit_sim(n_iter = 60, burn_in = 30, thin = 5)
  expect_equal(r$fit$n_samples, 6)
  expect_equal(unname(rowSums(r$fit$phi)), rep(1, 3))
  expect_equal(unname(rowSums(r$fit$theta)), rep(1, 80))
  expect_equal(sum(r$fit$topic_size), sum(r$sim$x))
  expect_equal(fit_sim(n_iter = 20)$fit$n_samples, 1)
})

test_that("min_shape bounds the Beta shapes", {
  r <- fit_sim(n_iter = 40, min_shape = 2)
  expect_true(all(r$fit$psi >= 2 - 1e-8))
  expect_no_error(fit_sim(n_iter = 10, min_shape = 0))
})

test_that("invalid new arguments are rejected", {
  sim <- simulate_tot(n_docs = 20, n_topics = 2)
  expect_error(tot(sim$x, sim$time, k = 2, n_iter = 10, burn_in = 10), "burn_in")
  expect_error(tot(sim$x, sim$time, k = 2, thin = 0), "thin")
  expect_error(tot(sim$x, sim$time, k = 2, loglik_every = -1), "loglik_every")
  expect_error(tot(sim$x, sim$time, k = 2, min_shape = -1), "min_shape")
})

test_that("terms() ranks by prob, frex and lift", {
  fit <- fit_sim(n_iter = 30)$fit
  for (type in c("prob", "frex", "lift")) {
    expect_equal(dim(terms(fit, 4, type = type)), c(4, 3))
  }
  expect_equal(dim(terms(fit, 1)), c(1, 3))
  top <- terms(fit, 1, type = "lift")[1, "Topic1"]
  lift <- fit$phi[1, ] / colSums(fit$phi * fit$topic_size / sum(fit$topic_size))
  expect_equal(unname(top), names(which.max(lift)))
  expect_error(terms(fit, 3, type = "frex", frex_weight = 2), "frex_weight")
})

test_that("top_docs() returns the documents with the largest topic share", {
  fit <- fit_sim(n_iter = 30)$fit
  td <- top_docs(fit, 2)
  expect_equal(dim(td), c(2, 3))
  expect_equal(unname(td[1, "Topic2"]), rownames(fit$theta)[which.max(fit$theta[, 2])])
})

test_that("tidy() and augment() return long data frames", {
  skip_if_not_installed("generics")
  r <- fit_sim(n_iter = 30)
  fit <- r$fit
  phi <- generics::tidy(fit)
  expect_named(phi, c("topic", "term", "phi"))
  expect_equal(nrow(phi), 3 * 100)
  expect_equal(as.vector(tapply(phi$phi, phi$topic, sum)), rep(1, 3))
  expect_equal(nrow(generics::tidy(fit, matrix = "theta")), 80 * 3)
  psi <- generics::tidy(fit, matrix = "psi")
  expect_named(psi, c("topic", "shape1", "shape2", "mean_time"))
  expect_true(all(psi$mean_time >= min(r$sim$time) & psi$mean_time <= max(r$sim$time)))

  aug <- generics::augment(fit, r$sim$x)
  expect_named(aug, c("document", "term", "count", ".topic"))
  expect_equal(nrow(aug), Matrix::nnzero(r$sim$x))
  expect_true(all(aug$.topic %in% rownames(fit$phi)))
  expect_error(generics::augment(fit, r$sim$x[1:5, ]), "fitted on")
})

test_that("tidy(matrix = 'psi') keeps Date timestamps", {
  skip_if_not_installed("generics")
  set.seed(4)
  sim <- simulate_tot(n_docs = 40, n_topics = 2)
  dates <- as.Date("2000-01-01") + round(sim$time * 3650)
  fit <- tot(sim$x, dates, k = 2, n_iter = 10)
  expect_s3_class(generics::tidy(fit, matrix = "psi")$mean_time, "Date")
})

test_that("autoplot() builds ggplot objects", {
  skip_if_not_installed("ggplot2")
  fit <- fit_sim(n_iter = 20)$fit
  for (type in c("trend", "density", "loglik")) {
    expect_s3_class(ggplot2::autoplot(fit, type = type), "ggplot")
  }
})

test_that("citation() includes the original paper", {
  path <- system.file("CITATION", package = "topicsovertime")
  skip_if(path == "")
  cit <- utils::readCitationFile(path, meta = list(Version = "0.0.0"))
  expect_length(cit, 2)
  expect_match(format(cit[1], style = "text"), "McCallum")
})
