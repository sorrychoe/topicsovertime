// Collapsed Gibbs sampler for Topics over Time (Wang & McCallum, KDD 2006).
//
// Layout: tokens are grouped by document; doc d owns words[doc_ptr[d] .. doc_ptr[d+1]).
// Count matrices are K x D (ndk) and K x V (nkw) so the inner loop over topics
// reads contiguous memory. All randomness comes from R's RNG (set.seed()).

#include <Rcpp.h>
#include <algorithm>
#include <cmath>
#include <vector>
using namespace Rcpp;

static inline double log_beta_pdf(double t, double a, double b) {
  return (a - 1.0) * std::log(t) + (b - 1.0) * std::log1p(-t) - R::lbeta(a, b);
}

// Per-document time factor p(t_d | psi_k)^w, rescaled by its max for stability.
// The timestamp is shared by every token of a document, so this is computed
// once per document instead of once per token.
static void time_factor(double t, double w, const NumericMatrix& psi,
                        std::vector<double>& tf) {
  const int K = tf.size();
  double mx = R_NegInf;
  for (int k = 0; k < K; ++k) {
    tf[k] = w * log_beta_pdf(t, psi(k, 0), psi(k, 1));
    if (tf[k] > mx) mx = tf[k];
  }
  for (int k = 0; k < K; ++k) tf[k] = std::exp(tf[k] - mx);
}

static inline int draw(const std::vector<double>& cum) {
  const int K = cum.size();
  const double u = unif_rand() * cum[K - 1];
  int k = 0;
  while (k < K - 1 && cum[k] < u) ++k;
  return k;
}

// Method-of-moments Beta update (paper, Section 2). Topics with too few
// tokens or zero spread keep their previous parameters.
static void update_psi(const IntegerMatrix& ndk, const NumericVector& time,
                       NumericMatrix& psi) {
  const int K = ndk.nrow(), D = ndk.ncol();
  for (int k = 0; k < K; ++k) {
    double s0 = 0, s1 = 0, s2 = 0;
    for (int d = 0; d < D; ++d) {
      const double c = ndk(k, d), t = time[d];
      s0 += c; s1 += c * t; s2 += c * t * t;
    }
    if (s0 < 2) continue;
    const double m = s1 / s0, v = s2 / s0 - m * m;
    if (v <= 1e-10) continue;
    double common = m * (1 - m) / v - 1;
    if (common <= 0) continue;
    // Keep both shapes >= 1 (same mean, higher concentration). Shapes below 1
    // put infinite density at t = 0 or 1, which lets a topic lock onto the
    // earliest/latest documents (seen on real corpora: Beta(3e-4, 1e-4)).
    common = std::max(common, 1.0 / std::min(m, 1 - m));
    psi(k, 0) = m * common;
    psi(k, 1) = (1 - m) * common;
  }
}

// [[Rcpp::export]]
List tot_gibbs_cpp(IntegerVector words, IntegerVector doc_ptr,
                   NumericVector time, NumericVector weight, int K, int V,
                   double alpha, double beta, int n_iter, bool verbose) {
  const int D = time.size(), N = words.size();
  const double Vbeta = V * beta;
  IntegerVector z(N), nk(K);
  IntegerMatrix ndk(K, D), nkw(K, V);
  NumericMatrix psi(K, 2);
  std::fill(psi.begin(), psi.end(), 1.0);  // Beta(1, 1): first sweep is plain LDA

  for (int d = 0; d < D; ++d) {
    for (int i = doc_ptr[d]; i < doc_ptr[d + 1]; ++i) {
      int k = static_cast<int>(unif_rand() * K);
      if (k == K) k = K - 1;
      z[i] = k; ndk(k, d)++; nkw(k, words[i])++; nk[k]++;
    }
  }

  std::vector<double> tf(K), cum(K);
  for (int it = 0; it < n_iter; ++it) {
    for (int d = 0; d < D; ++d) {
      time_factor(time[d], weight[d], psi, tf);
      for (int i = doc_ptr[d]; i < doc_ptr[d + 1]; ++i) {
        const int w = words[i];
        int k = z[i];
        ndk(k, d)--; nkw(k, w)--; nk[k]--;
        double s = 0;
        for (int j = 0; j < K; ++j) {
          s += (ndk(j, d) + alpha) * (nkw(j, w) + beta) / (nk[j] + Vbeta) * tf[j];
          cum[j] = s;
        }
        k = draw(cum);
        z[i] = k; ndk(k, d)++; nkw(k, w)++; nk[k]++;
      }
    }
    update_psi(ndk, time, psi);
    checkUserInterrupt();
    if (verbose && ((it + 1) % 50 == 0 || it + 1 == n_iter))
      Rcout << "iteration " << it + 1 << " / " << n_iter << "\n";
  }
  return List::create(_["ndk"] = ndk, _["nkw"] = nkw, _["nk"] = nk,
                      _["psi"] = psi);
}

// Fold-in for new documents: phi and psi stay fixed, only z is resampled.
// time[d] = NA drops the time factor (used when predicting the timestamp).
// [[Rcpp::export]]
IntegerMatrix tot_foldin_cpp(IntegerVector words, IntegerVector doc_ptr,
                             NumericVector time, NumericVector weight,
                             NumericMatrix phi, NumericMatrix psi,
                             double alpha, int n_iter) {
  const int K = phi.nrow(), D = time.size(), N = words.size();
  IntegerVector z(N);
  IntegerMatrix ndk(K, D);
  for (int d = 0; d < D; ++d) {
    for (int i = doc_ptr[d]; i < doc_ptr[d + 1]; ++i) {
      int k = static_cast<int>(unif_rand() * K);
      if (k == K) k = K - 1;
      z[i] = k; ndk(k, d)++;
    }
  }
  std::vector<double> tf(K, 1.0), cum(K);
  for (int it = 0; it < n_iter; ++it) {
    for (int d = 0; d < D; ++d) {
      if (ISNAN(time[d])) std::fill(tf.begin(), tf.end(), 1.0);
      else time_factor(time[d], weight[d], psi, tf);
      for (int i = doc_ptr[d]; i < doc_ptr[d + 1]; ++i) {
        const int w = words[i];
        int k = z[i];
        ndk(k, d)--;
        double s = 0;
        for (int j = 0; j < K; ++j) {
          s += (ndk(j, d) + alpha) * phi(j, w) * tf[j];
          cum[j] = s;
        }
        k = draw(cum);
        z[i] = k; ndk(k, d)++;
      }
    }
    checkUserInterrupt();
  }
  return ndk;
}
