#include <Rcpp.h>
using namespace Rcpp;

// Average the four loss derivatives over the selected original observations.
// No intermediate response-gradient matrix is allocated.
// [[Rcpp::export]]
NumericVector dose_gradient_cpp(NumericVector theta, NumericVector z,
                                NumericVector y, IntegerVector index) {
  if (theta.size() != 4 || z.size() != y.size() || index.size() == 0) {
    stop("Invalid dimensions");
  }
  NumericVector g(4);
  double difference = theta[3] - theta[2];
  for (int k = 0; k < index.size(); ++k) {
    int i = index[k] - 1; // Convert R indexing to C++ indexing.
    if (i < 0 || i >= z.size()) stop("Invalid index");
    double eta = theta[0] - theta[1] * z[i];
    double e = std::exp(-std::abs(eta));
    double q = eta >= 0 ? 1.0 / (1.0 + e) : e / (1.0 + e);
    double w = e / ((1.0 + e) * (1.0 + e));
    double residual = theta[2] + difference * q - y[i];
    double scale = 2 * residual / index.size();
    g[0] += scale * difference * w;
    g[1] -= scale * difference * z[i] * w;
    g[2] += scale * (1 - q);
    g[3] += scale * q;
  }
  return g;
}

// Complete mini-batch SGD loop for the controlled benchmark.
// Each column of orders is one precomputed epoch permutation (R indices).
// Sampling, compilation and trajectory monitoring are excluded in BOTH versions.
// [[Rcpp::export]]
NumericVector dose_sgd_cpp(NumericVector start, NumericVector z, NumericVector y,
                           IntegerMatrix orders, NumericVector rates,
                           int batch_size = 1) {
  if (start.size() != 4 || z.size() != y.size() || z.size() == 0 ||
      orders.nrow() != z.size() || orders.ncol() != rates.size() ||
      batch_size < 1 || batch_size > z.size())
    stop("Invalid dimensions or batch size");
  // Clone once: never modify the starting vector supplied by R.
  NumericVector theta = clone(start);
  for (int epoch = 0; epoch < orders.ncol(); ++epoch) {
    for (int first = 0; first < orders.nrow(); first += batch_size) {
      int last = std::min(first + batch_size, orders.nrow());
      // Four scalar accumulators: no batch-sized temporary objects.
      double g0 = 0, g1 = 0, g2 = 0, g3 = 0;
      double difference = theta[3] - theta[2];
      for (int k = first; k < last; ++k) {
        int i = orders(k, epoch) - 1;
        if (i < 0 || i >= z.size()) stop("Invalid index");
        double eta = theta[0] - theta[1] * z[i];
        double e = std::exp(-std::abs(eta));
        double q = eta >= 0 ? 1.0 / (1.0 + e) : e / (1.0 + e);
        double w = e / ((1.0 + e) * (1.0 + e));
        // Average over the ACTUAL batch size, including the final short batch.
        double scale = 2 * (theta[2] + difference * q - y[i]) / (last - first);
        g0 += scale * difference * w;
        g1 -= scale * difference * z[i] * w;
        g2 += scale * (1 - q);
        g3 += scale * q;
      }
      // Update only after all observations in this batch have been evaluated.
      theta[0] -= rates[epoch] * g0;
      theta[1] -= rates[epoch] * g1;
      theta[2] -= rates[epoch] * g2;
      theta[3] -= rates[epoch] * g3;
    }
    checkUserInterrupt();
  }
  return theta;
}
