#include "../include/Algorithms.h"
#include "../include/Operations.h"
#include <complex>
#include <stdexcept>
#include <vector>

#ifdef HAS_ACCELERATE
#include <Accelerate/Accelerate.h>
#endif

// Threshold to switch to naive multiplication
const int BLOCK_SIZE = 16;
const int RECURSION_BASE_SIZE = 1;

template <typename T>
Matrix<T> naive_multiply(const Matrix<T> &A, const Matrix<T> &B) {
  if (A.cols != B.rows) {
    throw std::invalid_argument("Matrix dimensions mismatch");
  }
  Matrix<T> C(A.rows, B.cols);

  // Initialize C to 0
  for (int i = 0; i < C.rows * C.cols; ++i)
    C.data[i] = T(0);

  for (int i = 0; i < A.rows; ++i) {
    for (int k = 0; k < A.cols; ++k) {
      for (int j = 0; j < B.cols; ++j) {
        C(i, j) += A(i, k) * B(k, j);
        stats.multiplications++;
        if (k > 0)
          stats.additions++;
      }
    }
  }
  return C;
}

template <typename T>
Matrix<T> strassen_multiply(const Matrix<T> &A, const Matrix<T> &B) {
  int n = A.rows;
  if (n <= RECURSION_BASE_SIZE) {
    return naive_multiply(A, B);
  }

  Matrix<T> A11(n / 2, n / 2), A12(n / 2, n / 2), A21(n / 2, n / 2),
      A22(n / 2, n / 2);
  Matrix<T> B11(n / 2, n / 2), B12(n / 2, n / 2), B21(n / 2, n / 2),
      B22(n / 2, n / 2);

  A.split(A11, A12, A21, A22);
  B.split(B11, B12, B21, B22);

  Matrix<T> M1 = strassen_multiply(A11 + A22, B11 + B22);
  Matrix<T> M2 = strassen_multiply(A21 + A22, B11);
  Matrix<T> M3 = strassen_multiply(A11, B12 - B22);
  Matrix<T> M4 = strassen_multiply(A22, B21 - B11);
  Matrix<T> M5 = strassen_multiply(A11 + A12, B22);
  Matrix<T> M6 = strassen_multiply(A21 - A11, B11 + B12);
  Matrix<T> M7 = strassen_multiply(A12 - A22, B21 + B22);

  Matrix<T> C11 = M1 + M4 - M5 + M7;
  Matrix<T> C12 = M3 + M5;
  Matrix<T> C21 = M2 + M4;
  Matrix<T> C22 = M1 - M2 + M3 + M6;

  Matrix<T> C(n, n);
  C.join(C11, C12, C21, C22);
  return C;
}

template <typename T>
Matrix<T> winograd_multiply(const Matrix<T> &A, const Matrix<T> &B) {
  if (A.cols != B.rows) {
    throw std::invalid_argument("Matrix dimensions mismatch");
  }
  int n = A.rows;
  int d = n / 2;
  std::vector<T> row_factors(n);
  std::vector<T> col_factors(n);

  for (int i = 0; i < n; ++i) {
    row_factors[i] = T(0);
    for (int j = 0; j < d; ++j) {
      row_factors[i] += A(i, 2 * j) * A(i, 2 * j + 1);
      stats.multiplications++;
      stats.additions++;
    }
  }

  for (int i = 0; i < n; ++i) {
    col_factors[i] = T(0);
    for (int j = 0; j < d; ++j) {
      col_factors[i] += B(2 * j, i) * B(2 * j + 1, i);
      stats.multiplications++;
      stats.additions++;
    }
  }

  Matrix<T> C(n, n);
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j) {
      T sum = -row_factors[i] - col_factors[j];
      stats.additions += 2;

      for (int k = 0; k < d; ++k) {
        sum +=
            (A(i, 2 * k) + B(2 * k + 1, j)) * (A(i, 2 * k + 1) + B(2 * k, j));
        stats.multiplications++;
        stats.additions += 3;
      }
      C(i, j) = sum;
    }
  }

  if (n % 2 != 0) {
    for (int i = 0; i < n; ++i) {
      for (int j = 0; j < n; ++j) {
        C(i, j) += A(i, n - 1) * B(n - 1, j);
        stats.multiplications++;
        stats.additions++;
      }
    }
  }
  return C;
}

// Base template: Naive fallback
template <typename T>
Matrix<T> blas_multiply(const Matrix<T> &A, const Matrix<T> &B) {
  return naive_multiply(A, B);
}

// Specialization for int
template <>
Matrix<int> blas_multiply(const Matrix<int> &A, const Matrix<int> &B) {
#ifdef HAS_ACCELERATE
  if (A.cols != B.rows)
    throw std::invalid_argument("Mismatch");
  int m = A.rows;
  int k = A.cols;
  int n = B.cols;
  std::vector<double> Ad(m * k), Bd(k * n), Cd(m * n);
  for (int i = 0; i < m * k; ++i)
    Ad[i] = (double)A.data[i];
  for (int i = 0; i < k * n; ++i)
    Bd[i] = (double)B.data[i];

  cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, m, n, k, 1.0,
              Ad.data(), k, Bd.data(), n, 0.0, Cd.data(), n);

  Matrix<int> C(m, n);
  for (int i = 0; i < m * n; ++i)
    C.data[i] = (int)Cd[i];
  stats.multiplications += (long long)m * n * k;
  stats.additions += (long long)m * n * k;
  return C;
#else
  return naive_multiply(A, B);
#endif
}

// Specialization for double
template <>
Matrix<double> blas_multiply(const Matrix<double> &A, const Matrix<double> &B) {
#ifdef HAS_ACCELERATE
  if (A.cols != B.rows)
    throw std::invalid_argument("Mismatch");
  int m = A.rows;
  int k = A.cols;
  int n = B.cols;

  Matrix<double> C(m, n);
  cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, m, n, k, 1.0,
              A.data.data(), k, B.data.data(), n, 0.0, C.data.data(), n);

  stats.multiplications += (long long)m * n * k;
  stats.additions += (long long)m * n * k;
  return C;
#else
  return naive_multiply(A, B);
#endif
}

// Specialization for complex<double>
template <>
Matrix<std::complex<double>>
blas_multiply(const Matrix<std::complex<double>> &A,
              const Matrix<std::complex<double>> &B) {
#ifdef HAS_ACCELERATE
  if (A.cols != B.rows)
    throw std::invalid_argument("Mismatch");
  int m = A.rows;
  int k = A.cols;
  int n = B.cols;
  Matrix<std::complex<double>> C(m, n);
  std::complex<double> alpha(1.0, 0.0), beta(0.0, 0.0);

  cblas_zgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, m, n, k, &alpha,
              A.data.data(), k, B.data.data(), n, &beta, C.data.data(), n);

  stats.multiplications += (long long)m * n * k;
  stats.additions += (long long)m * n * k;
  return C;
#else
  return naive_multiply(A, B);
#endif
}

// Explicit Instantiation
template Matrix<int> naive_multiply(const Matrix<int> &, const Matrix<int> &);
template Matrix<int> strassen_multiply(const Matrix<int> &,
                                       const Matrix<int> &);
template Matrix<int> winograd_multiply(const Matrix<int> &,
                                       const Matrix<int> &);

template Matrix<double> naive_multiply(const Matrix<double> &,
                                       const Matrix<double> &);
template Matrix<double> strassen_multiply(const Matrix<double> &,
                                          const Matrix<double> &);
template Matrix<double> winograd_multiply(const Matrix<double> &,
                                          const Matrix<double> &);

template Matrix<std::complex<double>>
naive_multiply(const Matrix<std::complex<double>> &,
               const Matrix<std::complex<double>> &);
template Matrix<std::complex<double>>
strassen_multiply(const Matrix<std::complex<double>> &,
                  const Matrix<std::complex<double>> &);
template Matrix<std::complex<double>>
winograd_multiply(const Matrix<std::complex<double>> &,
                  const Matrix<std::complex<double>> &);
