#include "../include/Matrix.h"
#include "../include/Operations.h"
#include <ctime>
#include <random>

// Global stats instance
OperationStats stats;

#include <complex>

template <typename T> Matrix<T> Matrix<T>::random(int r, int c) {
  Matrix<T> m(r, c);
  static std::mt19937 gen(time(0));

  if constexpr (std::is_same_v<T, int>) {
    std::uniform_int_distribution<> dis(0, 9);
    for (int i = 0; i < r * c; ++i) {
      m.data[i] = dis(gen);
    }
  } else if constexpr (std::is_same_v<T, std::complex<double>>) {
    std::uniform_real_distribution<> dis(-10.0, 10.0);
    for (int i = 0; i < r * c; ++i) {
      m.data[i] = std::complex<double>(dis(gen), dis(gen));
    }
  }
  return m;
}

template <typename T> Matrix<T> Matrix<T>::symmetric(int n) {
  Matrix<T> m(n, n);
  static std::mt19937 gen(time(0));

  if constexpr (std::is_same_v<T, int>) {
    std::uniform_int_distribution<> dis(0, 9);
    for (int i = 0; i < n; ++i) {
      for (int j = i; j < n; ++j) {
        int val = dis(gen);
        m(i, j) = val;
        m(j, i) = val;
      }
    }
  } else if constexpr (std::is_same_v<T, std::complex<double>>) {
    std::uniform_real_distribution<> dis(-10.0, 10.0);
    for (int i = 0; i < n; ++i) {
      for (int j = i; j < n; ++j) {
        std::complex<double> val(dis(gen), dis(gen));
        m(i, j) = val;
        m(j, i) = val;
      }
    }
  }
  return m;
}

template <typename T> Matrix<T> Matrix<T>::identity(int n) {
  Matrix<T> m(n, n);
  // Default init to 0
  for (int i = 0; i < n * n; ++i)
    m.data[i] = T(0);
  for (int i = 0; i < n; ++i) {
    m(i, i) = T(1);
  }
  return m;
}

template <typename T>
Matrix<T> Matrix<T>::operator+(const Matrix<T> &other) const {
  Matrix<T> res(rows, cols);
  for (int i = 0; i < rows * cols; ++i) {
    res.data[i] = data[i] + other.data[i];
    stats.additions++; // Count additions
  }
  return res;
}

template <typename T>
Matrix<T> Matrix<T>::operator-(const Matrix<T> &other) const {
  Matrix<T> res(rows, cols);
  for (int i = 0; i < rows * cols; ++i) {
    res.data[i] = data[i] - other.data[i];
    stats.additions++; // Count subtractions as additions
  }
  return res;
}

template <typename T>
void Matrix<T>::split(Matrix<T> &A11, Matrix<T> &A12, Matrix<T> &A21,
                      Matrix<T> &A22) const {
  int newSize = rows / 2;
  for (int i = 0; i < newSize; ++i) {
    for (int j = 0; j < newSize; ++j) {
      A11(i, j) = (*this)(i, j);
      A12(i, j) = (*this)(i, j + newSize);
      A21(i, j) = (*this)(i + newSize, j);
      A22(i, j) = (*this)(i + newSize, j + newSize);
    }
  }
}

template <typename T>
void Matrix<T>::join(const Matrix<T> &A11, const Matrix<T> &A12,
                     const Matrix<T> &A21, const Matrix<T> &A22) {
  int newSize = A11.rows;
  for (int i = 0; i < newSize; ++i) {
    for (int j = 0; j < newSize; ++j) {
      (*this)(i, j) = A11(i, j);
      (*this)(i, j + newSize) = A12(i, j);
      (*this)(i + newSize, j) = A21(i, j);
      (*this)(i + newSize, j + newSize) = A22(i, j);
    }
  }
}

template <typename T> void Matrix<T>::print() const {
  for (int i = 0; i < rows; ++i) {
    for (int j = 0; j < cols; ++j) {
      std::cout << (*this)(i, j) << " ";
    }
    std::cout << std::endl;
  }
}

// Explicit instantiation
template class Matrix<int>;
template class Matrix<double>;
template class Matrix<std::complex<double>>;
