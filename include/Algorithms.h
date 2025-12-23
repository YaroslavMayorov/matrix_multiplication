#ifndef ALGORITHMS_H
#define ALGORITHMS_H

#include "Matrix.h"

template <typename T>
Matrix<T> naive_multiply(const Matrix<T> &A, const Matrix<T> &B);

template <typename T>
Matrix<T> strassen_multiply(const Matrix<T> &A, const Matrix<T> &B);

template <typename T>
Matrix<T> winograd_multiply(const Matrix<T> &A, const Matrix<T> &B);

template <typename T>
Matrix<T> blas_multiply(const Matrix<T> &A, const Matrix<T> &B);

#endif // ALGORITHMS_H
