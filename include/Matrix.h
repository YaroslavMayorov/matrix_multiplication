#ifndef MATRIX_H
#define MATRIX_H

#include "Operations.h"
#include <iostream>
#include <stdexcept>
#include <vector>

template <typename T> class Matrix {
public:
  int rows;
  int cols;
  std::vector<T> data;

  Matrix(int r, int c) : rows(r), cols(c), data(r * c) {}

  // Accessor
  T &operator()(int r, int c) { return data[r * cols + c]; }
  const T &operator()(int r, int c) const { return data[r * cols + c]; }

  // Static helper to create a matrix with random values
  static Matrix random(int r, int c);
  static Matrix symmetric(int n);
  static Matrix identity(int n);

  // Matrix addition
  Matrix operator+(const Matrix &other) const;
  Matrix operator-(const Matrix &other) const;

  // Split matrix into 4 quadrants
  void split(Matrix &A11, Matrix &A12, Matrix &A21, Matrix &A22) const;

  // Join 4 quadrants into this matrix
  void join(const Matrix &A11, const Matrix &A12, const Matrix &A21,
            const Matrix &A22);

  void print() const;
};

#endif // MATRIX_H
