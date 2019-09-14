#pragma once

#include <iostream>
#include <fstream>
#include <stdexcept>
#include <utility>
#include <vector>

template<class T>
class Matrix {
 private:
  std::vector<std::vector<T>> elements_;

  const int n_{};
  const int m_{};

 public:
  explicit Matrix(int n, int m)
      : elements_(std::vector<std::vector<T>>(n, std::vector<T>(m)))
      , n_(n)
      , m_(m){}

  Matrix() = default;

  explicit Matrix(std::vector<std::vector<T>> elements)
      : elements_(std::move(elements))
      , n_(elements.size())
      , m_((n_ == 0) ? 0 : elements[0].size()) {}

  T &At(int i, int j) {
    return elements_.at(i).at(j);
  }

  const T &At(int i, int j) const {
    return elements_.at(i).at(j);
  }

  [[nodiscard]] int GetNumRows() const {
    return n_;
  }

  [[nodiscard]] int GetNumColumns() const {
    return m_;
  }
};

template<class T>
bool operator==(const Matrix<T> &one, const Matrix<T> &two) {
  if (one.GetNumRows() != two.GetNumRows() ||
      one.GetNumColumns() != two.GetNumColumns()) {
    return false;
  }

  for (int i = 0; i < one.GetNumRows(); ++i) {
    for (int j = 0; j < one.GetNumColumns(); ++j) {
      if (one.At(i, j) != two.At(i, j)) {
        return false;
      }
    }
  }

  return true;
}

template<class T>
Matrix<T> operator+(const Matrix<T> &one, const Matrix<T> &two) {
  if (one.GetNumRows() != two.GetNumRows()) {
    throw std::invalid_argument("Mismatched number of rows");
  }

  if (one.GetNumColumns() != two.GetNumColumns()) {
    throw std::invalid_argument("Mismatched number of columns");
  }

  Matrix<T> result(one.GetNumRows(), one.GetNumColumns());
  for (int i = 0; i < result.GetNumRows(); ++i) {
    for (int j = 0; j < result.GetNumColumns(); ++j) {
      result.At(i, j) = one.At(i, j) + two.At(i, j);
    }
  }

  return result;
}

template<class T>
Matrix<T> operator*(const Matrix<T> &one, const Matrix<T> &two) {
  if (one.GetNumColumns() != two.GetNumRows()) {
    throw std::invalid_argument("Mismatched number of rows");
  }

  Matrix<T> result(one.GetNumRows(), two.GetNumColumns());
  for (int i = 0; i < result.GetNumRows(); ++i) {
    for (int j = 0; j < result.GetNumColumns(); ++j) {
      for (int k = 0; k < one.GetNumColumns(); ++k) {
        result.At(i, j) += one.At(i, k) + two.At(k, j);
      }
    }
  }

  return result;
}

template<class T>
std::istream& operator>>(std::istream& in, Matrix<T>& matrix) {
  int n, m;
  in >> n >> m;

  matrix.Reset(n, m);
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < m; ++j) {
      in >> matrix.At(i, j);
    }
  }

  return in;
}

template<class T>
std::ostream &operator<<(std::ostream &out, const Matrix<T> &matrix) {
  for (int i = 0; i < matrix.GetNumRows(); ++i) {
    for (int j = 0; j < matrix.GetNumColumns(); ++j) {
      if (j > 0) {
        out << " ";
      }
      out << matrix.At(row, column);
    }
    out << std::endl;
  }

  return out;
}