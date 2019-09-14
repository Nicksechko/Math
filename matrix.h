#pragma once

#include <iostream>
#include <fstream>
#include <stdexcept>
#include <utility>
#include <vector>
#include <numeric>
#include <cassert>

enum class ChoiceType {
  Without,
  Row,
  Column,
  Submatrix
};

template<class T>
struct LUPMatrix;

template<class T>
class LinearSystem;

template<class T>
class Matrix {
 private:
  int n_{};
  int m_{};

  std::vector<std::vector<T>> elements_;

 public:
  explicit Matrix(int n, int m)
      : n_(n), m_(m),
        elements_(std::vector<std::vector<T>>(n, std::vector<T>(m))) {}

  Matrix() = default;

  Matrix(Matrix&& matrix) noexcept
    : n_(std::move(matrix.n_)),
      m_(std::move(matrix.m_)),
      elements_(std::move(matrix.elements_)) {}

  Matrix(const Matrix& matrix)
      : n_(matrix.n_),
        m_(matrix.m_),
        elements_(matrix.elements_) {}

  Matrix& operator=(Matrix&& matrix) noexcept {
    n_ = std::move(matrix.n_);
    m_ = std::move(matrix.m_);
    elements_ = std::move(matrix.elements_);
  }

  Matrix& operator=(const Matrix& matrix) {
    n_ = matrix.n_;
    m_ = matrix.m_;
    elements_ = matrix.elements_;
  }

  explicit Matrix(std::vector<std::vector<T>> elements)
      : n_(elements.size()), m_((n_ == 0) ? 0 : elements[0].size()), elements_(std::move(elements)) {}

  Matrix ToTranspose() const {
    Matrix result(m_, n_);
    for (int i = 0; i < m_; ++i) {
      for (int j = 0; j < n_; ++j) {
        result.At(i, j) = this->At(j, i);
      }
    }

    return result;
  }

  T& At(int i, int j) {
    return elements_.at(i).at(j);
  }

  const T& At(int i, int j) const {
    return elements_.at(i).at(j);
  }

  [[nodiscard]] int GetNumRows() const {
    return n_;
  }

  [[nodiscard]] int GetNumColumns() const {
    return m_;
  }

  void MultiplyRow(int row, T coefficient);

  void DivideRow(int row, T coefficient);

  void MultiplyColumn(int column, T coefficient);

  void DivideColumn(int column, T coefficient);

  void AddRows(int first_row, int second_row, T coefficient = T(1));

  void SubtractRows(int minuend_row, int subtrahend_row, T coefficient = T(1));

  void SwapRows(int first_row, int second_row);

  void AddColumns(int first_column, int second_column, T coefficient = T(1));

  void SubtractColumns(int minuend_column, int subtrahend_column, T coefficient = T(1));

  void SwapColumns(int first_column, int second_column);

  bool Inverse();

  bool LUPTransform(LUPMatrix<T>& lup_matrix);

  static Matrix<T> Identity(int n);
};

template<class T>
struct LUPMatrix {
  LUPMatrix() = default;

  LUPMatrix(Matrix<T> L, Matrix<T> U, std::vector<int> row_permutation,
            std::vector<int> column_permutation);

  LUPMatrix<T>& operator=(LUPMatrix<T>&& lup_matrix) noexcept {
    L = std::move(lup_matrix.L);
    U = std::move(lup_matrix.U);
    row_permutation = std::move(lup_matrix.row_permutation);
    column_permutation = std::move(lup_matrix.column_permutation);
  }

  LUPMatrix<T>& operator=(const LUPMatrix<T>& lup_matrix) {
    L = lup_matrix.L;
    U = lup_matrix.U;
    row_permutation = lup_matrix.row_permutation;
    column_permutation = lup_matrix.column_permutation;
  }

  Matrix<T> L, U;
  std::vector<int> row_permutation, column_permutation;
};

template<class T>
class LinearSystem {
 public:
  explicit LinearSystem(Matrix<T> A);

  LinearSystem(Matrix<T> A, Matrix<T> B);

  int RunDirectGauss(ChoiceType choice_type = ChoiceType::Submatrix);

  void RunReverseGauss();

  int RunGauss(ChoiceType choice_type = ChoiceType::Submatrix);

  Matrix<T> GetSolutionMatrix();

  LUPMatrix<T> ToLUPMatrix();

 private:
  void Choice(int pos, ChoiceType choice_type = ChoiceType::Submatrix);

  int rank_;
  Matrix<T> A_, B_;
  std::vector<int> row_permutation_, column_permutation_;
};

template<class T>
bool operator==(const Matrix<T>& one, const Matrix<T>& two) {
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
Matrix<T> operator+(const Matrix<T>& one, const Matrix<T>& two) {
  if (one.GetNumRows() != two.GetNumRows()) {
    throw std::invalid_argument("[SUM] Mismatched number of rows");
  }

  if (one.GetNumColumns() != two.GetNumColumns()) {
    throw std::invalid_argument("[SUM] Mismatched number of columns");
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
Matrix<T> operator*(const Matrix<T>& one, const Matrix<T>& two) {
  if (one.GetNumColumns() != two.GetNumRows()) {
    throw std::invalid_argument("[MUL] Mismatched number of rows");
  }

  Matrix<T> result(one.GetNumRows(), two.GetNumColumns());
  for (int i = 0; i < result.GetNumRows(); ++i) {
    for (int j = 0; j < result.GetNumColumns(); ++j) {
      for (int k = 0; k < one.GetNumColumns(); ++k) {
        result.At(i, j) += one.At(i, k) * two.At(k, j);
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
std::ostream& operator<<(std::ostream& out, const Matrix<T>& matrix) {
  for (int i = 0; i < matrix.GetNumRows(); ++i) {
    for (int j = 0; j < matrix.GetNumColumns(); ++j) {
      if (j > 0) {
        out << " ";
      }
      out << matrix.At(i, j);
    }
    out << std::endl;
  }

  return out;
}

template<class T>
void Matrix<T>::AddRows(int first_row, int second_row, T coefficient) {
  for (int column = 0; column < m_; ++column) {
    elements_[first_row][column] += coefficient * elements_[second_row][column];
  }
}

template<class T>
void Matrix<T>::SubtractRows(int minuend_row, int subtrahend_row, T coefficient) {
  AddRows(minuend_row, subtrahend_row, -coefficient);
}

template<class T>
void Matrix<T>::SwapRows(int first_row, int second_row) {
  for (int column = 0; column < m_; ++column) {
    std::swap(elements_[first_row][column], elements_[second_row][column]);
  }
}

template<class T>
void Matrix<T>::AddColumns(int first_column, int second_column, T coefficient) {
  for (int row = 0; row < n_; ++row) {
    elements_[row][first_column] += coefficient * elements_[row][second_column];
  }
}

template<class T>
void Matrix<T>::SubtractColumns(int minuend_column, int subtrahend_column, T coefficient) {
  AddColumns(minuend_column, subtrahend_column, -coefficient);
}

template<class T>
void Matrix<T>::SwapColumns(int first_column, int second_column) {
  for (int row = 0; row < n_; ++row) {
    std::swap(elements_[row][first_column], elements_[row][second_column]);
  }
}
template<class T>
void Matrix<T>::MultiplyRow(int row, T coefficient) {
  for (int column = 0; column < m_; ++column) {
    elements_[row][column] *= coefficient;
  }
}

template<class T>
void Matrix<T>::DivideRow(int row, T coefficient) {
  MultiplyRow(row, 1 / coefficient);
}

template<class T>
void Matrix<T>::MultiplyColumn(int column, T coefficient) {
  for (int row = 0; row < n_; ++row) {
    elements_[row][column] *= coefficient;
  }
}

template<class T>
void Matrix<T>::DivideColumn(int column, T coefficient) {
  assert(coefficient == T(0));
  MultiplyColumn(column, 1 / coefficient);
}

template<class T>
bool Matrix<T>::Inverse() {
  assert(n_ == m_);
  LinearSystem<T> inverser(*this, Identity(n_));
  int rank = inverser.RunGauss();
  if (rank != n_) {
    return false;
  }
  operator=(inverser.GetSolutionMatrix());
  return true;
}

template<class T>
bool Matrix<T>::LUPTransform(LUPMatrix<T>& lup_matrix) {
  assert(n_ == m_);
  LinearSystem<T> transformer(*this, Identity(n_));
  int rank = transformer.RunDirectGauss(ChoiceType::Row);
  if (rank < n_) {
    return false;
  }
  lup_matrix = transformer.ToLUPMatrix();
  assert(lup_matrix.L.Inverse());
  return true;
}

template<class T>
Matrix<T> Matrix<T>::Identity(int n) {
  Matrix<T> identity(n, n);
  for (int i = 0; i < n; i++) {
    identity.At(i, i) = T(1);
  }
}

template<class T>
LUPMatrix<T>::LUPMatrix(Matrix<T> L,
                        Matrix<T> U,
                        std::vector<int> row_permutation,
                        std::vector<int> column_permutation)
    : L(L),
      U(U),
      row_permutation(std::move(row_permutation)),
      column_permutation(std::move(column_permutation)) {}

template<class T>
LinearSystem<T>::LinearSystem(Matrix<T> A) {
  LinearSystem(A, Matrix(A.GetNumRows(), 1));
}

template<class T>
LinearSystem<T>::LinearSystem(Matrix<T> A, Matrix<T> B)
    : rank_(-1),
      A_(A),
      B_(B),
      row_permutation_(A.GetNumRows()),
      column_permutation_(A.GetNumColumns()) {
  assert(A_.GetNumRows() == B_.GetNumRows());
  std::iota(row_permutation_.begin(), row_permutation_.end(), 0);
  std::iota(column_permutation_.begin(), column_permutation_.end(), 0);
}

template<class T>
int LinearSystem<T>::RunDirectGauss(ChoiceType choice_type) {
  int max_pos = std::min(A_.GetNumRows(), A_.GetNumColumns());
  for (int pos = 0; pos < max_pos; ++pos) {
    Choice(pos, choice_type);
    if (A_.At(pos, pos) == T(0)) {
      Choice(pos);
    }
    if (A_.At(pos, pos) == T(0)) {
      return rank_ = pos;
    }
    for (int row = pos + 1; row < A_.GetNumRows(); ++row) {
      T coefficient = A_.At(row, pos) / A_.At(pos, pos);
      A_.SubtractRows(row, pos, coefficient);
      B_.SubtractRows(row, pos, coefficient);
    }
  }
  return rank_ = max_pos;
}

template<class T>
void LinearSystem<T>::RunReverseGauss() {
  int max_pos = std::min(A_.GetNumRows(), A_.GetNumColumns());
  for (int pos = max_pos - 1; pos >= 0; --pos) {
    if (A_.At(pos, pos) == T(0)) {
      continue;
    }
    A_.DivideRow(pos, A_.At(pos, pos));
    for (int row = pos - 1; row >= 0; --row) {
      T coefficient = A_.At(row, pos);
      A_.SubtractRows(row, pos, coefficient);
      B_.SubtractRows(row, pos, coefficient);
    }
  }
}

template<class T>
int LinearSystem<T>::RunGauss(ChoiceType choice_type) {
  RunDirectGauss(choice_type);
  RunReverseGauss();
  return rank_;
}

template<class T>
Matrix<T> LinearSystem<T>::GetSolutionMatrix() {
  return B_;
}

template<class T>
LUPMatrix<T> LinearSystem<T>::ToLUPMatrix() {
  return LUPMatrix<T>(B_, A_, row_permutation_, column_permutation_);
}

template<class T>
void LinearSystem<T>::Choice(int pos, ChoiceType choice_type) {
  int max_row = A_.GetNumColumns(), max_column = A_.GetNumColumns();
  if (choice_type == ChoiceType::Without) {
    return;
  } else if (choice_type == ChoiceType::Row) {
    max_row = pos + 1;
  } else if (choice_type == ChoiceType::Column) {
    max_column = pos + 1;
  }

  int best_column = pos, best_row = pos;
  T best_item = A_.At(pos, pos);
  for (int row = pos; row < max_row; ++row) {
    for (int column = pos; column < max_column; ++column) {
      if (best_item < fabs(A_.At(row, column))) {
        best_row = row;
        best_column = column;
        best_item = A_.At(row, column);
      }
    }
  }

  if (pos != best_row) {
    std::swap(row_permutation_[pos], row_permutation_[best_row]);
    A_.SwapRows(pos, best_row);
    B_.SwapRows(pos, best_row);
  }
  if (pos != best_column) {
    std::swap(column_permutation_[pos], column_permutation_[best_column]);
    A_.SwapColumns(pos, best_column);
  }
}
