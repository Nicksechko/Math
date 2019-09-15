#pragma once

#include <iostream>
#include <fstream>
#include <stdexcept>
#include <utility>
#include <vector>
#include <numeric>
#include <cassert>

#include "fraction.h"

enum class ChoiceType {
  Without,
  Row,
  Column,
  Submatrix
};

struct LUPMatrix;

class LinearSystem;

class Matrix {
 private:
  int n_{};
  int m_{};

  std::vector<std::vector<Fraction>> elements_;

 public:
  explicit Matrix(int n, int m)
      : n_(n), m_(m),
        elements_(std::vector<std::vector<Fraction>>(n, std::vector<Fraction>(m))) {}

  Matrix() = default;

//  Matrix(Matrix&& matrix) noexcept = default;
//
//  Matrix(const Matrix& matrix) = default;
//
//  Matrix& operator=(Matrix&& matrix) noexcept = default;
//
//  Matrix& operator=(const Matrix& matrix) = default;
//
  explicit Matrix(std::vector<std::vector<Fraction>> elements)
      : n_(elements.size()), m_((n_ == 0) ? 0 : static_cast<int>(elements[0].size())),
        elements_(std::move(elements)) {}

  [[nodiscard]] Matrix ToTranspose() const {
    Matrix result(m_, n_);
    for (int i = 0; i < m_; ++i) {
      for (int j = 0; j < n_; ++j) {
        result.At(i, j) = At(j, i);
      }
    }

    return result;
  }

  Fraction& At(int i, int j) {
    return elements_.at(i).at(j);
  }

  [[nodiscard]] const Fraction& At(int i, int j) const {
    return elements_.at(i).at(j);
  }

  [[nodiscard]] int GetNumRows() const {
    return n_;
  }

  [[nodiscard]] int GetNumColumns() const {
    return m_;
  }

  void MultiplyRow(int row, Fraction coefficient);

  void DivideRow(int row, Fraction coefficient);

  void MultiplyColumn(int column, Fraction coefficient);

  void DivideColumn(int column, Fraction coefficient);

  void AddRows(int first_row, int second_row, Fraction coefficient = Fraction(1));

  void SubtractRows(int minuend_row, int subtrahend_row, Fraction coefficient = Fraction(1));

  void SwapRows(int first_row, int second_row);

  void AddColumns(int first_column, int second_column, Fraction coefficient = Fraction(1));

  void SubtractColumns(int minuend_column, int subtrahend_column, Fraction coefficient = Fraction(1));

  void SwapColumns(int first_column, int second_column);

  void CopyRow(int source_row, int destination_row);

  void CopyColumn(int source_column, int destination_column);

  void CopyRow(const Matrix& source_matrix, int source_row, int destination_row);

  void CopyColumn(const Matrix& source_matrix, int source_column, int destination_column);

  bool Inverse();

  bool LUPTransform(LUPMatrix& lup_matrix);

  static Matrix Identity(int n);
};

struct LUPMatrix {
  LUPMatrix() = default;

  LUPMatrix(Matrix L, Matrix U, std::vector<int> row_permutation,
            std::vector<int> column_permutation);

  LUPMatrix& operator=(LUPMatrix&& lup_matrix) noexcept = default;

  LUPMatrix& operator=(const LUPMatrix& lup_matrix) = default;

  Matrix L, U;
  std::vector<int> row_permutation, column_permutation;
};

class LinearSystem {
 public:
  LinearSystem(Matrix A, Matrix B);

  int RunDirectGauss(ChoiceType choice_type = ChoiceType::Submatrix);

  void RunReverseGauss();

  int RunGauss(ChoiceType choice_type = ChoiceType::Submatrix);

  Matrix GetSolutionMatrix();

  LUPMatrix ToLUPMatrix();

 private:
  void Choice(int pos, ChoiceType choice_type = ChoiceType::Submatrix);

  int n_, rank_;
  Matrix A_, B_;
  std::vector<int> row_permutation_, column_permutation_;
};


bool operator==(const Matrix& one, const Matrix& two) {
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

Matrix operator+(const Matrix& one, const Matrix& two) {
  if (one.GetNumRows() != two.GetNumRows()) {
    throw std::invalid_argument("[SUM] Mismatched number of rows");
  }

  if (one.GetNumColumns() != two.GetNumColumns()) {
    throw std::invalid_argument("[SUM] Mismatched number of columns");
  }

  Matrix result(one.GetNumRows(), one.GetNumColumns());
  for (int i = 0; i < result.GetNumRows(); ++i) {
    for (int j = 0; j < result.GetNumColumns(); ++j) {
      result.At(i, j) = one.At(i, j) + two.At(i, j);
    }
  }

  return result;
}


Matrix operator*(const Matrix& one, const Matrix& two) {
  if (one.GetNumColumns() != two.GetNumRows()) {
    throw std::invalid_argument("[MUL] Mismatched number of rows");
  }

  Matrix result(one.GetNumRows(), two.GetNumColumns());
  for (int i = 0; i < result.GetNumRows(); ++i) {
    for (int j = 0; j < result.GetNumColumns(); ++j) {
      for (int k = 0; k < one.GetNumColumns(); ++k) {
        result.At(i, j) += one.At(i, k) * two.At(k, j);
      }
    }
  }

  return result;
}


std::istream& operator>>(std::istream& in, Matrix& matrix) {
  int n, m;
  in >> n >> m;

  matrix = Matrix(n, m);
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < m; ++j) {
      in >> matrix.At(i, j);
    }
  }

  return in;
}

std::ostream& operator<<(std::ostream& out, const Matrix& matrix) {
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

void Matrix::AddRows(int first_row, int second_row, Fraction coefficient) {
  for (int column = 0; column < m_; ++column) {
    elements_[first_row][column] += coefficient * elements_[second_row][column];
  }
}

void Matrix::SubtractRows(int minuend_row, int subtrahend_row, Fraction coefficient) {
  AddRows(minuend_row, subtrahend_row, -coefficient);
}

void Matrix::SwapRows(int first_row, int second_row) {
  for (int column = 0; column < m_; ++column) {
    std::swap(elements_[first_row][column], elements_[second_row][column]);
  }
}


void Matrix::AddColumns(int first_column, int second_column, Fraction coefficient) {
  for (int row = 0; row < n_; ++row) {
    elements_[row][first_column] += coefficient * elements_[row][second_column];
  }
}

void Matrix::SubtractColumns(int minuend_column, int subtrahend_column, Fraction coefficient) {
  AddColumns(minuend_column, subtrahend_column, -coefficient);
}

void Matrix::SwapColumns(int first_column, int second_column) {
  for (int row = 0; row < n_; ++row) {
    std::swap(elements_[row][first_column], elements_[row][second_column]);
  }
}

void Matrix::MultiplyRow(int row, Fraction coefficient) {
  for (int column = 0; column < m_; ++column) {
    elements_[row][column] *= coefficient;
  }
}

void Matrix::DivideRow(int row, Fraction coefficient) {
  MultiplyRow(row, 1 / coefficient);
}

void Matrix::MultiplyColumn(int column, Fraction coefficient) {
  for (int row = 0; row < n_; ++row) {
    elements_[row][column] *= coefficient;
  }
}

void Matrix::DivideColumn(int column, Fraction coefficient) {
  assert(coefficient == Fraction(0));
  MultiplyColumn(column, 1 / coefficient);
}

void Matrix::CopyRow(int source_row, int destination_row) {
  for (int column = 0; column < m_; ++column) {
    elements_[destination_row][column] = elements_[source_row][column];
  }
}

void Matrix::CopyColumn(int source_column, int destination_column) {
  for (int row = 0; row < n_; ++row) {
    elements_[row][destination_column] = elements_[row][source_column];
  }
}

void Matrix::CopyRow(const Matrix& source_matrix, int source_row, int destination_row) {
  for (int column = 0; column < m_; ++column) {
    elements_[destination_row][column] = source_matrix.At(source_row, column);
  }
}

void Matrix::CopyColumn(const Matrix& source_matrix, int source_column, int destination_column) {
  for (int row = 0; row < n_; ++row) {
    elements_[row][destination_column] = source_matrix.At(row, source_column);
  }
}

bool Matrix::Inverse() {
  assert(n_ == m_);
  LinearSystem inverser(*this, Identity(n_));
  operator<<(std::cout, *this);
  int rank = inverser.RunGauss();
  if (rank != n_) {
    return false;
  }
  operator=(inverser.GetSolutionMatrix());

  return true;
}

bool Matrix::LUPTransform(LUPMatrix& lup_matrix) {
  assert(n_ == m_);
  LinearSystem transformer(*this, Identity(n_));
  int rank = transformer.RunDirectGauss(ChoiceType::Row);
  if (rank < n_) {
    return false;
  }
  lup_matrix = transformer.ToLUPMatrix();
  assert(lup_matrix.L.Inverse());

  return true;
}

Matrix Matrix::Identity(int n) {
  Matrix identity(n, n);
  for (int i = 0; i < n; i++) {
    identity.At(i, i) = Fraction(1);
  }

  return identity;
}

LUPMatrix::LUPMatrix(Matrix L,
                        Matrix U,
                        std::vector<int> row_permutation,
                        std::vector<int> column_permutation)
    : L(std::move(L)),
      U(std::move(U)),
      row_permutation(std::move(row_permutation)),
      column_permutation(std::move(column_permutation)) {}

LinearSystem::LinearSystem(Matrix A, Matrix B)
    : n_(A.GetNumRows()),
      rank_(-1),
      A_(std::move(A)),
      B_(std::move(B)),
      row_permutation_(A.GetNumRows()),
      column_permutation_(A.GetNumColumns()) {
  assert(A_.GetNumRows() == B_.GetNumRows());
  assert(A_.GetNumRows() == A_.GetNumColumns());
  std::iota(row_permutation_.begin(), row_permutation_.end(), 0);
  std::iota(column_permutation_.begin(), column_permutation_.end(), 0);
}

int LinearSystem::RunDirectGauss(ChoiceType choice_type) {
  for (int pos = 0; pos < n_; ++pos) {
    Choice(pos, choice_type);
    if (A_.At(pos, pos) == Fraction(0)) {
      Choice(pos);
    }
    if (A_.At(pos, pos) == Fraction(0)) {
      return rank_ = pos;
    }
    for (int row = pos + 1; row < A_.GetNumRows(); ++row) {
      Fraction coefficient = A_.At(row, pos) / A_.At(pos, pos);
      A_.SubtractRows(row, pos, coefficient);
      B_.SubtractRows(row, pos, coefficient);
    }
  }

  return rank_ = n_;
}

void LinearSystem::RunReverseGauss() {
  for (int pos = n_ - 1; pos >= 0; --pos) {
    if (A_.At(pos, pos) == Fraction(0)) {
      continue;
    }
    B_.DivideRow(pos, A_.At(pos, pos));
    A_.DivideRow(pos, A_.At(pos, pos));
    for (int row = pos - 1; row >= 0; --row) {
      Fraction coefficient = A_.At(row, pos);
      A_.SubtractRows(row, pos, coefficient);
      B_.SubtractRows(row, pos, coefficient);
    }
  }
}

int LinearSystem::RunGauss(ChoiceType choice_type) {
  RunDirectGauss(choice_type);
  RunReverseGauss();

  return rank_;
}

Matrix LinearSystem::GetSolutionMatrix() {
  Matrix solution_matrix(B_.GetNumRows(), B_.GetNumColumns());
  for (int row = 0; row < static_cast<int>(column_permutation_.size()); row++) {
    solution_matrix.CopyRow(B_, row, column_permutation_[row]);
  }

  return solution_matrix;
}

LUPMatrix LinearSystem::ToLUPMatrix() {
  return LUPMatrix(B_, A_, row_permutation_, column_permutation_);
}

void LinearSystem::Choice(int pos, ChoiceType choice_type) {
  int max_row = n_, max_column = n_;
  if (choice_type == ChoiceType::Without) {
    return;
  } else if (choice_type == ChoiceType::Row) {
    max_row = pos + 1;
  } else if (choice_type == ChoiceType::Column) {
    max_column = pos + 1;
  }

  int best_column = pos, best_row = pos;
  Fraction best_item = A_.At(pos, pos);
  for (int row = pos; row < max_row; ++row) {
    for (int column = pos; column < max_column; ++column) {
      if (best_item < fabs(A_.At(row, column))) {
        best_row = row;
        best_column = column;
        best_item = fabs(A_.At(row, column));
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
