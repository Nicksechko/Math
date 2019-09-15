#pragma once

#include <iostream>
#include <fstream>
#include <stdexcept>
#include <utility>
#include <vector>
#include <numeric>
#include <cassert>
#include <sstream>

#include "fraction.h"

namespace Options {
enum class ChoiceType {
  Without,
  Row,
  Column,
  Submatrix
};

enum class StepByStepType {
  Without,
  MainSteps,
  AllSteps
};

enum class OutputType {
  Standard,
  LaTex
};

extern OutputType output_type;
extern StepByStepType step_by_step_type;
extern int latex_block_size;
}

struct LUPMatrix;

class LinearSystem;

class Matrix {
 public:
  explicit Matrix(int n, int m)
      : n_(n), m_(m),
        elements_(std::vector<std::vector<Fraction>>(n, std::vector<Fraction>(m))) {}

  Matrix() = default;

  explicit Matrix(std::vector<std::vector<Fraction>> elements)
      : n_(elements.size()), m_((n_ == 0) ? 0 : static_cast<int>(elements[0].size())),
        elements_(std::move(elements)) {}

  [[nodiscard]] Matrix ToTranspose() const;

  std::string ToLatex();

  Fraction& At(int i, int j);

  [[nodiscard]] const Fraction& At(int i, int j) const;

  [[nodiscard]] int GetNumRows() const;

  [[nodiscard]] int GetNumColumns() const;

  friend bool operator==(const Matrix& lhs, const Matrix& rhs);

  friend Matrix operator+(const Matrix& lhs, const Matrix& rhs);

  friend Matrix operator*(const Matrix& lhs, const Matrix& rhs);

  friend std::istream& operator>>(std::istream& in, Matrix& matrix);

  friend std::ostream& operator<<(std::ostream& out, const Matrix& matrix);

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

 private:
  int n_;
  int m_;

  std::vector<std::vector<Fraction>> elements_;
};

struct LUPMatrix {
  LUPMatrix() = default;

  LUPMatrix(Matrix L, Matrix U, std::vector<int> row_permutation,
            std::vector<int> column_permutation);

  std::string ToLatex();

  Matrix L, U;
  std::vector<int> row_permutation, column_permutation;
};

class LinearSystem {
 public:
  LinearSystem(Matrix A, Matrix B);

  int RunDirectGauss(Options::ChoiceType choice_type = Options::ChoiceType::Submatrix);

  void RunReverseGauss();

  int RunGauss(Options::ChoiceType choice_type = Options::ChoiceType::Submatrix);

  void OutputSystem();

  Matrix GetSolutionMatrix();

  LUPMatrix ToLUPMatrix();

 private:
  void Choice(int pos, Options::ChoiceType choice_type = Options::ChoiceType::Submatrix);

  bool first;
  int block_flag;
  int n_, rank_;
  Matrix A_, B_;
  std::vector<int> row_permutation_, column_permutation_;
};
