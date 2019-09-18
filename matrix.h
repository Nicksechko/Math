#pragma once

#include <iostream>
#include <stdexcept>
#include <vector>
#include <numeric>
#include <cassert>
#include <algorithm>
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

extern std::ostringstream writer;
extern OutputType output_type;
extern StepByStepType step_by_step_type;
extern int max_matrices_count_in_row;
}

class Permutation;

class Matrix;

class ExtendedMatrix;

class LUPMatrix;

class LinearSystem;

class Permutation {
 public:
  Permutation() = default;

  Permutation(const Permutation& permutation) = default;

  explicit Permutation(int n);

  explicit Permutation(std::vector<int> permutation);

  bool Transposition(int i, int j);

  [[nodiscard]] int At(int pos) const;

  [[nodiscard]] int Size() const;

  [[nodiscard]] Matrix AsMatrix() const;

  [[nodiscard]] Permutation GetInverse() const;

  void Inverse();

  friend Permutation operator*(const Permutation& lhs, const Permutation& rhs);

  [[nodiscard]] std::string ToString() const;

  [[nodiscard]] std::string ToLaTex() const;

  friend std::ostream& operator<<(std::ostream& out, const Permutation& permutation);

 private:
  std::vector<int> permutation_;
};

class Matrix {
 public:
  Matrix() = default;

  Matrix(int n, int m, int item = 0);

  explicit Matrix(std::vector<std::vector<Fraction>> elements);

  explicit Matrix(const Permutation& permutation);

  [[nodiscard]] Matrix GetTranspose() const;

  Fraction& At(int i, int j);

  [[nodiscard]] const Fraction& At(int i, int j) const;

  [[nodiscard]] int GetNumRows() const;

  [[nodiscard]] int GetNumColumns() const;

  Matrix& operator*=(const Permutation& permutation);

  friend bool operator==(const Matrix& lhs, const Fraction& rhs);

  friend bool operator==(const Fraction& lhs, const Matrix& rhs);

  friend bool operator==(const Matrix& lhs, const Matrix& rhs);

  friend Matrix operator+(const Matrix& lhs, const Matrix& rhs);

  friend Matrix operator*(const Matrix& lhs, const Matrix& rhs);

  friend Matrix operator*(const Permutation& lhs, const Matrix& rhs);

  friend Matrix operator*(const Matrix& lhs, const Permutation& rhs);

  friend std::istream& operator>>(std::istream& in, Matrix& matrix);

  [[nodiscard]] std::string ToString() const;

  [[nodiscard]] std::string ToLaTex() const;

  friend std::ostream& operator<<(std::ostream& out, const Matrix& matrix);

  bool AddRow(int first_row, int second_row, Fraction coefficient = 1);

  bool SubtractRow(int minuend_row, int subtrahend_row, Fraction coefficient = 1);

  bool MultiplyRow(int row, Fraction coefficient);

  bool DivideRow(int row, Fraction coefficient);

  bool SwapRows(int first_row, int second_row);

  bool CopyRow(int source_row, int destination_row);

  bool CopyRow(const Matrix& source_matrix, int source_row, int destination_row);

  void ApplyRowPermutation(const Permutation& permutation);

  bool AddColumn(int first_column, int second_column, Fraction coefficient = 1);

  bool SubtractColumn(int minuend_column, int subtrahend_column, Fraction coefficient = 1);

  bool MultiplyColumn(int column, Fraction coefficient);

  bool DivideColumn(int column, Fraction coefficient);

  bool SwapColumns(int first_column, int second_column);

  bool CopyColumn(int source_column, int destination_column);

  bool CopyColumn(const Matrix& source_matrix, int source_column, int destination_column);

  void ApplyColumnPermutation(const Permutation& permutation);

  [[nodiscard]] bool IsLowerTriangular() const;

  [[nodiscard]] bool IsUpperTriangular() const;

  [[nodiscard]] bool IsDiagonal() const;

  [[nodiscard]] Matrix GetInverse() const;

  bool Inverse();

  [[nodiscard]] LUPMatrix GetLUPTransform(
      Options::ChoiceType choice_type = Options::ChoiceType::Column) const;

  static Matrix Identity(int n);

 private:
  int n_{}, m_{};
  std::vector<std::vector<Fraction>> elements_;
};

class LUPMatrix {
 public:
  LUPMatrix() = default;

  LUPMatrix(Matrix L, Matrix U, const Permutation& row_permutation, const Permutation& column_permutation);

  explicit LUPMatrix(const ExtendedMatrix& system);

  [[nodiscard]] const Matrix& GetLowerMatrix() const;

  [[nodiscard]] const Matrix& GetUpperMatrix() const;

  [[nodiscard]] const Permutation& GetRowPermutation() const;

  [[nodiscard]] const Permutation& GetColumnPermutation() const;

  [[nodiscard]] Matrix GetMatrix() const;

  [[nodiscard]] std::string ToString() const;

  [[nodiscard]] std::string ToLaTex() const;

  friend std::ostream& operator<<(std::ostream& out, const LUPMatrix& lup_matrix);

 private:
  Matrix l_, u_;
  Permutation row_permutation_, column_permutation_;
};

class ExtendedMatrix {
 public:
  ExtendedMatrix(Matrix a, Matrix b);

  bool MultiplyRow(int row, Fraction coefficient);

  bool AddRow(int first_row, int second_row, Fraction coefficient = 1);

  bool SwapRows(int first_row, int second_row);

  bool SwapColumns(int first_column, int second_column);

  [[nodiscard]] const Matrix& GetMainMatrix() const;

  [[nodiscard]] const Matrix& GetExtensionMatrix() const;

  [[nodiscard]] const Permutation& GetRowPermutation() const;

  [[nodiscard]] const Permutation& GetMainColumnPermutation() const;

  [[nodiscard]] const Permutation& GetExtensionColumnPermutation() const;

  [[nodiscard]] const Fraction& At(int row, int column) const;

  [[nodiscard]] int GetNumRows() const;

  [[nodiscard]] int GetNumColumns() const;

  [[nodiscard]] std::string ToString() const;

  [[nodiscard]] std::string ToLaTex() const;

  friend std::ostream& operator<<(std::ostream& out, const ExtendedMatrix& matrix);

 private:
  Matrix main_, extension_;
  Permutation row_permutation_;
  Permutation main_column_permutation_, extension_column_permutation_;
};

class LinearSystem {
 public:
  enum class Mode {
    Standard,
    LU,
  };

  class SwapMatrix {
   public:
    SwapMatrix();

    void SetRowSwap(int first_row, int second_row);

    void SetMainColumnSwap(int first_column, int second_column);

    void SetExtensionColumnSwap(int first_column, int second_column);

    [[nodiscard]] bool IsEmpty() const;

    void Clear();

    [[nodiscard]] std::string ToString() const;

    [[nodiscard]] std::string ToLaTex() const;

    friend std::ostream& operator<<(std::ostream& out, const SwapMatrix& swaps);

   private:
    int first_row_, second_row_;
    int first_main_column_, second_main_column_;
    int first_extension_column_, second_extension_column_;
  };

  LinearSystem(const Matrix& a, const Matrix& b, Mode mode = Mode::Standard);

  explicit LinearSystem(const ExtendedMatrix& system, Mode mode = Mode::Standard);

  bool RunDirectGauss(Options::ChoiceType choice_type = Options::ChoiceType::Column);

  void RunReverseGauss();

  bool RunGauss(Options::ChoiceType choice_type = Options::ChoiceType::Column);

  [[nodiscard]] Matrix GetSolutionMatrix() const;

  [[nodiscard]] const ExtendedMatrix& GetExtendedMatrix() const;

  [[nodiscard]] LUPMatrix GetLUPMatrix() const;

  [[nodiscard]] std::string ToString() const;

  [[nodiscard]] std::string ToLaTex() const;

  friend std::ostream& operator<<(std::ostream& out, const LinearSystem& system);

 private:
  enum class OperationType {
    AddRow,
    MultiplyRow,
    SwapRows,
    SwapColumns,
    MainOutput
  };

  struct Operation {
    Operation() = default;

    Operation(OperationType operation_type, int lhs, int rhs, Fraction coefficient);

    OperationType operation_type = OperationType::MainOutput;
    int lhs = 0, rhs = 0;
    Fraction coefficient = 0;
  };

  bool AddMainOutput();

  bool PerformOperation(OperationType operation_type, int lhs, int rhs, Fraction coefficient);

  bool Choice(int pos, Options::ChoiceType choice_type = Options::ChoiceType::Column);

  bool ApplyOperation(OperationType operation_type, int lhs, int rhs, Fraction coefficient);

  bool ApplyOperation(const Operation& operation);

  Mode mode_;
  int n_, m_, rank_;
  const ExtendedMatrix start_system_;
  ExtendedMatrix system_;
  std::vector<Operation> history_;
};
