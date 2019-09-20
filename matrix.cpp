#include "matrix.h"

Matrix::Matrix(int n, int m, int item)
    : n_(n), m_(m),
      elements_(std::vector<std::vector<Fraction>>(n, std::vector<Fraction>(m, item))) {}

Matrix::Matrix(std::vector<std::vector<Fraction>> elements)
    : n_(elements.size()), m_((n_ == 0) ? 0 : static_cast<int>(elements[0].size())),
      elements_(std::move(elements)) {

}

Matrix::Matrix(const Permutation& permutation)
    : n_(permutation.Size()), m_(permutation.Size()),
      elements_(std::vector<std::vector<Fraction>>(
          permutation.Size(), std::vector<Fraction>(permutation.Size()))) {
  for (int i = 0; i < n_; i++) {
    elements_[i][permutation.At(i)] = 1;
  }
}

Matrix Matrix::GetTranspose() const {
  Matrix result(m_, n_);
  for (int i = 0; i < m_; ++i) {
    for (int j = 0; j < n_; ++j) {
      result.At(i, j) = At(j, i);
    }
  }

  return result;
}

Fraction& Matrix::At(int i, int j) {
  return elements_.at(i).at(j);
}

const Fraction& Matrix::At(int i, int j) const {
  return elements_.at(i).at(j);
}

int Matrix::GetNumRows() const {
  return n_;
}

int Matrix::GetNumColumns() const {
  return m_;
}

Matrix& Matrix::operator*=(const Permutation& permutation) {
  ApplyColumnPermutation(permutation.GetInverse());
  return *this;
}

bool operator==(const Matrix& lhs, const Fraction& rhs) {
  for (int i = 0; i < lhs.GetNumRows(); ++i) {
    for (int j = 0; j < lhs.GetNumRows(); ++j) {
      if (lhs.At(i, j) != rhs) {
        return false;
      }
    }
  }

  return true;
}

bool operator==(const Fraction& lhs, const Matrix& rhs) {
  return operator==(rhs, lhs);
}

bool operator==(const Matrix& lhs, const Matrix& rhs) {
  if (lhs.GetNumRows() != rhs.GetNumRows() ||
      lhs.GetNumColumns() != rhs.GetNumColumns()) {
    return false;
  }

  for (int i = 0; i < lhs.GetNumRows(); ++i) {
    for (int j = 0; j < lhs.GetNumColumns(); ++j) {
      if (lhs.At(i, j) != rhs.At(i, j)) {
        return false;
      }
    }
  }

  return true;
}

Matrix operator+(const Matrix& lhs, const Matrix& rhs) {
  if (lhs.GetNumRows() != rhs.GetNumRows()) {
    throw std::invalid_argument("[SUM] Mismatched number of rows");
  }

  if (lhs.GetNumColumns() != rhs.GetNumColumns()) {
    throw std::invalid_argument("[SUM] Mismatched number of columns");
  }

  Matrix result(lhs.GetNumRows(), lhs.GetNumColumns());
  for (int i = 0; i < result.GetNumRows(); ++i) {
    for (int j = 0; j < result.GetNumColumns(); ++j) {
      result.At(i, j) = lhs.At(i, j) + rhs.At(i, j);
    }
  }

  return result;
}

Matrix operator*(const Matrix& lhs, const Matrix& rhs) {
  if (lhs.GetNumColumns() != rhs.GetNumRows()) {
    throw std::invalid_argument("[MUL] Mismatched number of rows");
  }

  Matrix result(lhs.GetNumRows(), rhs.GetNumColumns());
  for (int i = 0; i < result.GetNumRows(); ++i) {
    for (int j = 0; j < result.GetNumColumns(); ++j) {
      for (int k = 0; k < lhs.GetNumColumns(); ++k) {
        result.At(i, j) += lhs.At(i, k) * rhs.At(k, j);
      }
    }
  }

  return result;
}

Matrix operator*(const Permutation& lhs, const Matrix& rhs) {
  Matrix result(rhs);
  result.ApplyRowPermutation(lhs);
  return result;
}

Matrix operator*(const Matrix& lhs, const Permutation& rhs) {
  Matrix result(lhs);
  result *= rhs;
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

std::string Matrix::ToString() const {
  std::ostringstream out;
  for (int i = 0; i < GetNumRows(); ++i) {
    for (int j = 0; j < GetNumColumns(); ++j) {
      if (j > 0) {
        out << " ";
      }
      out << elements_[i][j];
    }
    if (i + 1 != GetNumRows()) {
      out << std::endl;
    }
  }

  return out.str();
}

std::string Matrix::ToLaTex() const {
  std::ostringstream out;
  out << "\\begin{bmatrix}" << std::endl;
  for (int i = 0; i < n_; ++i) {
    for (int j = 0; j < m_; ++j) {
      out << elements_[i][j].ToLaTex();
      if (j + 1 != m_) {
        out << " & ";
      }
    }
    if (i + 1 != n_) {
      out << " \\\\";
    }
    out << std::endl;
  }
  out << "\\end{bmatrix}";

  return out.str();
}

std::ostream& operator<<(std::ostream& out, const Matrix& matrix) {
  if (Options::output_type == Options::OutputType::Standard) {
    out << matrix.ToString();
  } else if (Options::output_type == Options::OutputType::LaTex) {
    out << matrix.ToLaTex();
  }

  return out;
}

bool Matrix::AddRow(int first_row, int second_row, Fraction coefficient) {
  if (coefficient == 0) {
    return false;
  }

  bool result = false;
  for (int column = 0; column < m_; ++column) {
    if (elements_[second_row][column] != 0) {
      result = true;
      elements_[first_row][column] += coefficient * elements_[second_row][column];
    }
  }

  return result;
}

bool Matrix::SubtractRow(int minuend_row, int subtrahend_row, Fraction coefficient) {
  return AddRow(minuend_row, subtrahend_row, -coefficient);
}

bool Matrix::MultiplyRow(int row, Fraction coefficient) {
  if (coefficient == 1) {
    return false;
  }

  bool result = false;
  for (int column = 0; column < m_; ++column) {
    if (elements_[row][column] != 0) {
      result = true;
      elements_[row][column] *= coefficient;
    }
  }

  return result;
}

bool Matrix::DivideRow(int row, Fraction coefficient) {
  return MultiplyRow(row, 1 / coefficient);
}

bool Matrix::SwapRows(int first_row, int second_row) {
  if (first_row == second_row) {
    return false;
  }

  bool result = false;
  for (int column = 0; column < m_; ++column) {
    if (elements_[first_row][column] != elements_[second_row][column]) {
      result = true;
      std::swap(elements_[first_row][column], elements_[second_row][column]);
    }
  }

  return result;
}

bool Matrix::CopyRow(int source_row, int destination_row) {
  if (source_row == destination_row) {
    return false;
  }

  bool result = false;
  for (int column = 0; column < m_; ++column) {
    if (elements_[destination_row][column] != elements_[source_row][column]) {
      result = true;
      elements_[destination_row][column] = elements_[source_row][column];
    }
  }

  return result;
}

bool Matrix::CopyRow(const Matrix& source_matrix, int source_row, int destination_row) {
  bool result = false;
  for (int column = 0; column < m_; ++column) {
    if (elements_[destination_row][column] != source_matrix.At(source_row, column)) {
      result = true;
      elements_[destination_row][column] = source_matrix.At(source_row, column);
    }
  }

  return result;
}

bool Matrix::ApplyRowPermutation(const Permutation& permutation) {
  assert(permutation.Size() == n_);
  Matrix result(elements_);
  bool answer = false;
  for (int i = 0; i < n_; ++i) {
    answer |= result.CopyRow(*this, permutation.At(i), i);
  }
  operator=(result);

  return answer;
}

bool Matrix::AddColumn(int first_column, int second_column, Fraction coefficient) {
  if (coefficient == 0) {
    return false;
  }

  bool result = false;
  for (int row = 0; row < n_; ++row) {
    if (elements_[row][second_column] != 0) {
      result = true;
      elements_[row][first_column] += coefficient * elements_[row][second_column];
    }
  }

  return result;
}

bool Matrix::SubtractColumn(int minuend_column, int subtrahend_column, Fraction coefficient) {
  return AddColumn(minuend_column, subtrahend_column, -coefficient);
}

bool Matrix::MultiplyColumn(int column, Fraction coefficient) {
  if (coefficient == 1) {
    return false;
  }

  int result = false;
  for (int row = 0; row < n_; ++row) {
    if (elements_[row][column] != 0) {
      result = true;
      elements_[row][column] *= coefficient;
    }
  }

  return result;
}

bool Matrix::DivideColumn(int column, Fraction coefficient) {
  return MultiplyColumn(column, 1 / coefficient);
}

bool Matrix::SwapColumns(int first_column, int second_column) {
  if (first_column == second_column) {
    return false;
  }

  int result = false;
  for (int row = 0; row < n_; ++row) {
    if (elements_[row][first_column] != elements_[row][second_column]) {
      result = true;
      std::swap(elements_[row][first_column], elements_[row][second_column]);
    }
  }

  return result;
}

bool Matrix::CopyColumn(int source_column, int destination_column) {
  if (source_column == destination_column) {
    return false;
  }

  bool result = false;
  for (int row = 0; row < n_; ++row) {
    if (elements_[row][destination_column] != elements_[row][source_column]) {
      result = true;
      elements_[row][destination_column] = elements_[row][source_column];
    }
  }

  return result;
}

bool Matrix::CopyColumn(const Matrix& source_matrix, int source_column, int destination_column) {
  bool result = false;
  for (int row = 0; row < n_; ++row) {
    if (elements_[row][destination_column] != source_matrix.At(row, source_column)) {
      result = true;
      elements_[row][destination_column] = source_matrix.At(row, source_column);
    }
  }

  return result;
}

bool Matrix::ApplyColumnPermutation(const Permutation& permutation) {
  assert(permutation.Size() == m_);
  Matrix result(elements_);
  bool answer = false;
  for (int column = 0; column < m_; ++column) {
    answer |= result.CopyColumn(*this, permutation.At(column), column);
  }
  operator=(result);

  return answer;
}

bool Matrix::IsLowerTriangular() const {
  for (int row = 0; row < n_; ++row) {
    for (int column = row + 1; column < m_; ++column) {
      if (elements_[row][column] != 0) {
        return false;
      }
    }
  }

  return true;
}

bool Matrix::IsUpperTriangular() const {
  for (int column = 0; column < m_; ++column) {
    for (int row = column + 1; row < n_; ++row) {
      if (elements_[row][column] != 0) {
        return false;
      }
    }
  }

  return true;
}

bool Matrix::IsDiagonal() const {
  for (int row = 0; row < n_; ++row) {
    for (int column = 0; column < m_; ++column) {
      if (row != column && elements_[row][column] != 0) {
        return false;
      }
    }
  }

  return true;
}

Matrix Matrix::GetInverse(Options::ChoiceType choice_type) const {
  assert(n_ == m_);
  LinearSystem inverser(*this, Identity(n_));
  bool linear_independent = inverser.RunGauss(choice_type);
  Options::writer << inverser << std::endl;
  if (!linear_independent) {
    return Matrix(n_, n_);
  }

  return inverser.GetSolutionMatrix();
}

bool Matrix::Inverse(Options::ChoiceType choice_type) {
  Matrix result = GetInverse(choice_type);
  if (result == 0) {
    return false;
  } else {
    operator=(result);
    return true;
  }
}

LUPMatrix Matrix::GetLUPTransform(Options::ChoiceType choice_type) const {
  assert(n_ == m_);
  LinearSystem transformer(*this, Identity(n_), LinearSystem::Mode::LU);
  bool linear_independent = transformer.RunDirectGauss(choice_type);
  Options::writer << transformer << std::endl;
  if (!linear_independent) {
    return LUPMatrix();
  }

  return transformer.GetLUPMatrix();
}

Matrix Matrix::Identity(int n) {
  Matrix identity(n, n);
  for (int i = 0; i < n; i++) {
    identity.At(i, i) = 1;
  }

  return identity;
}

Fraction Matrix::GetOctahedralNorm() const {
    Fraction result = 0;
    for (int column = 0; column < m_; ++column) {
        Fraction sum = 0;
        for (int row = 0; row < n_; ++row) {
            sum += fabs(elements_[row][column]);
        }
        result = std::max(result, sum);
    }

    return result;
}

Fraction Matrix::GetCubicNorm() const {
    Fraction result;
    for (int row = 0; row < m_; ++row) {
        Fraction sum = 0;
        for (int column = 0; column < n_; ++column) {
            sum += fabs(elements_[row][column]);
        }
        result = std::max(result, sum);
    }

    return result;
}

double Matrix::GetFrobeniusNorm() const {
    Fraction result;
    for (const auto& row : elements_) {
        for (const auto& item : row) {
            result += item * item;
        }
    }

    return sqrt(result.ToDouble());
}
