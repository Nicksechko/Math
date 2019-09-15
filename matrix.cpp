#include "matrix.h"

Matrix Matrix::ToTranspose() const {
  Matrix result(m_, n_);
  for (int i = 0; i < m_; ++i) {
    for (int j = 0; j < n_; ++j) {
      result.At(i, j) = At(j, i);
    }
  }

  return result;
}

std::string Matrix::ToLatex() {
  std::ostringstream out;
  out << "\\begin{bmatrix}" << std::endl;
  for (int i = 0; i < n_; ++i) {
    for (int j = 0; j < m_; ++j) {
      out << elements_[i][j].ToLaTex();
      if (j + 1 != m_) {
        out << " & ";
      }
    }
    out << " \\\\" << std::endl;
  }
  out << "\\end{bmatrix}" << std::endl;
  return out.str();
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
  if (Options::step_by_step_type != Options::StepByStepType::Without) {
    std::cout << "Inverse Matrix" << std::endl << std::endl;
  }
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
  int rank = transformer.RunDirectGauss(Options::ChoiceType::Row);
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
