#include "matrix.h"

ExtendedMatrix::ExtendedMatrix(Matrix main, Matrix extension)
  : main_(std::move(main)),
    extension_(std::move(extension)),
    row_permutation_(main.GetNumRows()),
    main_column_permutation_(main.GetNumColumns()),
    extension_column_permutation_(extension.GetNumColumns()){
  assert(main_.GetNumRows() == extension_.GetNumRows());
}

bool ExtendedMatrix::AddRow(int first_row, int second_row, Fraction coefficient) {
  bool result = false;
  result |= main_.AddRow(first_row, second_row, coefficient);
  result |= extension_.AddRow(first_row, second_row, coefficient);

  return result;
}

bool ExtendedMatrix::MultiplyRow(int row, Fraction coefficient) {
  bool result = false;
  result |= main_.MultiplyRow(row, coefficient);
  result |= extension_.MultiplyRow(row, coefficient);

  return result;
}

bool ExtendedMatrix::SwapRows(int first_row, int second_row) {
  bool result = false;
  result |= main_.SwapRows(first_row, second_row);
  result |= extension_.SwapRows(first_row, second_row);

  if (result) {
    row_permutation_.Transposition(first_row, second_row);
  }

  return result;
}

bool ExtendedMatrix::SwapColumns(int first_column, int second_column) {
  int m = main_.GetNumColumns();
  bool result = false;
  if (first_column < m && second_column < m) {
    result = main_.SwapColumns(first_column, second_column);
    if (result) {
      main_column_permutation_.Transposition(first_column, second_column);
    }
  } else if (first_column >= m && second_column >= m) {
    result = extension_.SwapColumns(first_column - m, second_column - m);
    if (result) {
      extension_column_permutation_.Transposition(first_column - m, second_column - m);
    }
  }

  return result;
}

const Matrix& ExtendedMatrix::GetMainMatrix() const {
  return main_;
}

const Matrix& ExtendedMatrix::GetExtensionMatrix() const {
  return extension_;
}

const Permutation& ExtendedMatrix::GetRowPermutation() const {
  return row_permutation_;
}

const Permutation& ExtendedMatrix::GetMainColumnPermutation() const {
  return main_column_permutation_;
}

const Permutation& ExtendedMatrix::GetExtensionColumnPermutation() const {
  return extension_column_permutation_;
}

const Fraction& ExtendedMatrix::At(int row, int column) const {
  if (column < main_.GetNumColumns()) {
    return main_.At(row, column);
  } else {
    return extension_.At(row, column - main_.GetNumColumns());
  }
}

int ExtendedMatrix::GetNumRows() const {
  return main_.GetNumRows();
}

int ExtendedMatrix::GetNumColumns() const {
  return main_.GetNumColumns() + extension_.GetNumColumns();
}

std::ostream& operator<<(std::ostream& out, const ExtendedMatrix& matrix) {
  if (Options::output_type == Options::OutputType::Standard) {
    out << matrix.ToString();
  } else if (Options::output_type == Options::OutputType::LaTex){
    out << matrix.ToLaTex();
  }

  return out;
}

std::string ExtendedMatrix::ToString() const {
  std::ostringstream out;
  out << GetNumRows() << ' ' << GetNumColumns() << std::endl;
  for (int row = 0; row < GetNumRows(); ++row) {
    for (int column = 0; column < GetNumColumns(); ++column) {
      out << At(row, column) << ' ';
    }
    if (row + 1 != GetNumRows()) {
      out << std::endl;
    }
  }

  return out.str();
}

std::string ExtendedMatrix::ToLaTex() const {
  std::ostringstream out;
  out << "\\begin{bmatrix}" << std::endl;
  for (int row = 0; row < GetNumRows(); ++row) {
    for (int column = 0; column < main_.GetNumColumns(); ++column) {
      out << main_.At(row, column).ToLaTex() << " & ";
    }
    out << "\\vline & ";
    for (int column = 0; column < extension_.GetNumColumns(); ++column) {
      out << extension_.At(row, column).ToLaTex();
      if (column + 1 != extension_.GetNumColumns()) {
        out << " & ";
      }
    }
    if (row + 1 != GetNumRows()) {
      out << "\\\\";
    }
    out << std::endl;
  }
  out << "\\end{bmatrix}";

  return out.str();
}
