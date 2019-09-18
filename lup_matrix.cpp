#include "matrix.h"

LUPMatrix::LUPMatrix(Matrix l, Matrix u,
                     const Permutation& row_permutation,
                     const Permutation& column_permutation)
    : l_(std::move(l)),
      u_(std::move(u)),
      row_permutation_(row_permutation),
      column_permutation_(column_permutation) {
  assert(l_.IsLowerTriangular());
  assert(u_.IsUpperTriangular());
}

LUPMatrix::LUPMatrix(const ExtendedMatrix& system)
  : LUPMatrix(system.GetExtensionMatrix().GetInverse(), system.GetMainMatrix(),
              system.GetRowPermutation(), system.GetMainColumnPermutation()) {}

const Matrix& LUPMatrix::GetLowerMatrix() const {
  return l_;
}

const Matrix& LUPMatrix::GetUpperMatrix() const {
  return u_;
}

const Permutation& LUPMatrix::GetRowPermutation() const {
  return row_permutation_;
}

const Permutation& LUPMatrix::GetColumnPermutation() const {
  return column_permutation_;
}

Matrix LUPMatrix::GetMatrix() const {
  return row_permutation_.GetInverse() * l_ * u_ * column_permutation_;
}

std::string LUPMatrix::ToString() const {
  std::ostringstream out;
  out << "L: " << std::endl;
  out << l_ << std::endl;
  out << "U: " << std::endl;
  out << u_;
  out << "Row Permutation: " << std::endl;
  out << row_permutation_ << std::endl;
  out << "Column Permutation: " << std::endl;
  out << column_permutation_;

  return out.str();
}

std::string LUPMatrix::ToLaTex() const {
  std::ostringstream out;
  out << "\\[" << std::endl;
  out << "L: " << l_ << "\\quad" << std::endl;
  out << "U: " << u_;
  out << "\\]" << std::endl;
  out << "\\[" << std::endl;
  out << "Row Permutation: " << std::endl;
  out << row_permutation_ << "\\quad" << std::endl;
  out << "Column Permutation: " << std::endl;
  out << column_permutation_ << std::endl;
  out << "\\]" << std::endl;

  return out.str();
}

std::ostream& operator<<(std::ostream& out, const LUPMatrix& lup_matrix) {
  if (Options::output_type == Options::OutputType::Standard) {
    out << lup_matrix.ToString() << std::endl;
  } else if (Options::output_type == Options::OutputType::LaTex) {
    out << lup_matrix.ToLaTex() << std::endl;
  }

  return out;
}
