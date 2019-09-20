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
  : LUPMatrix(system.GetExtensionMatrix().GetInverse(Options::ChoiceType::Without), system.GetMainMatrix(),
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

Matrix LUPMatrix::SolveSystem(const Matrix& rhs) const {
    Matrix result = rhs;

    if (result.ApplyRowPermutation(row_permutation_)) {
        Options::writer << "\\[" << std::endl;
        Options::writer << row_permutation_.AsMatrix() << "\\ *\\" << std::endl;
        Options::writer << rhs << "\\ =\\" << std::endl;
        Options::writer << result << std::endl;
        Options::writer << "\\]" << std::endl;
    }

    LinearSystem lower_solver(l_, result);
    lower_solver.RunGauss(Options::ChoiceType::Without);
    Options::writer << lower_solver << std::endl;
    result = lower_solver.GetSolutionMatrix();

    LinearSystem upper_solver(u_, result);
    upper_solver.RunGauss(Options::ChoiceType::Without);
    Options::writer << upper_solver << std::endl;
    result = upper_solver.GetSolutionMatrix();

    Matrix answer = result;
    if (answer.ApplyRowPermutation(column_permutation_.GetInverse())) {
        Options::writer << "\\[" << std::endl;
        Options::writer << column_permutation_.GetInverse().AsMatrix() << "\\ *\\" << std::endl;
        Options::writer << result << "\\ =\\" << std::endl;
        Options::writer << answer << std::endl;
        Options::writer << "\\]";
    }

    return answer;
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
  if (!row_permutation_.IsIdentity()) {
      out << row_permutation_.GetInverse().AsMatrix() << " \\ *\\" << std::endl;
  }
  out << l_ << " \\ *\\" << std::endl;
  out << u_;
  if (!column_permutation_.IsIdentity()) {
      out << "\\ *\\" << std::endl;
      out << column_permutation_.AsMatrix();
  }

  return out.str();
}

std::ostream& operator<<(std::ostream& out, const LUPMatrix& lup_matrix) {
  if (Options::output_type == Options::OutputType::Standard) {
    out << lup_matrix.ToString();
  } else if (Options::output_type == Options::OutputType::LaTex) {
    out << lup_matrix.ToLaTex();
  }

  return out;
}
