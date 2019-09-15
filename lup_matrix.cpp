#include "matrix.h"

LUPMatrix::LUPMatrix(Matrix L,
                     Matrix U,
                     std::vector<int> row_permutation,
                     std::vector<int> column_permutation)
    : L(std::move(L)),
      U(std::move(U)),
      row_permutation(std::move(row_permutation)),
      column_permutation(std::move(column_permutation)) {}

std::string LUPMatrix::ToLatex() {
  std::ostringstream out;
  out << "$$" << std::endl;
  out << "L: " << L.ToLatex();
  out << "U: " << U.ToLatex();
  out << "$$" << std::endl;
  out << "$$" << std::endl;
  out << "Permutation: " << std::endl;
  for (int item : column_permutation) {
    out << item << "\\ ";
  }
  out << std::endl;
  out << "$$" << std::endl;
  return out.str();
}


