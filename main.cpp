#include <iostream>

#include "matrix.h"
#include "fraction.h"

int main() {
  // Fraction::Test();
  Matrix<Fraction> a({
    {1, 2},
    {3, 4},
  });

  LUPMatrix<Fraction> lup_a;
  if (a.LUPTransform(lup_a)) {
    std::cout << lup_a.L << std::endl;
    std::cout << lup_a.U << std::endl;
    for (int item : lup_a.row_permutation) {
      std::cout << item << " ";
    }
    std::cout << std::endl << std::endl;
    for (int item : lup_a.column_permutation) {
      std::cout << item << " ";
    }
    std::cout << std::endl << std::endl;
  }

  return 0;
}
