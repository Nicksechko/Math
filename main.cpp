#include <iostream>

#include "matrix.h"
#include "fraction.h"

int main() {
  // Fraction::Test();
  Matrix<Fraction> one({
    {1, 2, 4}
  });
  Matrix<Fraction> two({
    {1, 2, {1, 2}}
  });

  std::cout << one + two << std::endl;
  std::cout << one * two.ToTranspose() << std::endl;

  return 0;
}
