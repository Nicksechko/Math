#include <iostream>

#include "matrix.h"
#include "fraction.h"

int main() {
  // Fraction::Test();
  Matrix<Fraction> one({
    {{1, 2}}
  });
  Matrix<Fraction> two;

  cout << one + two << endl;

  return 0;
}
