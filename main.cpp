#include <iostream>
#include <vector>
#include <algorithm>
#include <numeric>

#include "matrix.h"

int main() {
  Options::step_by_step_type = Options::StepByStepType::MainSteps;
  Options::output_type = Options::OutputType::LaTex;

  int n;
  std::cin >> n;

//  Matrix a({{1, 2}, {3, 4}});

  Matrix a({{n + 1, n / 2, -n / 2, 1},
            {-n - 1, -(n + 1) / 2, (n + 1) / 2, -(n + 2) / 3},
            {-n + 1, (n + 1) / 2, -(n + 2) / 3, n - 1},
            {n / 3, -1, n, -n}});


  LUPMatrix lup_matrix = a.GetLUPTransform(Options::ChoiceType::Submatrix);

  Options::writer << lup_matrix << std::endl;

  std::cout << Options::writer.str() << std::endl;

//  std::cout << a << std::endl;
//
  std::cout << lup_matrix.GetMatrix() << std::endl;

  assert(a == lup_matrix.GetMatrix());

  return 0;
}
