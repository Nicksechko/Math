#include <iostream>

#include "matrix.h"

int main() {
  Options::step_by_step_type = Options::StepByStepType::MainSteps;
  Options::output_type = Options::OutputType::LaTex;

  // Fraction::Test();

  int n;
  std::cin >> n;

  Matrix a({{n + 1, n / 2, -n / 2, 1},
            {-n - 1, -(n + 1) / 2, (n + 1) / 2, -(n + 2) / 3},
            {-n + 1, (n + 1) / 2, -(n + 2) / 3, n - 1},
            {n / 3, -1, n, -n}});

  LUPMatrix lup_matrix;
  if (a.LUPTransform(lup_matrix)) {
    std::cout << lup_matrix.ToLatex() << std::endl;
  } else {
    std::cout << "LUP Transform doesn't exist";
  }

  std::cout << (lup_matrix.L * lup_matrix.U).ToLatex() << std::endl;

//  std::cout << Matrix::Identity(2);

  return 0;
}
