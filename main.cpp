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


  Options::writer << "n = " << n << std::endl;
  Options::writer << "\\[" << std::endl;
  Options::writer << a << std::endl;
  Options::writer << "\\]" << std::endl;
//  std::cout << a.GetOctahedralNorm() << std::endl;

//  std::cout << a.GetCubicNorm() << std::endl;

//  std::cout << a.GetFrobeniusNorm() << std::endl;

  LUPMatrix lup_matrix = a.GetLUPTransform(Options::ChoiceType::Row);

  Options::writer << "\\[" << std::endl;
  Options::writer << lup_matrix << std::endl;
  Options::writer << "\\]" << std::endl;

  Matrix inverse = lup_matrix.SolveSystem(Matrix::Identity(a.GetNumRows()));

  assert(inverse * a == Matrix::Identity(a.GetNumRows()));

  assert(a == lup_matrix.GetMatrix());

  Options::writer << std::endl;
  Options::writer << "\\[" << std::endl;
  Options::writer << "\\mu_{1} = " << a.GetCubicNorm() * inverse.GetCubicNorm() << std::endl;
  Options::writer << "\\]" << std::endl;
  Options::writer << "\\[" << std::endl;
  Options::writer << "\\mu_{\\infty} = " << a.GetOctahedralNorm() * inverse.GetOctahedralNorm() << std::endl;
  Options::writer << "\\]" << std::endl;
  Options::writer << "\\[" << std::endl;
  Options::writer << "\\mu_{F} = " << a.GetFrobeniusNorm() * inverse.GetFrobeniusNorm() << std::endl;
  Options::writer << "\\]" << std::endl;

  std::cout << Options::writer.str() << std::endl;

//  std::cout << a << std::endl;


  return 0;
}
