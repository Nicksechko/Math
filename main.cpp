#include <iostream>
#include <vector>
#include <algorithm>
#include <numeric>
#include <chrono>

#include "matrix.h"
#include <random>
#include <chrono>
#include <fstream>

void HW2(int n) {
    Matrix a({{n + 1, n / 2, -n / 2, 1},
              {-n - 1, -(n + 1) / 2, (n + 1) / 2, -(n + 2) / 3},
              {-n + 1, (n + 1) / 2, -(n + 2) / 3, n - 1},
              {n / 3, -1, n, -n}});


    Options::writer << "n = " << n << std::endl;
    Options::writer << "\\[" << std::endl;
    Options::writer << a << std::endl;
    Options::writer << "\\]" << std::endl;

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

    std::ofstream fout("output");
    fout << Options::writer.str() << std::endl;

}

void HW3() {
    std::mt19937 twister(std::chrono::steady_clock::now().time_since_epoch().count());
    std::uniform_int_distribution<int> dis(-3, 3);
//    Matrix a({{dis(twister), dis(twister), dis(twister)},
//                       {dis(twister), dis(twister), dis(twister)},
//                       {dis(twister), dis(twister), dis(twister)}});
    Matrix a({{0, -1, 0},
                        {-2, 3, 3},
                        {0, -1, 3}});


    Options::writer << "\\[" << std::endl;
    Options::writer << a << std::endl;
    Options::writer << "\\]" << std::endl;

    LUPMatrix lup_matrix = a.GetLUPTransform(Options::ChoiceType::Submatrix);

    Options::writer << "\\[" << std::endl;
    Options::writer << lup_matrix << std::endl;
    Options::writer << "\\]" << std::endl;

    Matrix inverse = lup_matrix.SolveSystem(Matrix::Identity(a.GetNumRows()));

    assert(inverse * a == Matrix::Identity(a.GetNumRows()));

    assert(a == lup_matrix.GetMatrix());

    Options::writer << std::endl;

    std::ofstream fout("output");
    fout << Options::writer.str() << std::endl;
}

int main() {
  Options::step_by_step_type = Options::StepByStepType::MainSteps;
  Options::output_type = Options::OutputType::LaTex;

  HW3();

  return 0;
}
