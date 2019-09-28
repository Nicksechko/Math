#include <iostream>
#include <vector>
#include <algorithm>
#include <numeric>
#include <chrono>

#include "matrix.h"
#include <random>
#include <chrono>
#include <fstream>
#include <unordered_set>
#include <set>

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
                        {0, -2, 3}});


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

void SolveLDLT(const Matrix& l, const Matrix& d, const Matrix& b) {
    LinearSystem z_solver(l, b);
    z_solver.RunGauss(Options::ChoiceType::Without);
    Options::writer << z_solver << std::endl;
    Matrix z = z_solver.GetSolutionMatrix();

    LinearSystem y_solver(d, z);
    y_solver.RunGauss(Options::ChoiceType::Without);
    Options::writer << y_solver << std::endl;
    Matrix y = y_solver.GetSolutionMatrix();

    LinearSystem solver(l.GetTranspose(), y);
    solver.RunGauss(Options::ChoiceType::Without);
    Options::writer << solver << std::endl;
}

void HW4() {
    std::mt19937 twister(8888);
//    std::random_device true_random;
//    std::mt19937 twister(true_random());
    std::uniform_int_distribution<int> dis(-4, 4);
    const int max_iterations = 10'000;

    Matrix a;

    std::ofstream fout("output");
    int iterations = max_iterations;
    while (--iterations) {
        Matrix l({{1, 0, 0, 0},
                  {dis(twister), 1, 0, 0},
                  {dis(twister), dis(twister), 1, 0},
                  {dis(twister), dis(twister), dis(twister), 1}});
        Matrix d({{dis(twister), 0, 0, 0},
                  {0, dis(twister), 0, 0},
                  {0, 0, dis(twister), 0},
                  {0, 0, 0, dis(twister)}});
        a = l * d * l.GetTranspose();

        if (a.GetFrobeniusNorm() < 10 || a.GetOctahedralNorm() < 10 || a.GetCubicNorm() < 10) {
            continue;
        }
        bool check = true;
        std::set<Fraction> used;
        for (int i = 0; i < a.GetNumRows() && check; ++i) {
            for (int j = 0; j <= i; ++j) {
                if (used.count(a.At(i, j)) || fabs(a.At(i, j)) <= 1) {
                    check = false;
                    break;
                }
                used.insert(a.At(i, j));
            }
        }

        if (!check) {
            continue;
        }

        break;
    }
    if (iterations == 0) {
        fout << "Nothing" << std::endl;
        return;
    }

    Matrix l(a.GetNumRows(), a.GetNumColumns());
    Matrix d = l;

    Options::writer << "\\[" << std::endl;
    Options::writer << "A = " << a << std::endl;
    Options::writer << "\\]" << std::endl;

    Options::writer << "\\[" << std::endl;
    Options::writer << R"(d_{i\, i} = a_{i\, i} - \sum\limits_{j = 0}^{i - 1} d_{j\, j})";
    Options::writer << "\\cdot l_{i\\, j}^{2}" << std::endl;
    Options::writer << "\\]" << std::endl;
    Options::writer << "\\[" << std::endl;
    Options::writer << R"(l_{i\, j} = \frac{a_{i\, j} - \sum\limits_{k = 0}^{j - 1} d_{k\, k})";
    Options::writer << R"(\cdot l_{i\, k} \cdot l_{j\, k}}{d_{j\, j}})" << std::endl;
    Options::writer << "\\]" << std::endl;

    for(int column = 0; column < a.GetNumColumns(); ++column) {
        l.At(column, column) = 1;
        Options::writer << "\\[" << std::endl;
        Options::writer << "d_{" << column << "\\, " << column << "} = " << a.At(column, column);
        d.At(column, column) += a.At(column, column);
        for (int row = 0; row < column; ++row) {
           Options::writer << " - \\left(" << d.At(row, row) << R"(\right) \cdot \left()";
           Options::writer << l.At(column, row) << "\\right)^{2}";
           d.At(column, column) -= d.At(row, row) * l.At(column, row) * l.At(column, row);
        }
        Options::writer << " = " << d.At(column, column) << std::endl;
        Options::writer << "\\]" << std::endl;
        for (int row = column + 1; row < a.GetNumRows(); ++row) {
            Options::writer << "\\[" << std::endl;
            Options::writer << "l_{" << row << "\\, " << column << "} = \\frac{" << a.At(row, column);
            l.At(row, column) += a.At(row, column);
            for (int k = 0; k < column; k++) {
                Options::writer << " - \\left(" << d.At(k, k) << R"(\right) \cdot \left()";
                Options::writer << l.At(row, k) << R"(\right) \cdot \left()" << l.At(column, k) << "\\right)";
                l.At(row, column) -= d.At(k, k) * l.At(row, k) * l.At(column, k);
            }
            Options::writer << "}{" << d.At(column, column) << "}";
            l.At(row, column) /= d.At(column, column);
            Options::writer << " = " << l.At(row, column) << std::endl;
            Options::writer << "\\]" << std::endl;
        }
    }

    Options::writer << "\\[" << std::endl;
    Options::writer << "A = " << a << std::endl;
    Options::writer << "\\ =\\ " << std::endl;
    Options::writer << l << R"(\ \cdot \)" << std::endl;
    Options::writer << d << R"(\ \cdot \)" << std::endl;
    Options::writer << l.GetTranspose() << std::endl;
    Options::writer << "\\]" << std::endl;

    const Matrix x_0 = Matrix({{1}, {0}, {-1}, {0}});
    const Matrix x_1 = Matrix({{-2}, {1}, {0}, {-1}});

    SolveLDLT(l, d, a * x_0);
    SolveLDLT(l, d, a * x_1);

    Matrix tri_diagonal;
    iterations = max_iterations;
    while (--iterations) {
        std::vector<int> x(5), y(5), magic(4);
        for (int& item : x) {
            item = dis(twister);
        }
        for (int& item : y) {
            item = dis(twister);
        }
        for (int& item : magic) {
            item = dis(twister);
        }

        tri_diagonal = Matrix({{x[0], y[0], 0, 0, 0},
                  {magic[0] * x[0], magic[0] * y[0] + x[1], y[1], 0, 0},
                  {0, magic[1] * x[1], magic[1] * y[1] + x[2] + magic[2] * x[3], magic[2] * y[3], 0},
                  {0, 0, x[3], magic[3] * x[4] + y[3], magic[3] * y[4]},
                  {0, 0, 0, x[4], y[4]}});

        bool check = true;
        std::set<Fraction> used;
        for (int i = 0; i < tri_diagonal.GetNumRows() && check; ++i) {
            for (int j = std::max(0, i - 1); j < std::min(i + 2, tri_diagonal.GetNumColumns()); ++j) {
                if (used.count(tri_diagonal.At(i, j))) {
                    check = false;
                    break;
                }
                used.insert(tri_diagonal.At(i, j));
            }
        }

        if (!check) {
            continue;
        }

        break;
    }
    if (iterations == 0) {
        Options::writer << "Nothing" << std::endl;
        return;
    }

    const Matrix x_2({{1}, {-1}, {1}, {2}, {1}});

    ExtendedMatrix system(tri_diagonal, tri_diagonal * x_2);

    Options::writer << tri_diagonal << std::endl;
    Options::writer << system << std::endl;

    fout << Options::writer.str() << std::endl;
}


int main() {
  Options::step_by_step_type = Options::StepByStepType::MainSteps;
  Options::output_type = Options::OutputType::LaTex;

  HW4();

  return 0;
}
