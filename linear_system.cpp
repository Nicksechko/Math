#include "matrix.h"

LinearSystem::LinearSystem(Matrix A, Matrix B)
    : first(true),
      block_flag(0),
      n_(A.GetNumRows()),
      rank_(-1),
      A_(std::move(A)),
      B_(std::move(B)),
      row_permutation_(A.GetNumRows()),
      column_permutation_(A.GetNumColumns()) {
  assert(A_.GetNumRows() == B_.GetNumRows());
  assert(A_.GetNumRows() == A_.GetNumColumns());
  std::iota(row_permutation_.begin(), row_permutation_.end(), 0);
  std::iota(column_permutation_.begin(), column_permutation_.end(), 0);
}

int LinearSystem::RunDirectGauss(Options::ChoiceType choice_type) {
  if (Options::step_by_step_type != Options::StepByStepType::Without) {
    std::cout << "Direct Gauss for system:" << std::endl << std::endl;
    OutputSystem();
  }
  for (int pos = 0; pos < n_; ++pos) {
    Choice(pos, choice_type);
    if (A_.At(pos, pos) == Fraction(0)) {
      Choice(pos);
    }
    if (A_.At(pos, pos) == Fraction(0)) {
      return rank_ = pos;
    }
    bool changed = false;
    for (int row = pos + 1; row < A_.GetNumRows(); ++row) {
      if (A_.At(row, pos) != 0) {
        changed = true;
        Fraction coefficient = A_.At(row, pos) / A_.At(pos, pos);
        A_.SubtractRows(row, pos, coefficient);
        B_.SubtractRows(row, pos, coefficient);
        if (Options::step_by_step_type == Options::StepByStepType::AllSteps) {
          OutputSystem();
        }
      }
    }
    if (Options::step_by_step_type == Options::StepByStepType::MainSteps && changed) {
      OutputSystem();
    }
  }
  if (Options::step_by_step_type != Options::StepByStepType::Without && block_flag != 0) {
    std::cout << "$$" << std::endl;
  }
  first = true;
  block_flag = 0;

  return rank_ = n_;
}

void LinearSystem::RunReverseGauss() {
  if (Options::step_by_step_type != Options::StepByStepType::Without) {
    std::cout << "Reverse Gauss for system:" << std::endl << std::endl;
    OutputSystem();
  }
  for (int pos = n_ - 1; pos >= 0; --pos) {
    if (A_.At(pos, pos) == Fraction(0)) {
      continue;
    }
    if (A_.At(pos, pos) != 1) {
      B_.DivideRow(pos, A_.At(pos, pos));
      A_.DivideRow(pos, A_.At(pos, pos));
      if (Options::step_by_step_type == Options::StepByStepType::AllSteps) {
        OutputSystem();
      }
    }
    bool changed = false;
    for (int row = pos - 1; row >= 0; --row) {
      if (A_.At(row, pos) != 0) {
        changed = true;
        Fraction coefficient = A_.At(row, pos);
        A_.SubtractRows(row, pos, coefficient);
        B_.SubtractRows(row, pos, coefficient);
        if (Options::step_by_step_type == Options::StepByStepType::AllSteps) {
          OutputSystem();
        }
      }
    }
    if (Options::step_by_step_type == Options::StepByStepType::MainSteps && changed) {
      OutputSystem();
    }
  }
  if (Options::step_by_step_type != Options::StepByStepType::Without && block_flag != 0) {
    std::cout << "$$" << std::endl;
  }
  first = true;
  block_flag = 0;
}

int LinearSystem::RunGauss(Options::ChoiceType choice_type) {
  if (Options::step_by_step_type != Options::StepByStepType::Without) {
    std::cout << "Gauss" << std::endl << std::endl;
  }
  RunDirectGauss(choice_type);
  RunReverseGauss();

  return rank_;
}

void LinearSystem::OutputSystem() {
  if (Options::output_type == Options::OutputType::Standard) {
    std::cout << n_ << ' ' << B_.GetNumColumns() << std::endl;
    for (int row = 0; row < n_; ++row) {
      for (int column = 0; column < n_; ++column) {
        std::cout << A_.At(row, column) << ' ';
      }
      std::cout << ' ';
      for (int column = 0; column < B_.GetNumColumns(); ++column) {
        std::cout << B_.At(row, column) << ' ';
      }
      std::cout << std::endl;
    }
  } else {
    if (block_flag == 0) {
      std::cout << "$$" << std::endl;
    }
    if (!first) {
      std::cout << "\\sim";
    }
    first = false;
    std::cout << "\\begin{bmatrix}" << std::endl;
    for (int row = 0; row < n_; ++row) {
      for (int column = 0; column < n_; ++column) {
        std::cout << A_.At(row, column).ToLaTex() << " & ";
      }
      std::cout << "\\vline & ";
      for (int column = 0; column < B_.GetNumColumns(); ++column) {
        std::cout << B_.At(row, column).ToLaTex();
        if (column + 1 != B_.GetNumColumns()) {
          std::cout << " & ";
        }
      }
      std::cout << "\\\\" << std::endl;
    }
    std::cout << "\\end{bmatrix}" << std::endl;
    if (block_flag + 1 == Options::latex_block_size) {
      std::cout << "\\sim" << std::endl;
      std::cout << "$$" << std::endl;
      block_flag = 0;
    } else {
      ++block_flag;
    }
  }
}

Matrix LinearSystem::GetSolutionMatrix() {
  Matrix solution_matrix(B_.GetNumRows(), B_.GetNumColumns());
  for (int row = 0; row < static_cast<int>(column_permutation_.size()); row++) {
    solution_matrix.CopyRow(B_, row, column_permutation_[row]);
  }

  return solution_matrix;
}

LUPMatrix LinearSystem::ToLUPMatrix() {
  return LUPMatrix(B_, A_, row_permutation_, column_permutation_);
}

void LinearSystem::Choice(int pos, Options::ChoiceType choice_type) {
  int max_row = n_, max_column = n_;
  if (choice_type == Options::ChoiceType::Without) {
    return;
  } else if (choice_type == Options::ChoiceType::Row) {
    max_row = pos + 1;
  } else if (choice_type == Options::ChoiceType::Column) {
    max_column = pos + 1;
  }

  int best_column = pos, best_row = pos;
  Fraction best_item = A_.At(pos, pos);
  for (int row = pos; row < max_row; ++row) {
    for (int column = pos; column < max_column; ++column) {
      if (best_item < fabs(A_.At(row, column))) {
        best_row = row;
        best_column = column;
        best_item = fabs(A_.At(row, column));
      }
    }
  }

  if (pos != best_row) {
    std::swap(row_permutation_[pos], row_permutation_[best_row]);
    A_.SwapRows(pos, best_row);
    B_.SwapRows(pos, best_row);
  }
  if (pos != best_column) {
    std::swap(column_permutation_[pos], column_permutation_[best_column]);
    A_.SwapColumns(pos, best_column);
  }
  if (Options::step_by_step_type != Options::StepByStepType::Without && (best_row != pos || best_column != pos)) {
    OutputSystem();
  }
}
