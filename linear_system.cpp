#include "matrix.h"

LinearSystem::SwapMatrix::SwapMatrix() :
    first_row_(0), second_row_(0),
    first_main_column_(0), second_main_column_(0),
    first_extension_column_(0), second_extension_column_(0) {}

void LinearSystem::SwapMatrix::SetRowSwap(int first_row, int second_row) {
  first_row_ = first_row + 1;
  second_row_ = second_row + 1;
}

void LinearSystem::SwapMatrix::SetMainColumnSwap(int first_column, int second_column) {
  first_main_column_ = first_column + 1;
  second_main_column_ = second_column + 1;
}

void LinearSystem::SwapMatrix::SetExtensionColumnSwap(int first_column, int second_column) {
  first_extension_column_ = first_column + 1;
  second_extension_column_ = second_column + 1;
}

bool LinearSystem::SwapMatrix::IsEmpty() const {
  return first_row_ == second_row_ &&
         first_main_column_ == second_main_column_ &&
         first_extension_column_ == second_extension_column_;
}

void LinearSystem::SwapMatrix::Clear() {
  first_row_ = second_row_ = 0;
  first_main_column_ = second_main_column_ = 0;
  first_extension_column_ = second_extension_column_ = 0;
}

std::string LinearSystem::SwapMatrix::ToString() const {
  if (IsEmpty()) {
    return "";
  }
  std::ostringstream out;
  bool started = false;
  if (first_row_ != second_row_) {
    started = true;
    out << "Swap rows: " << first_row_  << " " << second_row_;
  }
  if (first_main_column_ != second_main_column_) {
    if (started) {
      out << std::endl;
    }
    started = true;
    out << "Swap column: " << first_main_column_ << " " << second_main_column_;
  }
  if (first_extension_column_ != second_extension_column_) {
    if (started) {
      out << std::endl;
    }
    out << "Swap column: " << first_extension_column_ << " " << second_extension_column_;
  }

  return out.str();
}

std::string LinearSystem::SwapMatrix::ToLaTex() const {
  if (IsEmpty()) {
    return "";
  }
  std::ostringstream out;
  bool started = false;
  out << "\\begin{bmatrix}" << std::endl;
  if (first_row_ != second_row_) {
    started = true;
    out << "row: \\\\" << std::endl;
    out << first_row_  << " \\leftrightarrows " << second_row_;
  }
  bool is_column_writed = false;
  if (first_main_column_ != second_main_column_) {
    if (started) {
      out << " \\\\" << std::endl;
    }
    out << "col: \\\\" << std::endl;
    is_column_writed = true;
    started = true;
    out << first_main_column_ << " \\leftrightarrows " << second_main_column_;
  }
  if (first_extension_column_ != second_extension_column_) {
    if (started) {
      out << " \\\\" << std::endl;
    }
    if (!is_column_writed) {
        out << "col: \\\\" << std::endl;
    }
    out << first_extension_column_ << " \\leftrightarrows " << second_extension_column_;
  }
  out << std::endl;
  out << "\\end{bmatrix}";

  return out.str();
}

std::ostream& operator<<(std::ostream& out, const LinearSystem::SwapMatrix& swaps) {
  if (Options::output_type == Options::OutputType::Standard) {
    out << swaps.ToString();
  } else if (Options::output_type == Options::OutputType::LaTex){
    out << swaps.ToLaTex();
  }

  return out;
}

LinearSystem::LinearSystem(const ExtendedMatrix& system, Mode mode)
    : mode_(mode),
      n_(system.GetMainMatrix().GetNumRows()),
      m_(system.GetMainMatrix().GetNumColumns()),
      rank_(-1),
      start_system_(system),
      system_(system),
      history_() {}

LinearSystem::LinearSystem(const Matrix& a, const Matrix& b, Mode mode)
    : LinearSystem(ExtendedMatrix(a, b), mode) {}

bool LinearSystem::RunDirectGauss(Options::ChoiceType choice_type) {
  rank_ = 0;

  int max_pos = std::min(n_, m_);
  for (int pos = 0; pos < max_pos; ++pos) {
    if (!Choice(pos, choice_type) && !Choice(pos)) {
      continue;
    }

    ++rank_;

    PerformOperation(OperationType::MultiplyRow, pos, 0, 1 / system_.At(pos, pos));
    AddMainOutput();

    for (int row = pos + 1; row < n_; ++row) {
      Fraction coefficient = system_.At(row, pos) / system_.At(pos, pos);
      PerformOperation(OperationType::AddRow, row, pos, -coefficient);
    }

    AddMainOutput();
  }

  return rank_ = n_;
}

void LinearSystem::RunReverseGauss() {
  int max_pos = std::min(n_, m_);
  for (int pos = max_pos - 1; pos >= 0; --pos) {
    if (system_.At(pos, pos) == 0) {
      continue;
    }

    for (int row = pos - 1; row >= 0; --row) {
      Fraction coefficient = system_.At(row, pos);
      PerformOperation(OperationType::AddRow, row, pos, -coefficient);
    }

    AddMainOutput();
  }
}

bool LinearSystem::RunGauss(Options::ChoiceType choice_type) {
  RunDirectGauss(choice_type);
  RunReverseGauss();

  return rank_ == n_;
}

Matrix LinearSystem::GetSolutionMatrix() const {
  return system_.GetMainColumnPermutation().GetInverse() * system_.GetExtensionMatrix();
}

const ExtendedMatrix& LinearSystem::GetExtendedMatrix() const {
  return system_;
}

LUPMatrix LinearSystem::GetLUPMatrix() const {
  return LUPMatrix(system_);
}

LinearSystem::Operation::Operation(LinearSystem::OperationType operation_type,
                                   int lhs, int rhs, Fraction coefficient)
    : operation_type(operation_type), lhs(lhs), rhs(rhs), coefficient(coefficient) {}

bool LinearSystem::AddMainOutput() {
  if (!history_.empty() && history_.back().operation_type != OperationType::MainOutput) {
    history_.emplace_back();
    return true;
  }

  return false;
}

bool LinearSystem::PerformOperation(LinearSystem::OperationType operation_type, int lhs,
                                    int rhs, Fraction coefficient) {
  if (operation_type == OperationType::MainOutput) {
    return AddMainOutput();
  }

  Operation operation(operation_type, lhs, rhs, coefficient);
  if (ApplyOperation(operation)) {
    history_.push_back(operation);
    return true;
  }

  return false;
}

bool LinearSystem::ApplyOperation(OperationType operation_type, int lhs, int rhs,
                                  Fraction coefficient) {
  if (operation_type == OperationType::AddRow) {
    return system_.AddRow(lhs, rhs, coefficient);
  } else if (operation_type == OperationType::MultiplyRow) {
    return system_.MultiplyRow(lhs, coefficient);
  } else if (operation_type == OperationType::SwapRows) {
    bool result = false;
    result |= system_.SwapRows(lhs, rhs);
    if (mode_ == Mode::LU){
      int m = system_.GetMainMatrix().GetNumColumns();
      result |= system_.SwapColumns(lhs + m, rhs + m);
    }
    return result;
  } else if (operation_type == OperationType::SwapColumns) {
    return system_.SwapColumns(lhs, rhs);
  }

  return false;
}

bool LinearSystem::ApplyOperation(const Operation& operation) {
  return ApplyOperation(operation.operation_type, operation.lhs, operation.rhs,
      operation.coefficient);
}

bool LinearSystem::Choice(int pos, Options::ChoiceType choice_type) {
  int max_row = n_, max_column = m_;
  if (choice_type == Options::ChoiceType::Without) {
    return system_.At(pos, pos) != 0;
  } else if (choice_type == Options::ChoiceType::Row) {
    max_row = pos + 1;
  } else if (choice_type == Options::ChoiceType::Column) {
    max_column = pos + 1;
  }

  int best_column = pos, best_row = pos;
  Fraction best_item = system_.At(pos, pos);
  for (int row = pos; row < max_row; ++row) {
    for (int column = pos; column < max_column; ++column) {
      if (best_item < fabs(system_.At(row, column))) {
        best_row = row;
        best_column = column;
        best_item = fabs(system_.At(row, column));
      }
    }
  }

  PerformOperation(OperationType::SwapRows, pos, best_row, 0);
  PerformOperation(OperationType::SwapColumns, pos, best_column, 0);

  AddMainOutput();

  return best_item != 0;
}

std::ostream& operator<<(std::ostream& out, const LinearSystem& system) {
  if (Options::output_type == Options::OutputType::Standard) {
    out << system.ToString();
  } else if (Options::output_type == Options::OutputType::LaTex) {
    out << system.ToLaTex();
  }

  return out;
}

std::string LinearSystem::ToString() const {
  std::ostringstream out;
  LinearSystem current_system(start_system_, mode_);
  out << current_system.system_;
  SwapMatrix swaps;
  for (const auto& operation : history_) {
    current_system.ApplyOperation(operation);

    if (operation.operation_type == OperationType::SwapRows) {
      swaps.SetRowSwap(operation.lhs, operation.rhs);
      if (mode_ == Mode::LU) {
        int m = system_.GetMainMatrix().GetNumColumns();
        swaps.SetExtensionColumnSwap(operation.lhs + m, operation.rhs + m);
      }
    } else if (operation.operation_type == OperationType::SwapColumns) {
      swaps.SetMainColumnSwap(operation.lhs, operation.rhs);
    }

    if (Options::step_by_step_type != Options::StepByStepType::AllSteps &&
        (Options::step_by_step_type != Options::StepByStepType::MainSteps ||
            operation.operation_type != OperationType::MainOutput)) {
      continue;
    }

    out << std::endl;

    if (!swaps.IsEmpty()) {
      out << swaps << std::endl;
      swaps.Clear();
    }

    out << current_system.system_;
  }

  return out.str();
}

std::string LinearSystem::ToLaTex() const {
  std::ostringstream out;
  LinearSystem current_system(start_system_, mode_);

  out << "\\[" << std::endl;
  out << current_system.system_ << std::endl;

  SwapMatrix swaps;
  int matrices_count_in_row = 1;
  for (const auto& operation : history_) {
    current_system.ApplyOperation(operation);

    if (operation.operation_type == OperationType::SwapRows) {
      swaps.SetRowSwap(operation.lhs, operation.rhs);
      if (mode_ == Mode::LU) {
        int m = system_.GetMainMatrix().GetNumColumns();
        swaps.SetExtensionColumnSwap(operation.lhs + m, operation.rhs + m);
      }
    } else if (operation.operation_type == OperationType::SwapColumns) {
      swaps.SetMainColumnSwap(operation.lhs, operation.rhs);
    }

    if (Options::step_by_step_type != Options::StepByStepType::AllSteps &&
        (Options::step_by_step_type != Options::StepByStepType::MainSteps ||
         operation.operation_type != OperationType::MainOutput)) {
      continue;
    }

    if (matrices_count_in_row == Options::max_matrices_count_in_row) {
      out << "\\sim" << std::endl;
      out << "\\]" << std::endl;
      out << "\\[" << std::endl;
      matrices_count_in_row = 0;
    }

    if (!swaps.IsEmpty()) {
      out << "\\sim" << swaps;
      swaps.Clear();
    }

    out << "\\sim" << current_system.system_ << std::endl;
    ++matrices_count_in_row;
  }

  out << "\\]";

  return out.str();
}
