#include "matrix.h"

Permutation::Permutation(int n)
  : permutation_(n) {
  std::iota(permutation_.begin(), permutation_.end(), 0);
}

Permutation::Permutation(std::vector<int> permutation)
  : permutation_(std::move(permutation)) {
  std::vector<int> checker(permutation_.size());
  for (int item : permutation_) {
    assert(checker.size() <= item || item < 0 || checker[item] != 0);
    ++checker[item];
  }
}

bool Permutation::Transposition(int i, int j) {
  if (i == j) {
    return false;
  } else {
    std::swap(permutation_[i], permutation_[j]);
    return true;
  }
}

int Permutation::At(int pos) const {
  return permutation_[pos];
}

int Permutation::Size() const {
  return permutation_.size();
}

Matrix Permutation::AsMatrix() const {
  return Matrix(*this);
}

Permutation Permutation::GetInverse() const {
  Permutation result(permutation_.size());
  for (int i = 0; i < static_cast<int>(permutation_.size()); ++i) {
    result.permutation_[permutation_[i]] = i;
  }

  return result;
}

void Permutation::Inverse() {
  operator=(GetInverse());
}

Permutation operator*(const Permutation& lhs, const Permutation& rhs) {
  assert(lhs.Size() == rhs.Size());
  Permutation result(lhs.Size());
  for (int i = 0; i < lhs.Size(); ++i) {
    result.permutation_[i] = lhs.permutation_[rhs.permutation_[i]];
  }

  return result;
}

std::string Permutation::ToString() const {
  std::ostringstream out;
  for (int i = 0; i < static_cast<int>(permutation_.size()); ++i) {
    out << permutation_[i] + 1;
    if (i + 1 != permutation_.size()) {
      out << " ";
    }
  }

  return out.str();
}

std::string Permutation::ToLaTex() const {
  std::ostringstream out;
  for (int i = 0; i < static_cast<int>(permutation_.size()); ++i) {
    out << permutation_[i] + 1;
    if (i + 1 != permutation_.size()) {
      out << "\\ ";
    }
  }

  return out.str();
}

std::ostream& operator<<(std::ostream& out, const Permutation& permutation) {
  if (Options::output_type == Options::OutputType::Standard) {
    out << permutation.ToString();
  } else if (Options::output_type == Options::OutputType::LaTex) {
    out << permutation.ToLaTex();
  }

  return out;
}

bool Permutation::IsIdentity() const {
    for (int i = 0; i < static_cast<int>(permutation_.size()); ++i) {
        if (permutation_[i] != i) {
            return false;
        }
    }

    return true;
}

