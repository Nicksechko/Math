#include "fraction.h"

Fraction::Fraction() : numerator_(0), denominator_(1) {}

Fraction::Fraction(int numerator, int denominator) {
  assert(denominator != 0);
  numerator_ = numerator;
  denominator_ = denominator;
  Normalization();
}

double Fraction::ToDouble() const {
  return (double) numerator_ / denominator_;
}

std::string Fraction::ToLaTex() const {
  return "\frac{" + std::to_string(numerator_) + "}" +
      "{" + std::to_string(denominator_);
}

void Fraction::PrintFrac() const {
  cout << numerator_ << "/" << denominator_ << endl;
}

Fraction operator+(const Fraction& a, const Fraction& b) {
  Fraction current_frac = a;
  current_frac.denominator_ = lcm(a.denominator_, b.denominator_);
  current_frac.numerator_ *=
      current_frac.denominator_ / a.denominator_;
  current_frac.numerator_ +=
      b.numerator_ * (current_frac.denominator_ / b.denominator_);
  current_frac.Normalization();
  return current_frac;
}

Fraction operator-(const Fraction& a, const Fraction& b) {
  Fraction current_frac = a;
  current_frac.denominator_ = lcm(a.denominator_, b.denominator_);
  current_frac.numerator_ *=
      current_frac.denominator_ / a.denominator_;
  current_frac.numerator_ -=
      b.numerator_ * (current_frac.denominator_ / b.denominator_);
  current_frac.Normalization();
  return current_frac;
}

Fraction operator*(const Fraction& a, const Fraction& b) {
  Fraction current_frac = a;
  current_frac.numerator_ *= b.numerator_;
  current_frac.denominator_ *= b.denominator_;
  current_frac.Normalization();
  return current_frac;
}

Fraction operator/(const Fraction& a, const Fraction& b) {
  assert(b.numerator_ != 0);
  Fraction current_frac = a;
  current_frac.numerator_ *= b.denominator_;
  current_frac.denominator_ *= b.numerator_;
  current_frac.Normalization();
  return current_frac;
}

bool operator==(const Fraction& a, const Fraction& b) {
  return a.numerator_ == b.numerator_ && a.denominator_ == b.denominator_;
}

bool operator!=(const Fraction& a, const Fraction& b) {
  return a.numerator_ != b.numerator_ || a.denominator_ != b.denominator_;
}

bool operator<(const Fraction& a, const Fraction& b) {
  return a.numerator_ * b.denominator_ < a.denominator_ * b.numerator_;
}

bool operator>(const Fraction& a, const Fraction& b) {
  return a.numerator_ * b.denominator_ > a.denominator_ * b.numerator_;
}

bool operator>=(const Fraction& a, const Fraction& b) {
  return a.numerator_ * b.denominator_ >= a.denominator_ * b.numerator_;
}

bool operator<=(const Fraction& a, const Fraction& b) {
  return a.numerator_ * b.denominator_ <= a.denominator_ * b.numerator_;
}

Fraction operator+(const Fraction& a, const int& b) {
  Fraction current_frac = a;
  current_frac.numerator_ += a.denominator_ * b;
  current_frac.Normalization();
  return current_frac;
}

Fraction operator-(const Fraction& a, const int& b) {
  Fraction current_frac = a;
  current_frac.numerator_ -= a.denominator_ * b;
  current_frac.Normalization();
  return current_frac;
}

Fraction operator*(const Fraction& a, const int& b) {
  Fraction current_frac = a;
  current_frac.numerator_ *= b;
  current_frac.Normalization();
  return current_frac;
}

Fraction operator/(const Fraction& a, const int& b) {
  assert(b != 0);
  Fraction current_frac = a;
  current_frac.denominator_ *= b;
  current_frac.Normalization();
  return current_frac;
}

Fraction operator+(const int& a, const Fraction& b) {
  Fraction current_frac = b;
  current_frac.numerator_ += b.denominator_ * a;
  current_frac.Normalization();
  return current_frac;
}

Fraction operator-(const int& a, const Fraction& b) {
  Fraction current_frac = b;
  current_frac.numerator_ -= b.denominator_ * a;
  current_frac.numerator_ = -current_frac.numerator_;
  current_frac.Normalization();
  return current_frac;
}

Fraction operator*(const int& a, const Fraction& b) {
  Fraction current_frac = b;
  current_frac.numerator_ *= a;
  current_frac.Normalization();
  return current_frac;
}

Fraction operator/(const int& a, const Fraction& b) {
  assert(b.numerator_ != 0);
  Fraction current_frac(a, 1);
  current_frac.numerator_ *= b.denominator_;
  current_frac.denominator_ *= b.numerator_;
  current_frac.Normalization();
  return current_frac;
}

const Fraction& Fraction::operator+=(const Fraction& b) {
  int current_lcm = lcm(denominator_, b.denominator_);
  numerator_ *= current_lcm / denominator_;
  numerator_ += b.numerator_ * (current_lcm / b.denominator_);
  denominator_ = current_lcm;
  Normalization();
  return *this;
}

const Fraction& Fraction::operator-=(const Fraction& b) {
  int current_lcm = lcm(denominator_, b.denominator_);
  numerator_ *= current_lcm / denominator_;
  numerator_ -= b.numerator_ * (current_lcm / b.denominator_);
  denominator_ = current_lcm;
  Normalization();
  return *this;
}

const Fraction& Fraction::operator*=(const Fraction& b) {
  numerator_ *= b.numerator_;
  denominator_ *= b.denominator_;
  Normalization();
  return *this;
}

const Fraction& Fraction::operator/=(const Fraction& b) {
  assert(b.numerator_ != 0);
  numerator_ *= b.denominator_;
  denominator_ *= b.numerator_;
  Normalization();
  return *this;
}

const Fraction& Fraction::operator+=(const int& b) {
  numerator_ += denominator_ * b;
  Normalization();
  return *this;
}

const Fraction& Fraction::operator-=(const int& b) {
  numerator_ -= denominator_ * b;
  Normalization();
  return *this;
}

const Fraction& Fraction::operator*=(const int& b) {
  numerator_ *= b;
  Normalization();
  return *this;
}

const Fraction& Fraction::operator/=(const int& b) {
  assert(b != 0);
  denominator_ *= b;
  Normalization();
  return *this;
}

void Fraction::Test() {
  {//Test Case 1
    const double eps = 1e-9;
    assert(fabs(Fraction(0, 1).ToDouble() - (double) 0) < eps);
    assert(fabs(Fraction(2, 3).ToDouble() - (double) 2 / 3) < eps);
    assert(fabs(Fraction(-7, 3).ToDouble() - (double) (-7) / 3) < eps);
    assert(fabs(Fraction(5, -9).ToDouble() - (double) (-5) / 9) < eps);
    assert(fabs(Fraction(-3, -11).ToDouble() - (double) 3 / 11) < eps);
  }
  {//Test Case 2
    Fraction(0, 1).PrintFrac();
    Fraction(2, 3).PrintFrac();
    Fraction(-7, 3).PrintFrac();
    Fraction(5, -9).PrintFrac();
    Fraction(-3, -11).PrintFrac();
  }
  //Test Case 3
  {
    //TODO
  }
  {//Test Case 4
    assert(Fraction(0, 1) + Fraction(1, 3) == Fraction(1, 3));
    assert(Fraction(-5, 2) + Fraction(6, 5) == Fraction(-13, 10));
    assert(Fraction(-5, 2) - Fraction(0, 1) == Fraction(-5, 2));
    assert(Fraction(3, 7) - Fraction(-6, 25) == Fraction(117, 175));
    assert(Fraction(-5, 2) * Fraction(0, 5) == Fraction(0, 1));
    assert(Fraction(-5, 2) * Fraction(-2, 5) == Fraction(1, 1));
    assert(Fraction(0, 1) / Fraction(2, 5) == Fraction(0, 1));
    assert(Fraction(-5, 2) / Fraction(2, 5) == Fraction(-25, 4));
  }
  {//Test Case 5
    assert(Fraction(1, 3) == Fraction(1, 3));
    assert(Fraction(4, 2) == Fraction(2, 1));
    assert(Fraction(-5, 2) != Fraction(6, 5));
    assert(Fraction(5, 2) != Fraction(11, 4));
    assert(Fraction(-5, 2) < Fraction(0, 1));
    assert(Fraction(3, 5) < Fraction(2, 3));
    assert(Fraction(3, 7) > Fraction(-6, 25));
    assert(Fraction(5, 2) > Fraction(0, 1));
    assert(Fraction(-1, 2) <= Fraction(0, 1));
    assert(Fraction(-1, 2) <= Fraction(-4, 8));
    assert(Fraction(-2, 5) >= Fraction(-5, 2));
    assert(Fraction(-10, 4) >= Fraction(-5, 2));
  }
  {//Test Case 6
    assert(Fraction(0, 1) + 1 == Fraction(1, 1));
    assert(Fraction(-5, 2) + 5 == Fraction(5, 2));
    assert(Fraction(-5, 2) - 0 == Fraction(-5, 2));
    assert(Fraction(3, 7) - 25 == Fraction(-172, 7));
    assert(Fraction(-5, 2) * 0 == Fraction(0, 1));
    assert(Fraction(-5, 2) * (-2) == Fraction(5, 1));
    assert(Fraction(0, 1) / 5 == Fraction(0, 1));
    assert(Fraction(-5, 2) / 2 == Fraction(-5, 4));
  }
  {//Test Case 7
    assert(3 + Fraction(0, 1) == Fraction(3, 1));
    assert(6 + Fraction(-5, 2) == Fraction(7, 2));
    assert(1 - Fraction(-5, 2) == Fraction(7, 2));
    assert(25 - Fraction(3, 7) == Fraction(172, 7));
    assert(5 * Fraction(-5, 2) == Fraction(-25, 2));
    assert((-5) * Fraction(-5, 2) == Fraction(25, 2));
    assert(2 / Fraction(1, 1) == Fraction(2, 1));
    assert(2 / Fraction(-5, 2) == Fraction(-4, 5));
  }
  {//Test Case 8
    Fraction a;
    a = Fraction(0, 1);
    a += Fraction(1, 3);
    assert(a == Fraction(1, 3));

    a = Fraction(-5, 2);
    a += Fraction(6, 5);
    assert(a == Fraction(-13, 10));

    a = Fraction(-5, 2);
    a -= Fraction(0, 1);
    assert(a == Fraction(-5, 2));

    a = Fraction(3, 7);
    a -= Fraction(-6, 25);
    assert(a == Fraction(117, 175));

    a = Fraction(-5, 2);
    a *= Fraction(0, 5);
    assert(a == Fraction(0, 1));

    a = Fraction(-5, 2);
    a *= Fraction(-2, 5);
    assert(a == Fraction(1, 1));

    a = Fraction(0, 1);
    a /= Fraction(2, 5);
    assert(a == Fraction(0, 1));

    a = Fraction(-5, 2);
    a /= Fraction(2, 5);
    assert(a == Fraction(-25, 4));
  }

  {//Test Case 9
    Fraction a = Fraction(0, 1);
    a += 1;
    assert(a == Fraction(1, 1));

    a = Fraction(-5, 2);
    a += 5;
    assert(a == Fraction(5, 2));

    a = Fraction(-5, 2);
    a -= 0;
    assert(a == Fraction(-5, 2));

    a = Fraction(3, 7);
    a -= 25;
    assert(a == Fraction(-172, 7));

    a = Fraction(-5, 2);
    a *= 0;
    assert(a == Fraction(0, 1));

    a = Fraction(-5, 2);
    a *= (-2);
    assert(a == Fraction(5, 1));

    a = Fraction(0, 1);
    a /= 5;
    assert(a == Fraction(0, 1));

    a = Fraction(-5, 2);
    a /= 2;
    assert(a == Fraction(-5, 4));
  }
}

void Fraction::Normalization() {
  if (denominator_ < 0) {
    numerator_ = -numerator_;
    denominator_ = -denominator_;
  }
  int g = gcd(numerator_, denominator_);
  numerator_ /= g;
  denominator_ /= g;
}