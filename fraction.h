#ifndef MATH__FRACTION_H_
#define MATH__FRACTION_H_

#include <iostream>

#include <iostream>
#include <cassert>
#include <algorithm>
#include <numeric>
#include <cmath>

using namespace std;

class Fraction {
 public:
  Fraction();

  explicit Fraction(int numerator);

  explicit Fraction(int numerator, int denominator);

  [[nodiscard]] double ToDouble() const;

  [[nodiscard]] std::string ToLaTex() const;

  friend Fraction operator+(const Fraction& a, const Fraction& b);

  friend Fraction operator-(const Fraction& a, const Fraction& b);

  friend Fraction operator*(const Fraction& a, const Fraction& b);

  friend Fraction operator/(const Fraction& a, const Fraction& b);

  friend bool operator==(const Fraction& a, const Fraction& b);

  friend bool operator!=(const Fraction& a, const Fraction& b);

  friend bool operator<(const Fraction& a, const Fraction& b);

  friend bool operator>(const Fraction& a, const Fraction& b);

  friend bool operator>=(const Fraction& a, const Fraction& b);

  friend bool operator<=(const Fraction& a, const Fraction& b);

  friend Fraction operator+(const Fraction& a, const int& b);

  friend Fraction operator-(const Fraction& a, const int& b);

  friend Fraction operator*(const Fraction& a, const int& b);

  friend Fraction operator/(const Fraction& a, const int& b);

  friend Fraction operator+(const int& a, const Fraction& b);

  friend Fraction operator-(const int& a, const Fraction& b);

  friend Fraction operator*(const int& a, const Fraction& b);

  friend Fraction operator/(const int& a, const Fraction& b);

  const Fraction& operator+=(const Fraction& b);

  const Fraction& operator-=(const Fraction& b);

  const Fraction& operator*=(const Fraction& b);

  const Fraction& operator/=(const Fraction& b);

  const Fraction& operator+=(const int& b);

  const Fraction& operator-=(const int& b);

  const Fraction& operator*=(const int& b);

  const Fraction& operator/=(const int& b);

  friend istream& operator>>(istream& in, Fraction& f);

  friend ostream& operator<<(ostream& out, const Fraction& f);

  static void Test();

 private:
  void Normalization();

  int numerator_;
  int denominator_;
};

#endif //MATH__FRACTION_H_
