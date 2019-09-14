#ifndef MATH__FRACTION_H_
#define MATH__FRACTION_H_

#include <iostream>

#include <iostream>
#include <cassert>
#include <algorithm>
#include <numeric>
#include <cmath>

class Fraction {
 public:
  Fraction();

  Fraction(int numerator);

  Fraction(int numerator, int denominator);

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

  friend std::istream& operator>>(std::istream& in, Fraction& f);

  friend std::ostream& operator<<(std::ostream& out, const Fraction& f);

  static void Test();

 private:
  void Normalization();

  int numerator_;
  int denominator_;
};

#endif //MATH__FRACTION_H_
