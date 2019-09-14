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

  explicit Fraction(int numerator, int denominator);

  [[nodiscard]] double ToDouble() const;

  void PrintFrac() const;

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

  static void Test();

 private:
  void Normalization();

  int numerator_, denominator_;
};

#endif //MATH__FRACTION_H_
