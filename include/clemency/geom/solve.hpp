// clemency - a OpenCL marching cubes mesh generator
//
// copyright (C) 2013, Tobias Sargeant <tobias.sargeant@gmail.com>
// All rights reserved.
//
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 2 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License along
// with this program; if not, write to the Free Software Foundation, Inc.,
// 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.

#pragma once



#include <math.h>
#include <vector>



static inline std::vector<double> solve1(double a, double b, const double zero = 1e-10) {
    std::vector<double> r;

    if (::fabs(a) > zero) {
      r.push_back(-b / a);
    }

    return r;
}

static inline std::vector<double> solve2(double a, double b, double c, const double zero = 1e-10) {
  std::vector<double> r;

  if (::fabs(a) < zero) {
    return solve1(b, c, zero);
  }

  double D = b*b - 4.0*a*c;
  if (D > 0) {
    if (D < zero) {
      r.push_back(-b / (2.0 * a));
    } else {
      r.reserve(2);

      D = ::sqrt(D);
      r.push_back((-b - D) / (2.0 * a));
      r.push_back((-b + D) / (2.0 * a));
    }
  }

  return r;
}

static inline std::vector<double> solve3(double a, double b, double c, double d, const double zero = 1e-10) {
  std::vector<double> r;

  if (::fabs(a) < zero) {
    return solve2(b, c, d, zero);
  }

  b /= a;
  c /= a;
  d /= a;

  double p = c - b*b / 3.0;
  double q = b * (2.0*b*b - 9.0*c) / 27.0 + d;
  double p3 = p*p*p;
  double D = q*q + 4.0*p3 / 27.0;
  double offset = -b / 3.0;

  if (D > zero) {
    double z = ::sqrt(D);
    double u = (-q + z) / 2.0;
    double v = (-q - z) / 2.0;
    u = (u >= 0) ? ::pow(u, 1.0/3.0) : -::pow(-u, 1.0/3.0);
    v = (v >= 0) ? ::pow(v, 1.0/3.0) : -::pow(-v, 1.0/3.0);

    r.push_back(u + v + offset);
  } else if (D < -zero) {
    double u = 2.0 * ::sqrt(-p / 3.0);
    double v = ::acos(-::sqrt(-27.0 / p3) * q / 2.0) / 3.0;

    r.reserve(3);

    r.push_back(u * ::cos(v) + offset);
    r.push_back(u * ::cos(v + 2.0 * M_PI / 3.0) + offset);
    r.push_back(u * ::cos(v + 4.0 * M_PI / 3.0) + offset);
  } else {
    double u = (q < 0) ? ::pow(-q / 2.0, 1.0 / 3.0) : -::pow(q / 2.0, 1.0 / 3.0);

    r.reserve(2);

    r.push_back(2.0*u + offset);
    r.push_back(-u + offset);
  }

  return r;
}
