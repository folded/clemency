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

#include <clemency/config.h>
#include <clemency/geom.hpp>



double v3_t::dotcross(const v3_t &a, const v3_t &b, const v3_t &c) {
  return
    (a.x * b.y * c.z + a.y * b.z * c.x + a.z * b.x * c.y) -
    (a.x * c.y * b.z + a.y * c.z * b.x + a.z * c.x * b.y);
}



double v3_t::orient(const v3_t &a, const v3_t &b, const v3_t &c, const v3_t &d) {
  return dotcross((a - d), (b - d), (c - d));
}



double v3_t::tetrahedron_volume(const v3_t &a, const v3_t &b, const v3_t &c, const v3_t &d) {
  return dotcross((a - d), (b - d), (c - d)) / 6.0;
}



std::vector<double> solve3(double a, double b, double c, double d, const double zero) {
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
