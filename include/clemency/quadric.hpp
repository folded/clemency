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

#include <clemency/geom.hpp>

struct quadric_t {
  double coef[10];

  quadric_t() {}

  quadric_t(v3d_t n, v3d_t p, double d_off, double w = 1.0) {
    n = n.normalize();
    double d = -v3d_t::dot(n, p) + d_off;
    coef[0] = n.x * n.x * w;
    coef[1] = n.x * n.y * w;
    coef[2] = n.y * n.y * w;
    coef[3] = n.z * n.x * w;
    coef[4] = n.z * n.y * w;
    coef[5] = n.z * n.z * w;
    coef[6] = d * n.x * w;
    coef[7] = d * n.y * w;
    coef[8] = d * n.z * w;
    coef[9] = d * d * w;
  }

  m3d_t A() const {
    return m3d_t::init(coef[0], coef[1], coef[3], coef[1], coef[2], coef[4],
                       coef[3], coef[4], coef[5]);
  }

  v3d_t b() const { return v3d_t::init(coef[6], coef[7], coef[8]); }

  double c() const { return coef[9]; }

  bool minimize(const aabb3d_t& bbox, v3d_t& v) const {
    // TODO: implement
    return false;
  }

  bool minimize(v3d_t& v) const {
    m3d_t Ainv;
    double det;
    if (!A().invert(Ainv, det)) {
      return false;
    }

    v = Ainv.transform(b());

    return true;
  }

  double eval(const v3d_t& v) const {
    return coef[0] * v.x * v.x + 2 * coef[1] * v.y * v.x + coef[2] * v.y * v.y +
           2 * coef[3] * v.z * v.x + 2 * coef[4] * v.z * v.y +
           coef[5] * v.z * v.z + 2 * coef[6] * v.x + 2 * coef[7] * v.y +
           2 * coef[8] * v.z + coef[9];
  }

  quadric_t& operator+(const quadric_t& a) {
    for (size_t i = 0; i < 10; ++i)
      coef[i] += a.coef[i];
    return *this;
  }

  quadric_t& operator*=(double w) {
    for (size_t i = 0; i < 10; ++i)
      coef[i] *= w;
    return *this;
  }
};

static inline quadric_t operator*(double w, const quadric_t& a) {
  quadric_t r;
  for (size_t i = 0; i < 10; ++i)
    r.coef[i] = w * a.coef[i];
  return r;
}

static inline quadric_t operator+(const quadric_t& a, const quadric_t& b) {
  quadric_t r;
  for (size_t i = 0; i < 10; ++i)
    r.coef[i] = a.coef[i] + b.coef[i];
  return r;
}
