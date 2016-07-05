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



#include <clemency/geom/quat.hpp>
#include <cstdlib>
#include <cmath>



template<typename num_t>
class m3_t {
public:
  union {
    // .mCR
    struct {
      num_t
        m00, m01, m02,
        m10, m11, m12,
        m20, m21, m22;
    };
    // m[C][R]
    num_t m[3][3];
    // v[C*3+R]
    num_t v[9];
  };

  static m3_t ident() {
    m3_t r;
    r.m00 = num_t(1); r.m01 = num_t(0); r.m02 = num_t(0);
    r.m10 = num_t(0); r.m11 = num_t(1); r.m12 = num_t(0);
    r.m20 = num_t(0); r.m21 = num_t(0); r.m22 = num_t(1);
    return r;
  }

  static m3_t init(num_t m00, num_t m10, num_t m20,
                   num_t m01, num_t m11, num_t m21,
                   num_t m02, num_t m12, num_t m22) {
    m3_t r;
    r.m00 = m00; r.m01 = m01; r.m02 = m02;
    r.m10 = m10; r.m11 = m11; r.m12 = m12;
    r.m20 = m20; r.m21 = m21; r.m22 = m22;
    return r;
  }

  static m3_t rot(const quat_t<num_t> &q) {
    num_t w = q.a;
    num_t x = q.b.x;
    num_t y = q.b.y;
    num_t z = q.b.z;

    return m3_t::init(
      num_t(1) - num_t(2)*y*y - num_t(2)*z*z,            num_t(2)*x*y - num_t(2)*z*w,            num_t(2)*x*z + num_t(2)*y*w,
                 num_t(2)*x*y + num_t(2)*z*w, num_t(1) - num_t(2)*x*x - num_t(2)*z*z,            num_t(2)*y*z - num_t(2)*x*w,
                 num_t(2)*x*z - num_t(2)*y*w,            num_t(2)*y*z + num_t(2)*x*w, num_t(1) - num_t(2)*x*x - num_t(2)*y*y);
  }

  static m3_t scale(const v3_t<num_t> &v) {
    m3_t r;
    r.m00 =      v.x; r.m01 = num_t(0); r.m02 = num_t(0);
    r.m10 = num_t(0);      r.m11 = v.y; r.m12 = num_t(0);
    r.m20 = num_t(0); r.m21 = num_t(0);      r.m22 = v.z;
    return r;
  }

  m3_t transpose() const {
    m3_t r;
    r.m00 = m00; r.m01 = m10; r.m02 = m20;
    r.m10 = m01; r.m11 = m11; r.m12 = m21;
    r.m20 = m02; r.m21 = m12; r.m22 = m22;
    return r;
  }

  m3_t neg() const {
    m3_t r;
    for (size_t i = 0; i < 9; ++i) r.v[i] = -v[i];
    return r;
  }

  num_t operator()(int i, int j) const {
    return m[i][j];
  }

  num_t &operator()(int i, int j) {
    return m[i][j];
  }

  num_t mdet(int i, int j) const {
    return
      + m[i==0?1:0][j==0?1:0] * m[i==2?1:2][j==2?1:2]
      - m[i==2?1:2][j==0?1:0] * m[i==0?1:0][j==2?1:2];
  }

  num_t determinant() const {
    return
      +m00 * (m11 * m22 - m21 * m12)
      -m10 * (m01 * m22 - m21 * m02)
      +m20 * (m01 * m12 - m11 * m02);
  }

  bool invert(m3_t &inv, num_t &det) const {
    det = determinant();
    if (std::abs(det) < 1e-10) return false;
    num_t d = num_t(1)/det;

    inv.m00 = +d*mdet(0,0); inv.m01 = -d*mdet(1,0); inv.m02 = +d*mdet(2,0);
    inv.m10 = -d*mdet(0,1); inv.m11 = +d*mdet(1,1); inv.m12 = -d*mdet(2,1);
    inv.m20 = +d*mdet(0,2); inv.m21 = -d*mdet(1,2); inv.m22 = +d*mdet(2,2);

    return true;
  }

  v3_t<num_t> transform(const v3_t<num_t> &v) const {
    return v3_t<num_t>::init(v.x * m00 + v.y * m10 + v.z * m20,
                             v.x * m01 + v.y * m11 + v.z * m21,
                             v.x * m02 + v.y * m12 + v.z * m22);
  }

  bool exactly_eq(const m3_t &o) const {
    for (size_t i = 0; i < 9; ++i) if (v[i] != o.v[i]) return false;
    return true;
  }

  bool eq(const m3_t &o, num_t eps = 1e-10) const {
    for (size_t i = 0; i < 9; ++i) if (std::abs(v[i] - o.v[i]) > eps) return false;
    return true;
  }

  static m3_t sub(const m3_t &a, const m3_t &b) {
    m3_t r;
    for (size_t i = 0; i < 9; ++i) r.v[i] = a.v[i] - b.v[i];
    return r;
  }

  static m3_t add(const m3_t &a, const m3_t &b) {
    m3_t r;
    for (size_t i = 0; i < 9; ++i) r.v[i] = a.v[i] + b.v[i];
    return r;
  }

  static m3_t mul(const m3_t &a, const m3_t &b) {
    m3_t r;
    for (size_t i = 0; i < 3; ++i) {
      for (size_t j = 0; j < 3; ++j) {
        r.m[i][j] =
          a.m[i][0] * b.m[0][j] +
          a.m[i][1] * b.m[1][j] +
          a.m[i][2] * b.m[2][j];
      }
    }
    return r;
  }
};

template<typename num_t>
bool operator==(const m3_t<num_t> &a, const m3_t<num_t> &b) {
  return a.exactly_eq(b);
}

template<typename num_t>
m3_t<num_t> operator-(const m3_t<num_t> &a, const m3_t<num_t> &b) {
  return m3_t<num_t>::sub(a, b);
}

template<typename num_t>
m3_t<num_t> operator+(const m3_t<num_t> &a, const m3_t<num_t> &b) {
  return m3_t<num_t>::add(a, b);
}

template<typename num_t>
m3_t<num_t> operator-(const m3_t<num_t> &a) {
  return a.neg();
}

template<typename num_t>
m3_t<num_t> operator*(const m3_t<num_t> &a, const m3_t<num_t> &b) {
  return m3_t<num_t>::mul(a, b);
}

template<typename num_t>
std::ostream &operator<<(std::ostream &o, const m3_t<num_t> &m) {
  for (size_t i = 0; i < 3; ++i) {
    o << "[";
    for (size_t j = 0; j < 3; ++j) {
      o << " " << m.m[j][i];
    }
    o << " ]\n";
  }
  return o;
}



typedef m3_t<double> m3d_t;
typedef m3_t<float>  m3f_t;
