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



#include "v3.hpp"
#include "quat.hpp"



template<typename num_t>
class m4_t {
public:
  union {
    // column major storage.
    // .mCR
    struct {
      num_t
        m00, m01, m02, m03,
        m10, m11, m12, m13,
        m20, m21, m22, m23,
        m30, m31, m32, m33;
    };
    // m[C][R]
    num_t m[4][4];
    // v[C*4+R]
    num_t v[16];
  };

  static m4_t ident() {
    m4_t r;
    r.m00 = num_t(1); r.m01 = num_t(0); r.m02 = num_t(0); r.m03 = num_t(0);
    r.m10 = num_t(0); r.m11 = num_t(1); r.m12 = num_t(0); r.m13 = num_t(0);
    r.m20 = num_t(0); r.m21 = num_t(0); r.m22 = num_t(1); r.m23 = num_t(0);
    r.m30 = num_t(0); r.m31 = num_t(0); r.m32 = num_t(0); r.m33 = num_t(1);
    return r;
  }

  // note: arguments are provided in the natural orientation, but stored transposed.
  static m4_t init(num_t m00, num_t m10, num_t m20, num_t m30,
                   num_t m01, num_t m11, num_t m21, num_t m31,
                   num_t m02, num_t m12, num_t m22, num_t m32,
                   num_t m03, num_t m13, num_t m23, num_t m33) {
    m4_t r;
    r.m00 = m00; r.m01 = m01; r.m02 = m02; r.m03 = m03;
    r.m10 = m10; r.m11 = m11; r.m12 = m12; r.m13 = m13;
    r.m20 = m20; r.m21 = m21; r.m22 = m22; r.m23 = m23;
    r.m30 = m30; r.m31 = m31; r.m32 = m32; r.m33 = m33;
    return r;
  }

  static m4_t rot(const quat_t<num_t> &q) {
    num_t w = q.a;
    num_t x = q.b.x;
    num_t y = q.b.y;
    num_t z = q.b.z;

    return m4_t::init(
      num_t(1) - num_t(2)*y*y - num_t(2)*z*z,            num_t(2)*x*y - num_t(2)*z*w,            num_t(2)*x*z + num_t(2)*y*w, num_t(0),  
                 num_t(2)*x*y + num_t(2)*z*w, num_t(1) - num_t(2)*x*x - num_t(2)*z*z,            num_t(2)*y*z - num_t(2)*x*w, num_t(0),  
                 num_t(2)*x*z - num_t(2)*y*w,            num_t(2)*y*z + num_t(2)*x*w, num_t(1) - num_t(2)*x*x - num_t(2)*y*y, num_t(0),  
                                    num_t(0),                               num_t(0),                               num_t(0), num_t(1));
  }

  static m4_t scale(const v3_t<num_t> &v) {
    m4_t r;
    r.m00 =      v.x; r.m01 = num_t(0); r.m02 = num_t(0); r.m03 = num_t(0);
    r.m10 = num_t(0); r.m11 =      v.y; r.m12 = num_t(0); r.m13 = num_t(0);
    r.m20 = num_t(0); r.m21 = num_t(0);      r.m22 = v.z; r.m23 = num_t(0);
    r.m30 = num_t(0); r.m31 = num_t(0); r.m32 = num_t(0); r.m33 = num_t(1);
    return r;
  }

  static m4_t translate(const v3_t<num_t> &v) {
    m4_t r;
    r.m00 = num_t(1); r.m01 = num_t(0); r.m02 = num_t(0); r.m03 = num_t(0);
    r.m10 = num_t(0); r.m11 = num_t(1); r.m12 = num_t(0); r.m13 = num_t(0);
    r.m20 = num_t(0); r.m21 = num_t(0); r.m22 = num_t(1); r.m23 = num_t(0);
    r.m30 = v.x;      r.m31 = v.y;      r.m32 = v.z;      r.m33 = num_t(1);
    return r;
  }

  m4_t transpose() const {
    m4_t r;
    r.m00 = m00; r.m01 = m10; r.m02 = m20; r.m03 = m30;
    r.m10 = m01; r.m11 = m11; r.m12 = m21; r.m13 = m31;
    r.m20 = m02; r.m21 = m12; r.m22 = m22; r.m23 = m32;
    r.m30 = m03; r.m31 = m13; r.m32 = m23; r.m33 = m33;
    return r;
  }

  static m4_t mul(const m4_t &a, const m4_t &b) {
    m4_t r;
    for (size_t i = 0; i < 4; ++i) {
      for (size_t j = 0; j < 4; ++j) {
        r.m[i][j] = num_t(0);
        for (size_t k = 0; k < 4; ++k) {
          r.m[i][j] += a.m[k][j] * b.m[i][k];
        }
      }
    }
    return r;
  }

  num_t operator()(int col, int row) const {
    return m[col][row];
  }

  num_t &operator()(int col, int row) {
    return m[col][row];
  }

  v3_t<num_t> transform(const v3_t<num_t> &v) const {
    return v3_t<num_t>::init(v.x * m00 + v.y * m10 + v.z * m20 + m30,
                             v.x * m01 + v.y * m11 + v.z * m21 + m31,
                             v.x * m02 + v.y * m12 + v.z * m22 + m32);
  }
};

template<typename num_t>
m4_t<num_t> operator*(const m4_t<num_t> &a, const m4_t<num_t> &b) {
  return m4_t<num_t>::mul(a, b);
}



