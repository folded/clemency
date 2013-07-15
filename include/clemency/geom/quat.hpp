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



template<typename num_t>
class quat_t {
public:
  num_t a;
  v3_t<num_t> b;

  static quat_t zero() {
    return quat_t::init(0.0, v3_t<num_t>::zero());
  }

  static quat_t init(num_t a, const v3_t<num_t> &b) {
    quat_t result;
    result.a = a;
    result.b = b;
    return result;
  }

  static quat_t rot(num_t angle, const v3_t<num_t> &axis) {
    quat_t r;
    num_t len = axis.length();
    if (fabs(len) > 1e-10) {
      num_t omega = -0.5 * angle;
      num_t c = cos(omega);
      num_t s = sin(omega);

      r.a = c;
      r.b.x = axis.x * s / len;
      r.b.y = axis.y * s / len;
      r.b.z = axis.z * s / len;

      r.normalize();
    } else {
      r.a = 1.0;
      r.b.x = r.b.y = r.b.z = 0.0;
    }
    return r;
  }

  quat_t sq() const {
    return init(a * a - b.lengthsq(), b.scale(2.0 * a));
  }

  num_t length() const {
    return ::sqrt(quat_t::dot(*this, *this));
  }

  num_t lengthsq() const {
    return quat_t::dot(*this, *this);
  }

  quat_t scale(num_t s) const {
    return init(a * s, b.scale(s));
  }

  quat_t normalize() const {
    return scale(1.0/length());
  }

  static quat_t add(const quat_t &a, const quat_t &b) {
    return init(a.a + b.a, a.b + b.b);
  }

  static quat_t sub(const quat_t &a, const quat_t &b) {
    return init(a.a - b.a, a.b - b.b);
  }

  static quat_t add(const quat_t &a, const quat_t &b, const quat_t &c) {
    return init(a.a + b.a + c.a, a.b + b.b + c.b);
  }

  static num_t dot(const quat_t &a, const quat_t &b) {
    return a.a * b.a + v3_t<num_t>::dot(a.b, b.b);
  }

  static quat_t mul(const quat_t &a, const quat_t &b) {
    return init(a.a * b.a - v3_t<num_t>::dot(a.b, b.b),
                b.b.scale(a.a) + a.b.scale(b.a) + v3_t<num_t>::cross(a.b, b.b));
  }
};

template<typename num_t>
quat_t<num_t> operator*(const quat_t<num_t> &a, const quat_t<num_t> &b) {
  return quat_t<num_t>::mul(a, b);
}

template<typename num_t>
static  inline quat_t<num_t> operator*(num_t a, const quat_t<num_t> &b) {
  return b.scale(a);
}

template<typename num_t>
quat_t<num_t> operator+(const quat_t<num_t> &a, const quat_t<num_t> &b) {
  return quat_t<num_t>::add(a, b);
}

template<typename num_t>
quat_t<num_t> operator-(const quat_t<num_t> &a, const quat_t<num_t> &b) {
  return quat_t<num_t>::sub(a, b);
}

template<typename num_t>
std::ostream &operator<<(std::ostream &out, const quat_t<num_t> &q) {
  out << "<" << q.a << ";" << q.b.x << "," << q.b.y << "," << q.b.z << ">";
  return out;
}



typedef quat_t<double> quatd_t;
typedef quat_t<float>  quatf_t;
