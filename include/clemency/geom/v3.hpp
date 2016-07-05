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
#include <iostream>
#include <vector>
#include <functional>
#include <algorithm>

#include <clemency/geom/v2.hpp>

template<typename num_t>
struct v3_t {
  union {
    struct { num_t x,y,z; };
    num_t v[3];
  };

  static v3_t init(const num_t *_v) {
    v3_t result;
    std::copy(_v, _v+3, result.v);
    return result;
  }

  static v3_t init(const std::vector<num_t> &_v) {
    v3_t result;
    std::copy(_v.begin(), _v.begin()+3, result.v);
    return result;
  }

  static v3_t zero() {
    return v3_t::init(num_t(0), num_t(0), num_t(0));
  }

  static v3_t init(num_t x = num_t(0), num_t y = num_t(0), num_t z = num_t(0)) {
    v3_t result;
    result.x = x;
    result.y = y;
    result.z = z;
    return result;
  }

  num_t operator[](int i) const {
    return v[i];
  }

  num_t &operator[](int i) {
    return v[i];
  }

  num_t length() const {
    return ::sqrt(dot(*this, *this));
  }

  num_t lengthsq() const {
    return dot(*this, *this);
  }

  v3_t abs() const {
    return v3_t::init(std::abs(x), std::abs(y), std::abs(z));
  }

  v2_t<num_t> project(int axis) const {
    return v2_t<num_t>::init(v[(axis+1)%3], v[(axis+2)%3]);
  }

  int max_component() const {
    v3_t a = abs();
    return std::max_element(a.v, a.v+3) - a.v;
  }

  int min_component() const {
    v3_t a = abs();
    return std::min_element(a.v, a.v+3) - a.v;
  }

  static v3_t componentwise_min(const v3_t &a, const v3_t &b) {
    return init(std::min(a.x, b.x), std::min(a.y, b.y), std::min(a.z, b.z));
  }

  static v3_t componentwise_max(const v3_t &a, const v3_t &b) {
    return init(std::max(a.x, b.x), std::max(a.y, b.y), std::max(a.z, b.z));
  }

  v3_t scale(num_t s) const {
    return init(x * s, y * s, z * s);
  }

  v3_t normalize() const {
    num_t rl = num_t(1)/length();
    return init(x * rl, y * rl, z * rl);
  }

  static v3_t add(const v3_t &a, const v3_t &b) {
    return init(a.x + b.x, a.y + b.y, a.z + b.z);
  }

  static v3_t add(const v3_t &a, const v3_t &b, const v3_t &c) {
    return init(a.x + b.x + c.x, a.y + b.y + c.y, a.z + b.z + c.z);
  }

  v3_t &add(const v3_t &a) {
    x += a.x; y += a.y; z += a.z;
    return *this;
  }

  static v3_t sub(const v3_t &a, const v3_t &b) {
    return init(a.x - b.x, a.y - b.y, a.z - b.z);
  }

  v3_t &sub(const v3_t &a) {
    x -= a.x; y -= a.y; z -= a.z;
    return *this;
  }

  static v3_t neg(const v3_t &a) {
    return init(-a.x, -a.y, -a.z);
  }

  static num_t dot(const v3_t &a, const v3_t &b) {
    return a.x * b.x + a.y * b.y + a.z * b.z;
  }

  static v3_t cross(const v3_t &a, const v3_t &b) {
    return init(+(a.y * b.z - b.y * a.z),
                -(a.x * b.z - b.x * a.z),
                +(a.x * b.y - b.x * a.y));
  }

  static bool eq(const v3_t &a, const v3_t &b) {
    return std::equal(a.v, a.v+3, b.v);
  }

  static bool neq(const v3_t &a, const v3_t &b) {
    return !eq(a, b);
  }

  static bool lt(const v3_t &a, const v3_t &b) {
    return std::lexicographical_compare(a.v, a.v+3, b.v, b.v+3);
  }

  // Compute a . (b x c)
  static num_t dotcross(const v3_t &a, const v3_t &b, const v3_t &c);

  // test whether point d is above, below or on the plane formed by
  // the triangle a,b,c.
  // return: +ve = d is below a,b,c
  //         -ve = d is above a,b,c
  //           0 = d is on a,b,c
  static num_t orient(const v3_t &a, const v3_t &b, const v3_t &c, const v3_t &d);

  // Volume of a tetrahedron described by 4 points. Will be positive
  // if the anticlockwise normal of a,b,c is oriented out of the
  // tetrahedron.
  //
  // see: http://mathworld.wolfram.com/Tetrahedron.html
  static num_t tetrahedron_volume(const v3_t &a, const v3_t &b, const v3_t &c, const v3_t &d);
};

template<typename num_t>
static inline v3_t<num_t> operator*(num_t a, const v3_t<num_t> &b) {
  return b.scale(a);
}

template<typename num_t>
static  inline v3_t<num_t> operator/(const v3_t<num_t> &a, num_t b) {
  return a.scale(num_t(1)/b);
}

template<typename num_t>
static inline v3_t<num_t> operator+(const v3_t<num_t> &a, const v3_t<num_t> &b) {
  return v3_t<num_t>::add(a, b);
}

template<typename num_t>
static inline v3_t<num_t> operator-(const v3_t<num_t> &a, const v3_t<num_t> &b) {
  return v3_t<num_t>::sub(a, b);
}

template<typename num_t>
static inline v3_t<num_t> &operator+=(v3_t<num_t> &a, const v3_t<num_t> &b) {
  return a.add(b);
}

template<typename num_t>
static inline v3_t<num_t> &operator-=(v3_t<num_t> &a, const v3_t<num_t> &b) {
  return a.sub(b);
}

template<typename num_t>
static inline v3_t<num_t> operator-(const v3_t<num_t> &a) {
  return v3_t<num_t>::neg(a);
}

template<typename num_t>
static inline bool operator==(const v3_t<num_t> &a, const v3_t<num_t> &b) {
  return v3_t<num_t>::eq(a, b);
}

template<typename num_t>
static inline bool operator!=(const v3_t<num_t> &a, const v3_t<num_t> &b) {
  return !v3_t<num_t>::eq(a, b);
}

template<typename num_t>
static inline bool operator<(const v3_t<num_t> &a, const v3_t<num_t> &b) {
  return v3_t<num_t>::lt(a, b);
}

template<typename num_t>
static inline bool operator>=(const v3_t<num_t> &a, const v3_t<num_t> &b) {
  return !v3_t<num_t>::lt(a, b);
}

template<typename num_t>
static inline bool operator>(const v3_t<num_t> &a, const v3_t<num_t> &b) {
  return v3_t<num_t>::lt(b, a);
}

template<typename num_t>
static inline bool operator<=(const v3_t<num_t> &a, const v3_t<num_t> &b) {
  return !v3_t<num_t>::lt(b, a);
}

template<typename num_t>
static inline std::ostream &operator<<(std::ostream &out, const v3_t<num_t> &v) {
  out << "<" << v.x << "," << v.y << "," << v.z << ">";
  return out;
}

struct v3_hash_t {
  template<typename num_t>
  size_t operator()(const v3_t<num_t> &a) const {
    std::hash<num_t> h;
    size_t r = 0;
    r *= 131; r ^= h(a.x);
    r *= 131; r ^= h(a.y);
    r *= 131; r ^= h(a.z);
    return r;
  }

  template<typename num_t>
  size_t operator()(const std::pair<v3_t<num_t>, v3_t<num_t> > &a) const {
    std::hash<num_t> h;
    size_t r = 0;
    r *= 131; r ^= h(a.first.x);
    r *= 131; r ^= h(a.first.y);
    r *= 131; r ^= h(a.first.z);
    r *= 131; r ^= h(a.second.x);
    r *= 131; r ^= h(a.second.y);
    r *= 131; r ^= h(a.second.z);
    return r;
  }
};

template<typename num_t>
num_t v3_t<num_t>::dotcross(const v3_t<num_t> &a, const v3_t<num_t> &b, const v3_t<num_t> &c) {
  return
    (a.x * b.y * c.z + a.y * b.z * c.x + a.z * b.x * c.y) -
    (a.x * c.y * b.z + a.y * c.z * b.x + a.z * c.x * b.y);
}

template<typename num_t>
num_t v3_t<num_t>::orient(const v3_t<num_t> &a, const v3_t<num_t> &b, const v3_t<num_t> &c, const v3_t<num_t> &d) {
  return dotcross((a - d), (b - d), (c - d));
}

template<typename num_t>
num_t v3_t<num_t>::tetrahedron_volume(const v3_t<num_t> &a, const v3_t<num_t> &b, const v3_t<num_t> &c, const v3_t<num_t> &d) {
  return dotcross((a - d), (b - d), (c - d)) / num_t(6);
}





typedef v3_t<double> v3d_t;
typedef v3_t<float>  v3f_t;
