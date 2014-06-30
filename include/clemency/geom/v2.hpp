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
#include <tr1/functional>



template<typename _num_t>
struct v2_t {
  typedef _num_t num_t;

  union {
    struct { num_t x,y; };
    num_t v[2];
  };

  static v2_t init(const num_t *_v) {
    v2_t result;
    std::copy(_v, _v+2, result.v);
    return result;
  }

  static v2_t init(const std::vector<num_t> &_v) {
    v2_t result;
    std::copy(_v.begin(), _v.begin()+2, result.v);
    return result;
  }

  static v2_t zero() {
    return v2_t::init(num_t(0), num_t(0));
  }

  static v2_t init(num_t x = num_t(0), num_t y = num_t(0)) {
    v2_t result;
    result.x = x;
    result.y = y;
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

  v2_t abs() const {
    return v2_t::init(std::abs(x), std::abs(y));
  }

  v2_t scale(num_t s) const {
    return init(x * s, y * s);
  }

  v2_t normalize() const {
    num_t rl = num_t(1)/length();
    return init(x * rl, y * rl);
  }

  static v2_t add(const v2_t &a, const v2_t &b) {
    return init(a.x + b.x, a.y + b.y);
  }

  static v2_t add(const v2_t &a, const v2_t &b, const v2_t &c) {
    return init(a.x + b.x + c.x, a.y + b.y + c.y);
  }

  v2_t &add(const v2_t &a) {
    x += a.x; y += a.y;
    return *this;
  }

  static v2_t sub(const v2_t &a, const v2_t &b) {
    return init(a.x - b.x, a.y - b.y);
  }

  v2_t &sub(const v2_t &a) {
    x -= a.x; y -= a.y;
    return *this;
  }

  static v2_t neg(const v2_t &a) {
    return init(-a.x, -a.y);
  }

  static num_t dot(const v2_t &a, const v2_t &b) {
    return a.x * b.x + a.y * b.y;
  }

  static num_t cross(const v2_t &a, const v2_t &b) {
    return a.x * b.y - b.y * a.x;
  }

  static bool eq(const v2_t &a, const v2_t &b) {
    return std::equal(a.v, a.v+2, b.v);
  }

  static bool neq(const v2_t &a, const v2_t &b) {
    return !eq(a, b);
  }

  static bool lt(const v2_t &a, const v2_t &b) {
    return std::lexicographical_compare(a.v, a.v+2, b.v, b.v+2);
  }

  static v2_t cpow(const v2_t &a, num_t p) {
    num_t r2 = dot(a, a);
    num_t s = ::pow(r2, p / num_t(2));
    return init(s * ::cos(p * ::atan2(a.y, a.x)),
                s * ::sin(p * ::atan2(a.y, a.x)));
  }

  static v2_t csqr(const v2_t &a) {
    return init(a.x * a.x - a.y * a.y, num_t(2) * a.x * a.y);
  }

  static v2_t cinv(const v2_t &a) {
    num_t a_dot_a = dot(a, a);
    return init(a.x / a_dot_a, -a.y / a_dot_a);
  }

  static v2_t cmul(const v2_t &a, const v2_t &b) {
	return init(a.x * b.x - a.y * b.y,
                    a.x * b.y + a.y * b.x);
  }
  
  static v2_t cdiv(const v2_t &a, const v2_t &b) {
    return cmul(a, cinv(b));
  }

  /** 
   * \brief Return the orientation of c with respect to the ray defined by a->b.
   * 
   * @param[in] a 
   * @param[in] b 
   * @param[in] c 
   * 
   * @return positive, if c to the left of a->b.
   *         zero,     if c is colinear with a->b.
   *         negative, if c to the right of a->b.
   */
  static num_t orient(const v2_t &a, const v2_t &b, const v2_t &c) {
    num_t acx = a.x - c.x;
    num_t bcx = b.x - c.x;
    num_t acy = a.y - c.y;
    num_t bcy = b.y - c.y;
    return acx * bcy - acy * bcx;
  }
};

template<typename num_t>
v2_t<num_t> operator*(num_t a, const v2_t<num_t> &b) {
  return b.scale(a);
}

template<typename num_t>
static  inline v2_t<num_t> operator/(const v2_t<num_t> &a, num_t b) {
  return a.scale(num_t(1)/b);
}

template<typename num_t>
v2_t<num_t> operator+(const v2_t<num_t> &a, const v2_t<num_t> &b) {
  return v2_t<num_t>::add(a, b);
}

template<typename num_t>
v2_t<num_t> operator-(const v2_t<num_t> &a, const v2_t<num_t> &b) {
  return v2_t<num_t>::sub(a, b);
}

template<typename num_t>
v2_t<num_t> &operator+=(v2_t<num_t> &a, const v2_t<num_t> &b) {
  return a.add(b);
}

template<typename num_t>
v2_t<num_t> &operator-=(v2_t<num_t> &a, const v2_t<num_t> &b) {
  return a.sub(b);
}

template<typename num_t>
v2_t<num_t> operator-(const v2_t<num_t> &a) {
  return v2_t<num_t>::neg(a);
}

template<typename num_t>
bool operator==(const v2_t<num_t> &a, const v2_t<num_t> &b) {
  return v2_t<num_t>::eq(a, b);
}

template<typename num_t>
bool operator!=(const v2_t<num_t> &a, const v2_t<num_t> &b) {
  return !v2_t<num_t>::eq(a, b);
}

template<typename num_t>
bool operator<(const v2_t<num_t> &a, const v2_t<num_t> &b) {
  return v2_t<num_t>::lt(a, b);
}

template<typename num_t>
bool operator>=(const v2_t<num_t> &a, const v2_t<num_t> &b) {
  return !v2_t<num_t>::lt(a, b);
}

template<typename num_t>
bool operator>(const v2_t<num_t> &a, const v2_t<num_t> &b) {
  return v2_t<num_t>::lt(b, a);
}

template<typename num_t>
bool operator<=(const v2_t<num_t> &a, const v2_t<num_t> &b) {
  return !v2_t<num_t>::lt(b, a);
}

template<typename num_t>
std::ostream &operator<<(std::ostream &out, const v2_t<num_t> &v) {
  out << "<" << v.x << "," << v.y << ">";
  return out;
}

struct v2_hash_t {
  template<typename num_t>
  size_t operator()(const v2_t<num_t> &a) const {
    std::tr1::hash<num_t> h;
    size_t r = 0;
    r *= 131; r ^= h(a.x);
    r *= 131; r ^= h(a.y);
    return r;
  }

  template<typename num_t>
  size_t operator()(const std::pair<v2_t<num_t>, v2_t<num_t> > &a) const {
    std::tr1::hash<num_t> h;
    size_t r = 0;
    r *= 131; r ^= h(a.first.x);
    r *= 131; r ^= h(a.first.y);
    r *= 131; r ^= h(a.second.x);
    r *= 131; r ^= h(a.second.y);
    return r;
  }
};



typedef v2_t<double> v2d_t;
typedef v2_t<float>  v2f_t;
