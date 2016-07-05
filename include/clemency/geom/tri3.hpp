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



#include <clemency/geom/v3.hpp>



template<typename num_t>
class tri3_t {
public:
  union {
    v3_t<num_t> v[3];
    struct {
      v3_t<num_t> a,b,c;
    };
  };

  tri3_t(const v3_t<num_t> &_a, const v3_t<num_t> &_b, const v3_t<num_t> &_c) : a(_a), b(_b), c(_c) {
  }

  double orient(const v3_t<num_t> &p) const {
    return v3_t<num_t>::orient(a, b, c, p);
  }

  v3_t<num_t> closest_point(const v3_t<num_t> &p) const;

  double distancesq(const v3_t<num_t> &p) const {
    return (p - closest_point(p)).lengthsq();
  }

  double distance(const v3_t<num_t> &p) const {
    return (p - closest_point(p)).length();
  }
};

namespace {
  template<typename num_t>
  static inline num_t clamp(num_t v, num_t lo, num_t hi) {
    return std::min(std::max(v, lo), hi);
  }
}

template<typename num_t>
v3_t<num_t> tri3_t<num_t>::closest_point(const v3_t<num_t> &p) const {
  const v3_t<num_t> e0 = v[1] - v[0];
  const v3_t<num_t> e1 = v[2] - v[0];
  const v3_t<num_t> dp = v[0] - p;

  // triangle is defined by
  // T(s,t) = v[0] + e0.s + e1.t; s >= 0, t >= 0, s + t <= 1
  // distance from p to point on triangle plane:
  // Q(s, t) = a(s^2) + 2b(st) + c(t^2) + 2d(s) + 2e(t) + f
  // grad(Q) = 2(as + bt + d, bs + ct + e)

  const num_t a = v3_t<num_t>::dot(e0, e0);
  const num_t b = v3_t<num_t>::dot(e1, e0);
  const num_t c = v3_t<num_t>::dot(e1, e1);
  const num_t d = v3_t<num_t>::dot(e0, dp);
  const num_t e = v3_t<num_t>::dot(e1, dp);
  const num_t f = v3_t<num_t>::dot(dp, dp);

  const num_t det = a*c - b*b;

  num_t s = b*e - c*d;
  num_t t = b*d - a*e;

  /*      t             */
  /*      ^             */
  /*    \4|             */
  /*     \|             */
  /*      \             */
  /*      |\6           */
  /*     1|3\           */
  /*   ---+--\----> s   */
  /*     0|2  \5        */

  int edge;

  if (s+t <= det) { // regions 0, 1, 2, 3
    if (s < 0 && t < 0) {
      // region 0 - minimum on edge (1) or (2)
      if (d < 0) {
        // grad(Q)(0,0).(1,0) < 0
        edge = 1;
      } else {
        edge = 2;
      }
    } else if (s < 0) {
      // region 1 - minimum on edge (1)
      edge = 1;
    } else if (t < 0) {
      // region 2 - minimum on edge (2)
      edge = 2;
    } else {
      // region 3 - closest point within triangle bounds.
      edge = 0;
    }
  } else { // regions 4, 5, 6
    if (s < 0) {
      // region 4 - minimum on edge (1) or (3)
      if (-(c+e) < 0) {
        // grad(Q)(0,1).(0,-1) < 0
        edge = 1;
      } else {
        edge = 3;
      }
    } else if (t < 0) {
      // region 5 - minimum on edge (2) or (3)
      if (-(a+d) < 0) {
        // grad(Q)(1,0).(-1,0) < 0
        edge = 2;
      } else {
        edge = 3;
      }
    } else {
      // region 6 - minimum on edge (3)
      edge = 3;
    }
  }

  switch (edge) {
  case 0: {
    s /= det;
    t /= det;
    break;
  }
  case 1: {
    // edge (1)
    // s = 0, t = [0,1]
    // :: Q(t)  = c(t^2) + 2e(t) + f
    // :: dQ/dt = 2ct + 2e
    s = 0;
    t = clamp(-e/c, 0.0, 1.0);
    break;
  }
  case 2: {
    // edge (2)
    // t = 0, s = [0,1]
    // :: Q(s)  = a(s^2) + 2d(s) + f
    // :: dQ/ds = 2as + 2d
    s = clamp(-d/a, 0.0, 1.0);
    t = 0;
    break;
  }
  case 3: {
    // edge (3)
    // t = 1 - s, s = [0,1]
    // Q(s)  = a(s^2) + 2b(s(1-s)) + c((1-s)^2) + 2d(s) + 2e(1-s) + f
    //       = a(s^2) + 2b(s-s^2) + c(1-2s+s^2) + 2d(s) + 2e(1-s) + f
    //       = a(s^2) + 2b(s) - 2b(s^2) + c - 2c(s) + c(s^2) + 2d(s) + 2e - 2e(s) + f
    //       = (a - 2b + c)(s^2) + (2b - 2c + 2d - 2e)(s) + (c + 2e + f)
    // dQ/ds = 2(a - 2b + c)s + 2(b - c + d - e)
    s = clamp((c+e-b-d)/(a-2*b+c), 0.0, 1.0);
    t = 1 - s;
    break;
  }
  }

  const v3_t<num_t> closest = v[0] + s*e0 + t*e1;

  return closest;
}



typedef tri3_t<double> tri3d_t;
typedef tri3_t<float>  tri3f_t;
