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



#include "v2.hpp"
#include "solve.hpp"



template<typename num_t>
struct qbezier2_t {
  v2_t<num_t> p0, p1, p2;
  num_t A, B, C;

  void init(const v2_t<num_t> &_p0, const v2_t<num_t> &_p1, const v2_t<num_t> &_p2) {
    p0 = _p0;
    p1 = _p1 - _p0;
    p2 = _p2 - _p0;

    A = num_t(+4) * v2_t<num_t>::dot(p1, p1) - num_t(4) * v2_t<num_t>::dot(p1, p2) + v2_t<num_t>::dot(p2, p2);
    B = num_t(-6) * v2_t<num_t>::dot(p1, p1) + num_t(3) * v2_t<num_t>::dot(p1, p2);
    C = num_t(+2) * v2_t<num_t>::dot(p1, p1);
  }

  qbezier2_t() {
    init(v2_t<num_t>::zero(), v2_t<num_t>::zero(), v2_t<num_t>::zero());
  }

  qbezier2_t(const v2_t<num_t> &_p0, const v2_t<num_t> &_p1, const v2_t<num_t> &_p2) {
    init(_p0, _p1, _p2);
  }

  v2_t<num_t> get_relpos(num_t t) const {
    return num_t(2) * t * (num_t(1) - t) * p1 + t*t * p2;
  }

  v2_t<num_t> get_pos(num_t t) const {
    return p0 + get_relpos(t);
  }

  v2_t<num_t> get_tangent(num_t t) const {
    return num_t(2) * (num_t(1) - t) * p1 + num_t(2) * t * (p2 - p1);
  }

  std::pair<v2_t<num_t>, num_t> closest_point(v2_t<num_t> pos) const {
    v2_t<num_t> p = pos - p0;
    num_t p1p = v2_t<num_t>::dot(p1, p);
    num_t p2p = v2_t<num_t>::dot(p2, p);

    std::vector<num_t> sol = solve3(A,
                                    B,
                                    C - p2p + num_t(2) * p1p,
                                    -p1p);
    num_t d;
    num_t min_dist = p.lengthsq();
    v2_t<num_t> min_pos = v2_t<num_t>::zero();
    num_t min_t = num_t(0);

    d = (p - p2).lengthsq();
    if (d < min_dist) {
      min_dist = d;
      min_pos = p2;
      min_t = num_t(1);
    }

    for (size_t i = 0; i < sol.size(); ++i) {
      if (sol[i] > num_t(0) && sol[i] < num_t(1)) {
        v2_t<num_t> q = get_relpos(sol[i]);
        d = (p - q).lengthsq();
        if (d < min_dist) {
          min_dist = d;
          min_pos = q;
          min_t = sol[i];
        }
      }
    }

    min_dist = ::sqrt(min_dist);
    if (v2_t<num_t>::orient(min_pos, min_pos + get_tangent(min_t), p) < num_t(0)) {
      min_dist = -min_dist;
    }

    return std::make_pair(min_pos + p0, min_dist);
  }
};



typedef qbezier2_t<double> qbezier2d_t;
typedef qbezier2_t<float>  qbezier2f_t;
