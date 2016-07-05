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

#include <clemency/geom/v2.hpp>

template <typename num_t>
class tri2_t {
 public:
  union {
    v2_t<num_t> v[3];
    struct {
      v2_t<num_t> a, b, c;
    };
  };

  tri2_t(const v2_t<num_t>& _a, const v2_t<num_t>& _b, const v2_t<num_t>& _c)
      : a(_a), b(_b), c(_c) {}

  // not implemented
  bool contains(const v2_t<num_t>& p) const;

  // not implemented
  v2_t<num_t> closest_point(const v2_t<num_t>& p) const;

  double distancesq(const v2_t<num_t>& p) const {
    return (p - closest_point(p)).lengthsq();
  }

  double distance(const v2_t<num_t>& p) const {
    return (p - closest_point(p)).length();
  }
};

typedef tri2_t<double> tri2d_t;
typedef tri2_t<float> tri2f_t;
