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
#include "circle.hpp"



template<typename num_t>
class aabb2_t {
public:
  v2_t<num_t> mid, extent;

  static aabb2_t init(const v2_t<num_t> &mid, const v2_t<num_t> &extent) {
    aabb2_t result;
    result.mid = mid;
    result.extent = extent;
    return result;
  }

  static aabb2_t initWithPoints(const v2_t<num_t> &a, const v2_t<num_t> &b, const v2_t<num_t> &c) {
    aabb2_t result;
    result._set_dim<0>(a.x, b.x, c.x);
    result._set_dim<1>(a.y, b.y, c.y);
    return result;
  }

  static aabb2_t initWithPoints(const v2_t<num_t> &lo, const v2_t<num_t> &hi) {
    aabb2_t result;
    result._set_dim<0>(lo.x, hi.x);
    result._set_dim<1>(lo.y, hi.y);
    return result;
  }

  template<unsigned dim>
  void _set_dim(num_t a, num_t b, num_t c) {
    _set_dim<dim>(std::min(std::min(a, b), c), std::max(std::max(a, b), c));
  }

  template<unsigned dim>
  void _set_dim(num_t l, num_t h) {
    mid.v[dim] = (l + h) / num_t(2);
    extent.v[dim] = std::max(h - mid.v[dim], mid.v[dim] - l);
  }

  aabb2_t &aabb() { return *this; }
  const aabb2_t &aabb() const { return *this; }

  aabb2_t &expand(const aabb2_t &other) {
    if (other.xl() < xl() || other.xh() > xh()) {
      _set_dim<0>(std::min(other.xl(), xl()), std::max(other.xh(), xh()));
    }

    if (other.yl() < yl() || other.yh() > yh()) {
      _set_dim<1>(std::min(other.yl(), yl()), std::max(other.yh(), yh()));
    }

    return *this;
  }

  num_t _l(size_t dim) const { return mid.v[dim] - extent.v[dim]; }
  num_t _m(size_t dim) const { return mid.v[dim]; }
  num_t _h(size_t dim) const { return mid.v[dim] + extent.v[dim]; }

  num_t xl() const { return mid.x - extent.x; }
  num_t xm() const { return mid.x;            }
  num_t xh() const { return mid.x + extent.x; }

  num_t yl() const { return mid.y - extent.y; }
  num_t ym() const { return mid.y;            }
  num_t yh() const { return mid.y + extent.y; }

  num_t area() const {
    return (xh()-xl()) * (yh()-yl());
  }

  v2_t<num_t> corner(unsigned c) const {
    return v2_t<num_t>::init(mid.x + extent.x * (c & 1 ? num_t(+1) : num_t(-1)),
                             mid.y + extent.y * (c & 2 ? num_t(+1) : num_t(-1)));
  }

  bool touchesFace(const aabb2_t &bbox) const {
    num_t dx = std::abs(mid.x - bbox.mid.x) - extent.x + bbox.extent.x;
    num_t dy = std::abs(mid.y - bbox.mid.y) - extent.y + bbox.extent.y;
    return
      (dx == num_t(0) && dy < num_t(0)) ||
      (dy == num_t(0) && dx < num_t(0));
  }

  bool intersects(const aabb2_t &bbox) const {
    return std::abs(mid.x - bbox.mid.x) < extent.x + bbox.extent.x ||
           std::abs(mid.y - bbox.mid.y) < extent.y + bbox.extent.y;
  }

  bool contains(const aabb2_t &bbox) const {
    return std::abs(mid.x - bbox.mid.x) < extent.x - bbox.extent.x &&
           std::abs(mid.y - bbox.mid.y) < extent.y - bbox.extent.y;
  }

  bool contains(const v2_t<num_t> &p) const {
    return
      std::abs(p.x - mid.x) < extent.x &&
      std::abs(p.y - mid.y) < extent.y;
  }

  bool intersects(const v2_t<num_t> &p) const {
    return contains(p);
  }

  bool contains(const circle_t<num_t> &sp) const {
    return
      std::abs(sp.v.x - mid.x) < extent.x - sp.r &&
      std::abs(sp.v.y - mid.y) < extent.y - sp.r;
  }

  bool intersects(const circle_t<num_t> &sp) const {
    num_t r = num_t(0);
    for (unsigned i = 0; i < 2; ++i) {
      num_t t = std::abs(sp.v.v[i] - mid.v[i]) - extent.v[i]; if (t > num_t(0)) r += t*t;
    }
    return r <= sp.r * sp.r;
  }

  num_t distancesq(const v2_t<num_t> &p) const {
    num_t dist = num_t(0);
    for (size_t i = 0; i < 2; ++i) {
      num_t d = std::abs(p.v[i] - mid.v[i]) - extent.v[i];
      if (d > num_t(0)) dist += d*d;
    }
    return dist;
  }

  num_t distance(const v2_t<num_t> &p) const {
    return sqrt(distancesq(p));
  }
};

template<typename num_t>
std::ostream &operator<<(std::ostream &out, const aabb2_t<num_t> &aabb) {
  out << aabb.xl() << "," << aabb.yl() << ":" << aabb.xh() << "," << aabb.yh();
  return out;
}

template<typename rng_t, typename num_t>
v2_t<num_t> randomInAABB(rng_t &rng, const aabb2_t<num_t> &aabb) {
  boost::uniform_01<num_t> uniform;
  boost::variate_generator<rng_t &, boost::uniform_01<num_t> > gen(rng, uniform);
  return v2_t<num_t>::init(aabb.mid.x + (gen() - num_t(.5)) * aabb.extent.x,
                           aabb.mid.y + (gen() - num_t(.5)) * aabb.extent.y);
}
