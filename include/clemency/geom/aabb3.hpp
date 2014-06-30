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
#include "tri3.hpp"
#include "sphere.hpp"



template<typename num_t>
class aabb3_t {
public:
  v3_t<num_t> mid, extent;

  static aabb3_t init(const v3_t<num_t> &mid, const v3_t<num_t> &extent) {
    aabb3_t result;
    result.mid = mid;
    result.extent = extent;
    return result;
  }

  static aabb3_t initWithPoints(const v3_t<num_t> &a, const v3_t<num_t> &b, const v3_t<num_t> &c) {
    aabb3_t result;
    result._set_dim<0>(a.x, b.x, c.x);
    result._set_dim<1>(a.y, b.y, c.y);
    result._set_dim<2>(a.z, b.z, c.z);
    return result;
  }

  static aabb3_t initWithPoints(const v3_t<num_t> &a, const v3_t<num_t> &b) {
    aabb3_t result;
    result._set_dim<0>(a.x, b.x);
    result._set_dim<1>(a.y, b.y);
    result._set_dim<2>(a.z, b.z);
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

  aabb3_t &aabb() { return *this; }
  const aabb3_t &aabb() const { return *this; }

  aabb3_t &expand(const aabb3_t &other) {
    if (other.xl() < xl() || other.xh() > xh()) {
      _set_dim<0>(std::min(other.xl(), xl()), std::max(other.xh(), xh()));
    }

    if (other.yl() < yl() || other.yh() > yh()) {
      _set_dim<1>(std::min(other.yl(), yl()), std::max(other.yh(), yh()));
    }

    if (other.zl() < zl() || other.zh() > zh()) {
      _set_dim<2>(std::min(other.zl(), zl()), std::max(other.zh(), zh()));
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

  num_t zl() const { return mid.z - extent.z; }
  num_t zm() const { return mid.z;            }
  num_t zh() const { return mid.z + extent.z; }

  num_t volume() const {
    return (xh()-xl()) * (yh()-yl()) * (zh()-zl());
  }

  v3_t<num_t> corner(unsigned c) const {
    return v3_t<num_t>::init(mid.x + extent.x * (c & 1 ? num_t(+1) : num_t(-1)),
                      mid.y + extent.y * (c & 2 ? num_t(+1) : num_t(-1)),
                      mid.z + extent.z * (c & 4 ? num_t(+1) : num_t(-1)));
  }

  bool touchesFace(const aabb3_t &bbox) const {
    num_t dx = std::abs(mid.x - bbox.mid.x) - extent.x + bbox.extent.x;
    num_t dy = std::abs(mid.y - bbox.mid.y) - extent.y + bbox.extent.y;
    num_t dz = std::abs(mid.z - bbox.mid.z) - extent.z + bbox.extent.z;
    return
      (dx == num_t(0) && dy < num_t(0) && dz < num_t(0)) ||
      (dy == num_t(0) && dx < num_t(0) && dz < num_t(0)) ||
      (dz == num_t(0) && dx < num_t(0) && dy < num_t(0));
  }

  bool intersects(const aabb3_t &bbox) const {
    return std::abs(mid.x - bbox.mid.x) < extent.x + bbox.extent.x ||
           std::abs(mid.y - bbox.mid.y) < extent.y + bbox.extent.y ||
           std::abs(mid.z - bbox.mid.z) < extent.z + bbox.extent.z;
  }

  bool contains(const aabb3_t &bbox) const {
    return std::abs(mid.x - bbox.mid.x) < extent.x - bbox.extent.x &&
           std::abs(mid.y - bbox.mid.y) < extent.y - bbox.extent.y &&
           std::abs(mid.z - bbox.mid.z) < extent.z - bbox.extent.z;
  }

  bool contains(const v3_t<num_t> &p) const {
    return
      std::abs(p.x - mid.x) < extent.x &&
      std::abs(p.y - mid.y) < extent.y &&
      std::abs(p.z - mid.z) < extent.z;
  }

  bool intersects(const v3_t<num_t> &p) const {
    return contains(p);
  }

  bool contains(const tri3_t<num_t> &tri) const {
    return contains(tri.a) && contains(tri.b) && contains(tri.c);
  }

  bool intersects(tri3_t<num_t> tri) const;

  bool contains(const sphere_t<num_t> &sp) const {
    return
      std::abs(sp.v.x - mid.x) < extent.x - sp.r &&
      std::abs(sp.v.y - mid.y) < extent.y - sp.r &&
      std::abs(sp.v.z - mid.z) < extent.z - sp.r;
  }

  bool intersects(const sphere_t<num_t> &sp) const {
    num_t r = num_t(0);
    for (unsigned i = 0; i < 3; ++i) {
      num_t t = std::abs(sp.v.v[i] - mid.v[i]) - extent.v[i]; if (t > num_t(0)) r += t*t;
    }
    return r <= sp.r * sp.r;
  }

  num_t distancesq(const v3_t<num_t> &p) const {
    num_t dist = num_t(0);
    for (size_t i = 0; i < 3; ++i) {
      num_t d = std::abs(p.v[i] - mid.v[i]) - extent.v[i];
      if (d > num_t(0)) dist += d*d;
    }
    return dist;
  }

  num_t distance(const v3_t<num_t> &p) const {
    return sqrt(distancesq(p));
  }
};

template<typename num_t>
std::ostream &operator<<(std::ostream &out, const aabb3_t<num_t> &aabb) {
  out << aabb.xl() << "," << aabb.yl() << "," << aabb.zl() << ":" << aabb.xh() << "," << aabb.yh() << "," << aabb.zh();
  return out;
}

template<typename rng_t, typename num_t>
std::pair<v3_t<num_t>, v3_t<num_t> > randomLineSegmentInAABB(rng_t &rng, const aabb3_t<num_t> &aabb) {
  std::pair<v3_t<num_t>, v3_t<num_t> > result;
  result.first = randomOnAABB(rng, aabb);
  result.second = randomOnAABB(rng, aabb);
  return result;
}

template<typename rng_t, typename num_t>
v3_t<num_t> randomOnAABB(rng_t &rng, const aabb3_t<num_t> &aabb) {
  boost::random::uniform_01<num_t> uniform;
  boost::random::uniform_int_distribution<> face(0, 5);
  boost::random::variate_generator<rng_t &, boost::uniform_01<num_t> > gen(rng, uniform);
  v3_t<num_t> pos;
  int f = face(rng);

  if (f != 0 && f != 3) {
    pos.x = aabb.mid.x + (num_t(-1) + gen() * num_t(2)) * aabb.extent.x;
  } else if (f == 0) {
    pos.x = aabb.mid.x - aabb.extent.x;
  } else if (f == 3) {
    pos.x = aabb.mid.x + aabb.extent.x;
  }

  if (f != 1 && f != 4) {
    pos.y = aabb.mid.y + (num_t(-1) + gen() * num_t(2)) * aabb.extent.y;
  } else if (f == 1) {
    pos.y = aabb.mid.y - aabb.extent.y;
  } else if (f == 4) {
    pos.y = aabb.mid.y + aabb.extent.y;
  }

  if (f != 2 && f != 5) {
    pos.z = aabb.mid.z + (num_t(-1) + gen() * num_t(2)) * aabb.extent.z;
  } else if (f == 2) {
    pos.z = aabb.mid.z - aabb.extent.z;
  } else if (f == 5) {
    pos.z = aabb.mid.z + aabb.extent.z;
  }

  return pos;
}

template<typename rng_t, typename num_t>
v3_t<num_t> randomInAABB(rng_t &rng, const aabb3_t<num_t> &aabb) {
  boost::uniform_01<num_t> uniform;
  boost::variate_generator<rng_t &, boost::uniform_01<num_t> > gen(rng, uniform);
  return v3_t<num_t>::init(aabb.mid.x + (gen() - num_t(0.5)) * aabb.extent.x,
                           aabb.mid.y + (gen() - num_t(0.5)) * aabb.extent.y,
                           aabb.mid.z + (gen() - num_t(0.5)) * aabb.extent.z);
}

template<int Ax, int Ay, int Az, int c, typename num_t>
bool intersectsTriangle_axisTest_3(const aabb3_t<num_t> &aabb, const tri3_t<num_t> &tri) {
  const int d = (c+1) % 3, e = (c+2) % 3;
  const v3_t<num_t> a = v3_t<num_t>::cross(v3_t<num_t>::init(Ax, Ay, Az), tri.v[d] - tri.v[c]);
  num_t p1 = v3_t<num_t>::dot(a, tri.v[c]);
  num_t p2 = v3_t<num_t>::dot(a, tri.v[e]);
  if (p1 > p2) std::swap(p1, p2);
  const num_t r = v3_t<num_t>::dot(a.abs(), aabb.extent);
  return !(p1 > r || p2 < -r);
}

template<int c, typename num_t>
bool intersectsTriangle_axisTest_2(const aabb3_t<num_t> &aabb, const tri3_t<num_t> &tri) {
  num_t vmin = std::min(std::min(tri.v[0][c], tri.v[1][c]), tri.v[2][c]);
  num_t vmax = std::max(std::max(tri.v[0][c], tri.v[1][c]), tri.v[2][c]);
  return !(vmin > aabb.extent[c] || vmax < -aabb.extent[c]);
}

template<typename num_t>
bool intersectsTriangle_axisTest_1(const aabb3_t<num_t> &aabb, const tri3_t<num_t> &tri) {
  v3_t<num_t> n = v3_t<num_t>::cross(tri.v[1] - tri.v[0], tri.v[2] - tri.v[0]);
  num_t d1 = std::abs(v3_t<num_t>::dot(n, tri.v[0]));
  num_t d2 = v3_t<num_t>::dot(n.abs(), aabb.extent);
  return d1 <= d2;
}

template<typename num_t>
bool aabb3_t<num_t>::intersects(tri3_t<num_t> tri) const {
  tri.v[0] -= mid;
  tri.v[1] -= mid;
  tri.v[2] -= mid;

  if (!intersectsTriangle_axisTest_2<0>(*this, tri)) return false;
  if (!intersectsTriangle_axisTest_2<1>(*this, tri)) return false;
  if (!intersectsTriangle_axisTest_2<2>(*this, tri)) return false;

  if (!intersectsTriangle_axisTest_3<1,0,0,0>(*this, tri)) return false;
  if (!intersectsTriangle_axisTest_3<1,0,0,1>(*this, tri)) return false;
  if (!intersectsTriangle_axisTest_3<1,0,0,2>(*this, tri)) return false;
                                                          
  if (!intersectsTriangle_axisTest_3<0,1,0,0>(*this, tri)) return false;
  if (!intersectsTriangle_axisTest_3<0,1,0,1>(*this, tri)) return false;
  if (!intersectsTriangle_axisTest_3<0,1,0,2>(*this, tri)) return false;
                                                          
  if (!intersectsTriangle_axisTest_3<0,0,1,0>(*this, tri)) return false;
  if (!intersectsTriangle_axisTest_3<0,0,1,1>(*this, tri)) return false;
  if (!intersectsTriangle_axisTest_3<0,0,1,2>(*this, tri)) return false;

  if (!intersectsTriangle_axisTest_1(*this, tri)) return false;

  return true;
}



typedef aabb3_t<double> aabb3d_t;
typedef aabb3_t<float>  aabb3f_t;
