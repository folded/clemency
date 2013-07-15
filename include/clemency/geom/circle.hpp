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



#include <boost/random.hpp>



#include "v2.hpp"



template<typename num_t>
class circle_t {
public:
  v2_t<num_t> v;
  num_t r;

  static circle_t init(const v2_t<num_t> &v, num_t r) {
    circle_t result;
    result.v = v;
    result.r = r;
    return result;
  }

  bool contains(const v2_t<num_t> &p) const {
    return (p - v).lengthsq() <= r * r;
  }

  bool intersects(const v2_t<num_t> &p) const {
    return contains(p);
  }

  bool contains(const circle_t &sp) const {
    return (sp.v - v).lengthsq() <= (r - sp.r) * (r - sp.r);
  }

  bool intersects(const circle_t &sp) const {
    return (sp.v - v).lengthsq() <= (r + sp.r) * (r + sp.r);
  }
};



template<typename rng_t, typename num_t>
v2_t<num_t> randomOnCircle(rng_t &rng, const circle_t<num_t> &sp) {
  boost::uniform_on_sphere<num_t> distrib(2);
  boost::variate_generator<rng_t &, boost::uniform_on_sphere<num_t> > gen(rng, distrib);
  return v2_t<num_t>::add(v2_t<num_t>::init(gen()).scale(sp.r), sp.v);
}



template<typename rng_t, typename num_t>
v2_t<num_t> randomInCircle(rng_t &rng, const circle_t<num_t> &sp) {
  boost::uniform_on_sphere<num_t> distrib(2);
  boost::uniform_01<num_t> uniform;
  boost::variate_generator<rng_t &, boost::uniform_on_sphere<num_t> > gen(rng, distrib);
  boost::variate_generator<rng_t &, boost::uniform_01<num_t> > genu(rng, uniform);
  return v2_t<num_t>::add(v2_t<num_t>::init(gen()).scale(sp.r * pow(genu(), 1.0/2.0)), sp.v);
}
