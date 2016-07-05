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



#include <clemency/geom/v3.hpp>



template<typename num_t>
class sphere_t {
public:
  v3_t<num_t> v;
  num_t r;

  static sphere_t init(const v3_t<num_t> &v, num_t r) {
    sphere_t result;
    result.v = v;
    result.r = r;
    return result;
  }

  bool contains(const v3_t<num_t> &p) const {
    return v3_t<num_t>::sub(p, v).lengthsq() <= r * r;
  }

  bool intersects(const v3_t<num_t> &p) const {
    return contains(p);
  }

  bool contains(const sphere_t &sp) const {
    return v3_t<num_t>::sub(sp.v, v).lengthsq() <= (r - sp.r) * (r - sp.r);
  }

  bool intersects(const sphere_t &sp) const {
    return v3_t<num_t>::sub(sp.v, v).lengthsq() <= (r + sp.r) * (r + sp.r);
  }
};



template<typename rng_t, typename num_t>
v3_t<num_t> randomOnSphere(rng_t &rng, const sphere_t<num_t> &sp) {
  boost::uniform_on_sphere<num_t> distrib(3);
  boost::variate_generator<rng_t &, boost::uniform_on_sphere<num_t> > gen(rng, distrib);
  return v3_t<num_t>::add(v3_t<num_t>::init(gen()).scale(sp.r), sp.v);
}



template<typename rng_t, typename num_t>
v3_t<num_t> randomInSphere(rng_t &rng, const sphere_t<num_t> &sp) {
  boost::uniform_on_sphere<num_t> distrib(3);
  boost::uniform_01<num_t> uniform;
  boost::variate_generator<rng_t &, boost::uniform_on_sphere<num_t> > gen(rng, distrib);
  boost::variate_generator<rng_t &, boost::uniform_01<num_t> > genu(rng, uniform);
  return v3_t<num_t>::add(v3_t<num_t>::init(gen()).scale(sp.r * pow(genu(), 1.0/3.0)), sp.v);
}
