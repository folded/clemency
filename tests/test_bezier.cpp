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

#include <gtest/gtest.h>
#include <clemency/geom.hpp>

TEST(closest_point, closest_point) {
  // come up with more tests. these are kind of trivial.
  qbezier2d_t bezier(v2d_t::init(0.0, 1.0), v2d_t::init(0.0, 0.0),
                     v2d_t::init(1.0, 0.0));
  std::pair<v2d_t, double> result;

  result = bezier.closest_point(v2d_t::init(0.5, 0.5));
  EXPECT_EQ(result.first, v2d_t::init(0.25, 0.25));
  EXPECT_EQ(result.second, ::sqrt(2.0 * .25 * .25));

  result = bezier.closest_point(v2d_t::init(0.0, 0.0));
  EXPECT_EQ(result.first, v2d_t::init(0.25, 0.25));
  EXPECT_EQ(result.second, -::sqrt(2.0 * .25 * .25));

  result = bezier.closest_point(v2d_t::init(2.0, 0.0));
  EXPECT_EQ(result.first, v2d_t::init(1.0, 0.0));
  EXPECT_EQ(result.second, 1.0);

  result = bezier.closest_point(v2d_t::init(0.0, 2.0));
  EXPECT_EQ(result.first, v2d_t::init(0.0, 1.0));
  EXPECT_EQ(result.second, 1.0);
}
