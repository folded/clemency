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

TEST(matrix_det, matrix_det)
{
  m3d_t mat = m3d_t::init(+1, +3, +0,
                          +2, -2, +1,
                          -4, +1, -1);
  EXPECT_EQ(mat.determinant(), -5);
}

TEST(matrix_mdet, matrix_mdet)
{
  m3d_t mat = m3d_t::init(+1, +3, +0,
                          +2, -2, +1,
                          -4, +1, -1);
  EXPECT_EQ(mat.mdet(0,0), 1);
  EXPECT_EQ(mat.mdet(1,0), 2);
  EXPECT_EQ(mat.mdet(2,0), -6);

  EXPECT_EQ(mat.mdet(0,1), -3);
  EXPECT_EQ(mat.mdet(1,1), -1);
  EXPECT_EQ(mat.mdet(2,1), 13);

  EXPECT_EQ(mat.mdet(0,2), 3);
  EXPECT_EQ(mat.mdet(1,2), 1);
  EXPECT_EQ(mat.mdet(2,2), -8);
}

TEST(matrix_inv, matrix_inv)
{
  m3d_t mat = m3d_t::init(+1, +3, +0,
                          +2, -2, +1,
                          -4, +1, -1);
  m3d_t inv;
  double det;
  EXPECT_TRUE(mat.invert(inv, det));
  EXPECT_EQ(det, -5);
  EXPECT_TRUE(m3d_t::ident().eq(mat * inv));
}
