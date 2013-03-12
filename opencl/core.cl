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

__kernel void grid(
    __global const float *gx,
    __global const float *gy,
    __global const float *gz,
    __global uchar *out,
    size_t y_pitch,
    size_t z_pitch) {
  size_t ix = get_global_id(0);
  size_t iy = get_global_id(1);
  size_t iz = get_global_id(2);

  size_t sx = get_global_size(0);
  size_t sy = get_global_size(1);

  float3 pos = (float3)(gx[ix], gy[iy], gz[iz]);

  out[ix + iy * y_pitch + iz * z_pitch] = inside(pos) ? 1 : 0;
}

__kernel void chop(
    __global const float3 *av,
    __global const float3 *bv,
    __global float3 *cv,
  const size_t n_steps) {
  size_t i = get_global_id(0);

  float3 a = av[i];
  float3 b = bv[i];
  float3 c;

  bool ai = inside(a);
  bool bi = inside(b);

  for (size_t iter = 0; iter < n_steps; ++iter) {
      c = (a+b)/2.0f;
      bool ci = inside(c);
      if (ci == ai) { a = c; } else { b = c; }
  }

  cv[i] = (a+b)/2.0f;
}
