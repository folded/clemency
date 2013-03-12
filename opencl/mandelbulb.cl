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

#if defined(ITERATIONS)
__constant int Iterations = ITERATIONS;
#else
__constant int Iterations = 10;
#endif

__constant float Bailout = 5.0;
__constant float Power = 8.0;
__constant float Epsilon = 1e-5;

float de(float3 q) {
  float3 z = q;

  float dr = 1.0;
  float r = 0.0;
  float theta, phi, zr;

  for (int i = 0; i < Iterations; i++) {
    r = length(z);
    if (r > Bailout) break;

    theta = acos(z.z / r);
    phi = atan2(z.y, z.x);
    dr = pow(r, Power - 1.0f) * Power * dr + 1.0f;

    zr = pow(r, Power);
    theta = theta * Power;
    phi = phi * Power;

    z = zr * (float3)(sin(theta) * cos(phi), sin(phi) * sin(theta), cos(theta));

    z = z + q;
  }

  return 0.5 * log(r) * r / dr;
}

bool inside(float3 pos) {
  return de(pos) < Epsilon;
}

#include <core.cl>
