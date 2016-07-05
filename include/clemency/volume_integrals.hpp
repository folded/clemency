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

#include <clemency/mesh.hpp>

// This code is derived from
// http://www.astro.cornell.edu/~carcich/LRO/feature/test/homegrown/misc/hypspin/volInt.c

template <typename mesh_t>
struct volint_t {
  int A;  // alpha
  int B;  // beta
  int C;  // gamma

  // projection integrals
  double P1, Pa, Pb, Paa, Pab, Pbb, Paaa, Paab, Pabb, Pbbb;

  // face integrals
  double Fa, Fb, Fc, Faa, Fbb, Fcc, Faaa, Fbbb, Fccc, Faab, Fbbc, Fcca;

  // volume integrals
  double T0;
  v3d_t T1, T2, TP;

  static inline double SQR(double x) { return x * x; }
  static inline double CUBE(double x) { return x * x * x; }

  void _compute_projection_integrals(const mesh_t* mesh,
                                     const typename mesh_t::tri_t* t) {
    double a0, a1, da;
    double b0, b1, db;
    double a0_2, a0_3, a0_4, b0_2, b0_3, b0_4;
    double a1_2, a1_3, b1_2, b1_3;
    double C1, Ca, Caa, Caaa, Cb, Cbb, Cbbb;
    double Cab, Kab, Caab, Kaab, Cabb, Kabb;
    int i;

    P1 = Pa = Pb = Paa = Pab = Pbb = Paaa = Paab = Pabb = Pbbb = 0.0;

    for (i = 0; i < 3; i++) {
      a0 = mesh->vertices[t->v[i]].pos.v[A];
      b0 = mesh->vertices[t->v[i]].pos.v[B];
      a1 = mesh->vertices[t->v[(i + 1) % 3]].pos.v[A];
      b1 = mesh->vertices[t->v[(i + 1) % 3]].pos.v[B];
      da = a1 - a0;
      db = b1 - b0;
      a0_2 = a0 * a0;
      a0_3 = a0_2 * a0;
      a0_4 = a0_3 * a0;
      b0_2 = b0 * b0;
      b0_3 = b0_2 * b0;
      b0_4 = b0_3 * b0;
      a1_2 = a1 * a1;
      a1_3 = a1_2 * a1;
      b1_2 = b1 * b1;
      b1_3 = b1_2 * b1;

      C1 = a1 + a0;

      Ca = a1 * C1 + a0_2;
      Caa = a1 * Ca + a0_3;
      Caaa = a1 * Caa + a0_4;

      Cb = b1 * (b1 + b0) + b0_2;
      Cbb = b1 * Cb + b0_3;
      Cbbb = b1 * Cbb + b0_4;

      Cab = 3 * a1_2 + 2 * a1 * a0 + a0_2;
      Kab = a1_2 + 2 * a1 * a0 + 3 * a0_2;

      Caab = a0 * Cab + 4 * a1_3;
      Kaab = a1 * Kab + 4 * a0_3;

      Cabb = 4 * b1_3 + 3 * b1_2 * b0 + 2 * b1 * b0_2 + b0_3;
      Kabb = b1_3 + 2 * b1_2 * b0 + 3 * b1 * b0_2 + 4 * b0_3;

      P1 += db * C1;
      Pa += db * Ca;
      Paa += db * Caa;
      Paaa += db * Caaa;
      Pb += da * Cb;
      Pbb += da * Cbb;
      Pbbb += da * Cbbb;
      Pab += db * (b1 * Cab + b0 * Kab);
      Paab += db * (b1 * Caab + b0 * Kaab);
      Pabb += da * (a1 * Cabb + a0 * Kabb);
    }

    P1 /= 2.0;
    Pa /= 6.0;
    Paa /= 12.0;
    Paaa /= 20.0;
    Pb /= -6.0;
    Pbb /= -12.0;
    Pbbb /= -20.0;
    Pab /= 24.0;
    Paab /= 60.0;
    Pabb /= -60.0;
  }

  void _compute_face_integrals(const mesh_t* mesh,
                               const typename mesh_t::tri_t* t,
                               const v3d_t& norm,
                               const double w) {
    double k1, k2, k3, k4;

    _compute_projection_integrals(mesh, t);

    k1 = 1 / norm.v[C];
    k2 = k1 * k1;
    k3 = k2 * k1;
    k4 = k3 * k1;

    Fa = k1 * Pa;
    Fb = k1 * Pb;
    Fc = -k2 * (norm.v[A] * Pa + norm.v[B] * Pb + w * P1);

    Faa = k1 * Paa;
    Fbb = k1 * Pbb;
    Fcc = k3 * (SQR(norm.v[A]) * Paa + 2 * norm.v[A] * norm.v[B] * Pab +
                SQR(norm.v[B]) * Pbb +
                w * (2 * (norm.v[A] * Pa + norm.v[B] * Pb) + w * P1));

    Faaa = k1 * Paaa;
    Fbbb = k1 * Pbbb;
    Fccc =
        -k4 * (CUBE(norm.v[A]) * Paaa + 3 * SQR(norm.v[A]) * norm.v[B] * Paab +
               3 * norm.v[A] * SQR(norm.v[B]) * Pabb + CUBE(norm.v[B]) * Pbbb +
               3 * w * (SQR(norm.v[A]) * Paa + 2 * norm.v[A] * norm.v[B] * Pab +
                        SQR(norm.v[B]) * Pbb) +
               w * w * (3 * (norm.v[A] * Pa + norm.v[B] * Pb) + w * P1));

    Faab = k1 * Paab;
    Fbbc = -k2 * (norm.v[A] * Pabb + norm.v[B] * Pbbb + w * Pbb);
    Fcca = k3 * (SQR(norm.v[A]) * Paaa + 2 * norm.v[A] * norm.v[B] * Paab +
                 SQR(norm.v[B]) * Pabb +
                 w * (2 * (norm.v[A] * Paa + norm.v[B] * Pab) + w * Pa));
  }

  template <typename iter_t>
  void _compute_volume_integrals(const mesh_t* mesh,
                                 iter_t fbegin,
                                 iter_t fend) {
    T0 = 0.0;
    T1 = v3d_t::zero();
    T2 = v3d_t::zero();
    TP = v3d_t::zero();

    for (; fbegin != fend; ++fbegin) {
      const typename mesh_t::tri_t* t = *fbegin;
      v3d_t norm = mesh->normal(*t);

      C = norm.max_component();
      A = (C + 1) % 3;
      B = (C + 2) % 3;

      double w = -v3d_t::dot(norm, mesh->vertices[t->a].pos);

      _compute_face_integrals(mesh, t, norm, w);

      T0 += norm.x * ((A == 0) ? Fa : ((B == 0) ? Fb : Fc));

      T1[A] += norm.v[A] * Faa;
      T1[B] += norm.v[B] * Fbb;
      T1[C] += norm.v[C] * Fcc;

      T2[A] += norm.v[A] * Faaa;
      T2[B] += norm.v[B] * Fbbb;
      T2[C] += norm.v[C] * Fccc;

      TP[A] += norm.v[A] * Faab;
      TP[B] += norm.v[B] * Fbbc;
      TP[C] += norm.v[C] * Fcca;
    }

    T1.x /= 2;
    T1.y /= 2;
    T1.z /= 2;
    T2.x /= 3;
    T2.y /= 3;
    T2.z /= 3;
    TP.x /= 2;
    TP.y /= 2;
    TP.z /= 2;
  }

  template <typename iter_t>
  volint_t(const mesh_t* mesh, iter_t fbegin, iter_t fend) {
    _compute_volume_integrals(mesh, fbegin, fend);
  }

  double volume() const { return T0; }

  v3d_t centre_of_mass() const { return T1 / T0; }
};
