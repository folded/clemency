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

#include <math.h>
#include <iosfwd>

#include <algorithm>
#include <vector>

#include <tr1/functional>

#include <boost/random.hpp>

class v2_t {
public:
  union {
    struct { double x,y; };
    double v[2];
  };

  static v2_t init(const double *_v) {
    v2_t result;
    std::copy(_v, _v+2, result.v);
    return result;
  }

  static v2_t init(const std::vector<double> &_v) {
    v2_t result;
    std::copy(_v.begin(), _v.begin()+2, result.v);
    return result;
  }

  static v2_t zero() {
    return v2_t::init(0.0, 0.0);
  }

  static v2_t init(double x = 0.0, double y = 0.0) {
    v2_t result;
    result.x = x;
    result.y = y;
    return result;
  }

  double operator[](int i) const {
    return v[i];
  }

  double &operator[](int i) {
    return v[i];
  }

  double length() const {
    return ::sqrt(dot(*this, *this));
  }

  double lengthsq() const {
    return dot(*this, *this);
  }

  v2_t scale(double s) const {
    return init(x * s, y * s);
  }

  v2_t normalize() const {
    double rl = 1.0/length();
    return init(x * rl, y * rl);
  }

  static v2_t add(const v2_t &a, const v2_t &b) {
    return init(a.x + b.x, a.y + b.y);
  }

  static v2_t add(const v2_t &a, const v2_t &b, const v2_t &c) {
    return init(a.x + b.x + c.x, a.y + b.y + c.y);
  }

  v2_t &add(const v2_t &a) {
    x += a.x; y += a.y;
    return *this;
  }

  static v2_t sub(const v2_t &a, const v2_t &b) {
    return init(a.x - b.x, a.y - b.y);
  }

  v2_t &sub(const v2_t &a) {
    x -= a.x; y -= a.y;
    return *this;
  }

  static v2_t neg(const v2_t &a) {
    return init(-a.x, -a.y);
  }

  static double dot(const v2_t &a, const v2_t &b) {
    return a.x * b.x + a.y * b.y;
  }

  static double cross(const v2_t &a, const v2_t &b) {
    return a.x * b.y - b.y * a.x;
  }

  static bool eq(const v2_t &a, const v2_t &b) {
    return std::equal(a.v, a.v+2, b.v);
  }

  static bool neq(const v2_t &a, const v2_t &b) {
    return !eq(a, b);
  }

  static bool lt(const v2_t &a, const v2_t &b) {
    return std::lexicographical_compare(a.v, a.v+2, b.v, b.v+2);
  }

  static v2_t cpow(const v2_t &a, double p) {
    double r2 = dot(a, a);
    double s = ::pow(r2, p / 2.0);
    return init(s * ::cos(p * ::atan2(a.y, a.x)),
                s * ::sin(p * ::atan2(a.y, a.x)));
  }

  static v2_t csqr(const v2_t &a) {
    return init(a.x * a.x - a.y * a.y, 2.0 * a.x * a.y);
  }

  static v2_t cinv(const v2_t &a) {
    double a_dot_a = dot(a, a);
    return init(a.x / a_dot_a, -a.y / a_dot_a);
  }

  static v2_t cmul(const v2_t &a, const v2_t &b) {
	return init(a.x * b.x - a.y * b.y,
                    a.x * b.y + a.y * b.x);
  }
  
  static v2_t cdiv(const v2_t &a, const v2_t &b) {
    return cmul(a, cinv(b));
  }

  /** 
   * \brief Return the orientation of c with respect to the ray defined by a->b.
   * 
   * @param[in] a 
   * @param[in] b 
   * @param[in] c 
   * 
   * @return positive, if c to the left of a->b.
   *         zero,     if c is colinear with a->b.
   *         negative, if c to the right of a->b.
   */
  static double orient(const v2_t &a, const v2_t &b, const v2_t &c) {
    double acx = a.x - c.x;
    double bcx = b.x - c.x;
    double acy = a.y - c.y;
    double bcy = b.y - c.y;
    return acx * bcy - acy * bcx;
  }
};

static inline v2_t operator*(double a, const v2_t &b) {
  return b.scale(a);
}

static  inline v2_t operator/(const v2_t &a, double b) {
  return a.scale(1.0/b);
}

static inline v2_t operator+(const v2_t &a, const v2_t &b) {
  return v2_t::add(a, b);
}

static inline v2_t operator-(const v2_t &a, const v2_t &b) {
  return v2_t::sub(a, b);
}

static inline v2_t &operator+=(v2_t &a, const v2_t &b) {
  return a.add(b);
}

static inline v2_t &operator-=(v2_t &a, const v2_t &b) {
  return a.sub(b);
}

static inline v2_t operator-(const v2_t &a) {
  return v2_t::neg(a);
}

static inline bool operator==(const v2_t &a, const v2_t &b) {
  return v2_t::eq(a, b);
}

static inline bool operator!=(const v2_t &a, const v2_t &b) {
  return !v2_t::eq(a, b);
}

static inline bool operator<(const v2_t &a, const v2_t &b) {
  return v2_t::lt(a, b);
}

static inline bool operator>=(const v2_t &a, const v2_t &b) {
  return !v2_t::lt(a, b);
}

static inline bool operator>(const v2_t &a, const v2_t &b) {
  return v2_t::lt(b, a);
}

static inline bool operator<=(const v2_t &a, const v2_t &b) {
  return !v2_t::lt(b, a);
}



static inline std::ostream &operator<<(std::ostream &out, const v2_t &v) {
  out << "<" << v.x << "," << v.y << ">";
  return out;
}

static inline std::vector<double> solve1(double a, double b, const double zero = 1e-10) {
    std::vector<double> r;

    if (::fabs(a) > zero) {
      r.push_back(-b / a);
    }

    return r;
}

static inline std::vector<double> solve2(double a, double b, double c, const double zero = 1e-10) {
  std::vector<double> r;

  if (::fabs(a) < zero) {
    return solve1(b, c, zero);
  }

  double D = b*b - 4.0*a*c;
  if (D > 0) {
    if (D < zero) {
      r.push_back(-b / (2.0 * a));
    } else {
      r.reserve(2);

      D = ::sqrt(D);
      r.push_back((-b - D) / (2.0 * a));
      r.push_back((-b + D) / (2.0 * a));
    }
  }

  return r;
}

static inline std::vector<double> solve3(double a, double b, double c, double d, const double zero = 1e-10) {
  std::vector<double> r;

  if (::fabs(a) < zero) {
    return solve2(b, c, d, zero);
  }

  b /= a;
  c /= a;
  d /= a;

  double p = c - b*b / 3.0;
  double q = b * (2.0*b*b - 9.0*c) / 27.0 + d;
  double p3 = p*p*p;
  double D = q*q + 4.0*p3 / 27.0;
  double offset = -b / 3.0;

  if (D > zero) {
    double z = ::sqrt(D);
    double u = (-q + z) / 2.0;
    double v = (-q - z) / 2.0;
    u = (u >= 0) ? ::pow(u, 1.0/3.0) : -::pow(-u, 1.0/3.0);
    v = (v >= 0) ? ::pow(v, 1.0/3.0) : -::pow(-v, 1.0/3.0);

    r.push_back(u + v + offset);
  } else if (D < -zero) {
    double u = 2.0 * ::sqrt(-p / 3.0);
    double v = ::acos(-::sqrt(-27.0 / p3) * q / 2.0) / 3.0;

    r.reserve(3);

    r.push_back(u * ::cos(v) + offset);
    r.push_back(u * ::cos(v + 2.0 * M_PI / 3.0) + offset);
    r.push_back(u * ::cos(v + 4.0 * M_PI / 3.0) + offset);
  } else {
    double u = (q < 0) ? ::pow(-q / 2.0, 1.0 / 3.0) : -::pow(q / 2.0, 1.0 / 3.0);

    r.reserve(2);

    r.push_back(2.0*u + offset);
    r.push_back(-u + offset);
  }

  return r;
}



struct qbezier2_t {
  v2_t p0, p1, p2;
  double A, B, C;

  void init(const v2_t &_p0, const v2_t &_p1, const v2_t &_p2) {
    p0 = _p0;
    p1 = _p1 - _p0;
    p2 = _p2 - _p0;

    A = +4.0 * v2_t::dot(p1, p1) - 4.0 * v2_t::dot(p1, p2) + v2_t::dot(p2, p2);
    B = -6.0 * v2_t::dot(p1, p1) + 3.0 * v2_t::dot(p1, p2);
    C = +2.0 * v2_t::dot(p1, p1);
  }

  qbezier2_t() {
    init(v2_t::zero(), v2_t::zero(), v2_t::zero());
  }

  qbezier2_t(const v2_t &_p0, const v2_t &_p1, const v2_t &_p2) {
    init(_p0, _p1, _p2);
  }

  v2_t get_relpos(double t) const {
    return 2.0*t*(1.0-t) * p1 + t*t * p2;
  }

  v2_t get_pos(double t) const {
    return p0 + get_relpos(t);
  }

  v2_t get_tangent(double t) const {
    return 2.0*(1-t) * p1 + 2.0*t * (p2-p1);
  }

  std::pair<v2_t, double> closest_point(v2_t pos) const {
    v2_t p = pos - p0;
    double p1p = v2_t::dot(p1, p);
    double p2p = v2_t::dot(p2, p);

    std::vector<double> sol = solve3(A,
                                     B,
                                     C - p2p + 2.0 * p1p,
                                     -p1p);
    double d;
    double min_dist = p.lengthsq();
    v2_t min_pos = v2_t::zero();
    double min_t = 0.0;

    d = (p - p2).lengthsq();
    if (d < min_dist) { min_dist = d; min_pos = p2; min_t = 1.0; }

    for (size_t i = 0; i < sol.size(); ++i) {
      if (sol[i] > 0.0 && sol[i] < 1.0) {
        v2_t q = get_relpos(sol[i]);
        d = (p - q).lengthsq();
        if (d < min_dist) { min_dist = d; min_pos = q; min_t = sol[i]; }
      }
    }

    min_dist = ::sqrt(min_dist);
    if (v2_t::orient(min_pos, min_pos + get_tangent(min_t), p) < 0.0) {
      min_dist = -min_dist;
    }

    return std::make_pair(min_pos + p0, min_dist);
  }
};



class v3_t {
public:
  union {
    struct { double x,y,z; };
    double v[3];
  };

  static v3_t init(const double *_v) {
    v3_t result;
    std::copy(_v, _v+3, result.v);
    return result;
  }

  static v3_t init(const std::vector<double> &_v) {
    v3_t result;
    std::copy(_v.begin(), _v.begin()+3, result.v);
    return result;
  }

  static v3_t zero() {
    return v3_t::init(0.0, 0.0, 0.0);
  }

  static v3_t init(double x = 0.0, double y = 0.0, double z = 0.0) {
    v3_t result;
    result.x = x;
    result.y = y;
    result.z = z;
    return result;
  }

  double operator[](int i) const {
    return v[i];
  }

  double &operator[](int i) {
    return v[i];
  }

  double length() const {
    return ::sqrt(dot(*this, *this));
  }

  double lengthsq() const {
    return dot(*this, *this);
  }

  v3_t scale(double s) const {
    return init(x * s, y * s, z * s);
  }

  v3_t normalize() const {
    double rl = 1.0/length();
    return init(x * rl, y * rl, z * rl);
  }

  static v3_t add(const v3_t &a, const v3_t &b) {
    return init(a.x + b.x, a.y + b.y, a.z + b.z);
  }

  static v3_t add(const v3_t &a, const v3_t &b, const v3_t &c) {
    return init(a.x + b.x + c.x, a.y + b.y + c.y, a.z + b.z + c.z);
  }

  v3_t &add(const v3_t &a) {
    x += a.x; y += a.y; z += a.z;
    return *this;
  }

  static v3_t sub(const v3_t &a, const v3_t &b) {
    return init(a.x - b.x, a.y - b.y, a.z - b.z);
  }

  v3_t &sub(const v3_t &a) {
    x -= a.x; y -= a.y; z -= a.z;
    return *this;
  }

  static v3_t neg(const v3_t &a) {
    return init(-a.x, -a.y, -a.z);
  }

  static double dot(const v3_t &a, const v3_t &b) {
    return a.x * b.x + a.y * b.y + a.z * b.z;
  }

  static v3_t cross(const v3_t &a, const v3_t &b) {
    return init(+(a.y * b.z - b.y * a.z),
                -(a.x * b.z - b.x * a.z),
                +(a.x * b.y - b.x * a.y));
  }

  static bool eq(const v3_t &a, const v3_t &b) {
    return std::equal(a.v, a.v+3, b.v);
  }

  static bool neq(const v3_t &a, const v3_t &b) {
    return !eq(a, b);
  }

  static bool lt(const v3_t &a, const v3_t &b) {
    return std::lexicographical_compare(a.v, a.v+3, b.v, b.v+3);
  }
};

static inline v3_t operator*(double a, const v3_t &b) {
  return b.scale(a);
}

static  inline v3_t operator/(const v3_t &a, double b) {
  return a.scale(1.0/b);
}

static inline v3_t operator+(const v3_t &a, const v3_t &b) {
  return v3_t::add(a, b);
}

static inline v3_t operator-(const v3_t &a, const v3_t &b) {
  return v3_t::sub(a, b);
}

static inline v3_t &operator+=(v3_t &a, const v3_t &b) {
  return a.add(b);
}

static inline v3_t &operator-=(v3_t &a, const v3_t &b) {
  return a.sub(b);
}

static inline v3_t operator-(const v3_t &a) {
  return v3_t::neg(a);
}

static inline bool operator==(const v3_t &a, const v3_t &b) {
  return v3_t::eq(a, b);
}

static inline bool operator!=(const v3_t &a, const v3_t &b) {
  return !v3_t::eq(a, b);
}

static inline bool operator<(const v3_t &a, const v3_t &b) {
  return v3_t::lt(a, b);
}

static inline bool operator>=(const v3_t &a, const v3_t &b) {
  return !v3_t::lt(a, b);
}

static inline bool operator>(const v3_t &a, const v3_t &b) {
  return v3_t::lt(b, a);
}

static inline bool operator<=(const v3_t &a, const v3_t &b) {
  return !v3_t::lt(b, a);
}



static inline std::ostream &operator<<(std::ostream &out, const v3_t &v) {
  out << "<" << v.x << "," << v.y << "," << v.z << ">";
  return out;
}



struct v3_hash_t {
  size_t operator()(const v3_t &a) const {
    std::tr1::hash<double> h;
    size_t r = 0;
    r *= 131; r ^= h(a.x);
    r *= 131; r ^= h(a.y);
    r *= 131; r ^= h(a.z);
    return r;
  }

  size_t operator()(const std::pair<v3_t, v3_t> &a) const {
    std::tr1::hash<double> h;
    size_t r = 0;
    r *= 131; r ^= h(a.first.x);
    r *= 131; r ^= h(a.first.y);
    r *= 131; r ^= h(a.first.z);
    r *= 131; r ^= h(a.second.x);
    r *= 131; r ^= h(a.second.y);
    r *= 131; r ^= h(a.second.z);
    return r;
  }
};



class quat_t {
public:
  double a;
  v3_t b;

  static quat_t zero() {
    return quat_t::init(0.0, v3_t::zero());
  }

  static quat_t init(double a, const v3_t &b) {
    quat_t result;
    result.a = a;
    result.b = b;
    return result;
  }

  static quat_t rot(double angle, const v3_t &axis) {
    quat_t r;
    double len = axis.length();
    if (fabs(len) > 1e-10) {
      double omega = -0.5 * angle;
      double c = cos(omega);
      double s = sin(omega);

      r.a = c;
      r.b.x = axis.x * s / len;
      r.b.y = axis.y * s / len;
      r.b.z = axis.z * s / len;

      r.normalize();
    } else {
      r.a = 1.0;
      r.b.x = r.b.y = r.b.z = 0.0;
    }
    return r;
  }

  quat_t sq() const {
    return init(a * a - b.lengthsq(), b.scale(2.0 * a));
  }

  double length() const {
    return ::sqrt(quat_t::dot(*this, *this));
  }

  double lengthsq() const {
    return quat_t::dot(*this, *this);
  }

  quat_t scale(double s) const {
    return init(a * s, b.scale(s));
  }

  quat_t normalize() const {
    return scale(1.0/length());
  }

  static quat_t add(const quat_t &a, const quat_t &b) {
    return init(a.a + b.a, v3_t::add(a.b, b.b));
  }

  static quat_t sub(const quat_t &a, const quat_t &b) {
    return init(a.a - b.a, v3_t::sub(a.b, b.b));
  }

  static quat_t add(const quat_t &a, const quat_t &b, const quat_t &c) {
    return init(a.a + b.a + c.a, v3_t::add(a.b, b.b, c.b));
  }

  static double dot(const quat_t &a, const quat_t &b) {
    return a.a * b.a + v3_t::dot(a.b, b.b);
  }

  static quat_t mul(const quat_t &a, const quat_t &b) {
    return init(a.a * b.a - v3_t::dot(a.b, b.b),
                v3_t::add(b.b.scale(a.a),
                        a.b.scale(b.a),
                        v3_t::cross(a.b, b.b)));
  }
};

static inline quat_t operator*(const quat_t &a, const quat_t &b) {
  return quat_t::mul(a, b);
}

static  inline quat_t operator*(double a, const quat_t &b) {
  return b.scale(a);
}

static inline quat_t operator+(const quat_t &a, const quat_t &b) {
  return quat_t::add(a, b);
}

static inline quat_t operator-(const quat_t &a, const quat_t &b) {
  return quat_t::sub(a, b);
}



static inline std::ostream &operator<<(std::ostream &out, const quat_t &q) {
  out << "<" << q.a << ";" << q.b.x << "," << q.b.y << "," << q.b.z << ">";
  return out;
}



class m3_t {
public:
  union {
    // .mCR
    struct {
      double
        m00, m01, m02,
        m10, m11, m12,
        m20, m21, m22;
    };
    // m[C][R]
    double m[3][3];
    // v[C*3+R]
    double v[9];
  };

  static m3_t ident() {
    m3_t r;
    r.m00 = 1.0; r.m01 = 0.0; r.m02 = 0.0;
    r.m10 = 0.0; r.m11 = 1.0; r.m12 = 0.0;
    r.m20 = 0.0; r.m21 = 0.0; r.m22 = 1.0;
    return r;
  }

  static m3_t init(double m00, double m10, double m20,
                   double m01, double m11, double m21,
                   double m02, double m12, double m22) {
    m3_t r;
    r.m00 = m00; r.m01 = m01; r.m02 = m02;
    r.m10 = m10; r.m11 = m11; r.m12 = m12;
    r.m20 = m20; r.m21 = m21; r.m22 = m22;
    return r;
  }

  static m3_t rot(const quat_t &q) {
    double w = q.a;
    double x = q.b.x;
    double y = q.b.y;
    double z = q.b.z;

    return m3_t::init(1 - 2*y*y - 2*z*z,   2*x*y - 2*z*w,       2*x*z + 2*y*w,
                      2*x*y + 2*z*w,       1 - 2*x*x - 2*z*z,   2*y*z - 2*x*w,
                      2*x*z - 2*y*w,       2*y*z + 2*x*w,       1 - 2*x*x - 2*y*y);
  }

  static m3_t scale(const v3_t &v) {
    m3_t r;
    r.m00 = v.x; r.m01 = 0.0; r.m02 = 0.0;
    r.m10 = 0.0; r.m11 = v.y; r.m12 = 0.0;
    r.m20 = 0.0; r.m21 = 0.0; r.m22 = v.z;
    return r;
  }

  m3_t transpose() const {
    m3_t r;
    r.m00 = m00; r.m01 = m10; r.m02 = m20;
    r.m10 = m01; r.m11 = m11; r.m12 = m21;
    r.m20 = m02; r.m21 = m12; r.m22 = m22;
    return r;
  }

  m3_t neg() const {
    m3_t r;
    for (size_t i = 0; i < 9; ++i) r.v[i] = -v[i];
    return r;
  }

  double operator()(int i, int j) const {
    return m[i][j];
  }

  double &operator()(int i, int j) {
    return m[i][j];
  }

  double mdet(int i, int j) const {
    return
      + m[i==0?1:0][j==0?1:0] * m[i==2?1:2][j==2?1:2]
      - m[i==2?1:2][j==0?1:0] * m[i==0?1:0][j==2?1:2];
  }

  double determinant() const {
    return
      +m00 * (m11 * m22 - m21 * m12)
      -m10 * (m01 * m22 - m21 * m02)
      +m20 * (m01 * m12 - m11 * m02);
  }

  bool invert(m3_t &inv, double &det) const {
    det = determinant();
    if (fabs(det) < 1e-10) return false;
    double d = 1.0/det;

    inv.m00 = +d*mdet(0,0); inv.m01 = -d*mdet(1,0); inv.m02 = +d*mdet(2,0);
    inv.m10 = -d*mdet(0,1); inv.m11 = +d*mdet(1,1); inv.m12 = -d*mdet(2,1);
    inv.m20 = +d*mdet(0,2); inv.m21 = -d*mdet(1,2); inv.m22 = +d*mdet(2,2);

    return true;
  }

  v3_t transform(const v3_t &v) const {
    return v3_t::init(v.x * m00 + v.y * m10 + v.z * m20,
                      v.x * m01 + v.y * m11 + v.z * m21,
                      v.x * m02 + v.y * m12 + v.z * m22);
  }

  bool exactly_eq(const m3_t &o) const {
    for (size_t i = 0; i < 9; ++i) if (v[i] != o.v[i]) return false;
    return true;
  }

  bool eq(const m3_t &o, double eps = 1e-10) const {
    for (size_t i = 0; i < 9; ++i) if (fabs(v[i] - o.v[i]) > eps) return false;
    return true;
  }

  static m3_t sub(const m3_t &a, const m3_t &b) {
    m3_t r;
    for (size_t i = 0; i < 9; ++i) r.v[i] = a.v[i] - b.v[i];
    return r;
  }

  static m3_t add(const m3_t &a, const m3_t &b) {
    m3_t r;
    for (size_t i = 0; i < 9; ++i) r.v[i] = a.v[i] + b.v[i];
    return r;
  }

  static m3_t mul(const m3_t &a, const m3_t &b) {
    m3_t r;
    for (size_t i = 0; i < 3; ++i) {
      for (size_t j = 0; j < 3; ++j) {
        r.m[i][j] =
          a.m[i][0] * b.m[0][j] +
          a.m[i][1] * b.m[1][j] +
          a.m[i][2] * b.m[2][j];
      }
    }
    return r;
  }
};



static inline bool operator==(const m3_t &a, const m3_t &b) {
  return a.exactly_eq(b);
}

static inline m3_t operator-(const m3_t &a, const m3_t &b) {
  return m3_t::sub(a, b);
}

static inline m3_t operator+(const m3_t &a, const m3_t &b) {
  return m3_t::add(a, b);
}

static inline m3_t operator-(const m3_t &a) {
  return a.neg();
}

static inline m3_t operator*(const m3_t &a, const m3_t &b) {
  return m3_t::mul(a, b);
}



static inline std::ostream &operator<<(std::ostream &o, const m3_t &m) {
  for (size_t i = 0; i < 3; ++i) {
    o << "[";
    for (size_t j = 0; j < 3; ++j) {
      o << " " << m.m[j][i];
    }
    o << " ]\n";
  }
  return o;
}



class m4_t {
public:
  union {
    // column major storage.
    // .mCR
    struct {
      double
        m00, m01, m02, m03,
        m10, m11, m12, m13,
        m20, m21, m22, m23,
        m30, m31, m32, m33;
    };
    // m[C][R]
    double m[4][4];
    // v[C*4+R]
    double v[16];
  };

  static m4_t ident() {
    m4_t r;
    r.m00 = 1.0; r.m01 = 0.0; r.m02 = 0.0; r.m03 = 0.0;
    r.m10 = 0.0; r.m11 = 1.0; r.m12 = 0.0; r.m13 = 0.0;
    r.m20 = 0.0; r.m21 = 0.0; r.m22 = 1.0; r.m23 = 0.0;
    r.m30 = 0.0; r.m31 = 0.0; r.m32 = 0.0; r.m33 = 1.0;
    return r;
  }

  // note: arguments are provided in the natural orientation, but stored transposed.
  static m4_t init(double m00, double m10, double m20, double m30,
                   double m01, double m11, double m21, double m31,
                   double m02, double m12, double m22, double m32,
                   double m03, double m13, double m23, double m33) {
    m4_t r;
    r.m00 = m00; r.m01 = m01; r.m02 = m02; r.m03 = m03;
    r.m10 = m10; r.m11 = m11; r.m12 = m12; r.m13 = m13;
    r.m20 = m20; r.m21 = m21; r.m22 = m22; r.m23 = m23;
    r.m30 = m30; r.m31 = m31; r.m32 = m32; r.m33 = m33;
    return r;
  }

  static m4_t rot(const quat_t &q) {
    double w = q.a;
    double x = q.b.x;
    double y = q.b.y;
    double z = q.b.z;

    return m4_t::init(1 - 2*y*y - 2*z*z,   2*x*y - 2*z*w,       2*x*z + 2*y*w,       0.0,  
                      2*x*y + 2*z*w,       1 - 2*x*x - 2*z*z,   2*y*z - 2*x*w,       0.0,  
                      2*x*z - 2*y*w,       2*y*z + 2*x*w,       1 - 2*x*x - 2*y*y,   0.0,  
                      0.0,                 0.0,                 0.0,                 1.0);
  }

  static m4_t scale(const v3_t &v) {
    m4_t r;
    r.m00 = v.x; r.m01 = 0.0; r.m02 = 0.0; r.m03 = 0.0;
    r.m10 = 0.0; r.m11 = v.y; r.m12 = 0.0; r.m13 = 0.0;
    r.m20 = 0.0; r.m21 = 0.0; r.m22 = v.z; r.m23 = 0.0;
    r.m30 = 0.0; r.m31 = 0.0; r.m32 = 0.0; r.m33 = 1.0;
    return r;
  }

  static m4_t translate(const v3_t &v) {
    m4_t r;
    r.m00 = 1.0; r.m01 = 0.0; r.m02 = 0.0; r.m03 = 0.0;
    r.m10 = 0.0; r.m11 = 1.0; r.m12 = 0.0; r.m13 = 0.0;
    r.m20 = 0.0; r.m21 = 0.0; r.m22 = 1.0; r.m23 = 0.0;
    r.m30 = v.x; r.m31 = v.y; r.m32 = v.z; r.m33 = 1.0;
    return r;
  }

  m4_t transpose() const {
    m4_t r;
    r.m00 = m00; r.m01 = m10; r.m02 = m20; r.m03 = m30;
    r.m10 = m01; r.m11 = m11; r.m12 = m21; r.m13 = m31;
    r.m20 = m02; r.m21 = m12; r.m22 = m22; r.m23 = m32;
    r.m30 = m03; r.m31 = m13; r.m32 = m23; r.m33 = m33;
    return r;
  }

  static m4_t mul(const m4_t &a, const m4_t &b) {
    m4_t r;
    for (size_t i = 0; i < 4; ++i) {
      for (size_t j = 0; j < 4; ++j) {
        r.m[i][j] = 0.0;
        for (size_t k = 0; k < 4; ++k) {
          r.m[i][j] += a.m[k][j] * b.m[i][k];
        }
      }
    }
    return r;
  }

  double operator()(int col, int row) const {
    return m[col][row];
  }

  double &operator()(int col, int row) {
    return m[col][row];
  }

  v3_t transform(const v3_t &v) const {
    return v3_t::init(v.x * m00 + v.y * m10 + v.z * m20 + m30,
                      v.x * m01 + v.y * m11 + v.z * m21 + m31,
                      v.x * m02 + v.y * m12 + v.z * m22 + m32);
  }
};



static inline m4_t operator*(const m4_t &a, const m4_t &b) {
  return m4_t::mul(a, b);
}



class sphere_t {
public:
  v3_t v;
  double r;

  static sphere_t init(const v3_t &v, double r) {
    sphere_t result;
    result.v = v;
    result.r = r;
    return result;
  }

  bool contains(const v3_t &p) const {
    return v3_t::sub(p, v).lengthsq() <= r * r;
  }

  bool intersects(const v3_t &p) const {
    return contains(p);
  }

  bool contains(const sphere_t &sp) const {
    return v3_t::sub(sp.v, v).lengthsq() <= (r - sp.r) * (r - sp.r);
  }

  bool intersects(const sphere_t &sp) const {
    return v3_t::sub(sp.v, v).lengthsq() <= (r + sp.r) * (r + sp.r);
  }
};



template<typename rng_t>
v3_t randomOnSphere(rng_t &rng, const sphere_t &sp) {
  boost::uniform_on_sphere<double> distrib(3);
  boost::variate_generator<rng_t &, boost::uniform_on_sphere<double> > gen(rng, distrib);
  return v3_t::add(v3_t::init(gen()).scale(sp.r), sp.v);
}



template<typename rng_t>
v3_t randomInSphere(rng_t &rng, const sphere_t &sp) {
  boost::uniform_on_sphere<double> distrib(3);
  boost::uniform_01<double> uniform;
  boost::variate_generator<rng_t &, boost::uniform_on_sphere<double> > gen(rng, distrib);
  boost::variate_generator<rng_t &, boost::uniform_01<double> > genu(rng, uniform);
  return v3_t::add(v3_t::init(gen()).scale(sp.r * pow(genu(), 1.0/3.0)), sp.v);
}



class aabb_t {
public:
  v3_t mid, extent;

  static aabb_t init(const v3_t &mid, const v3_t &extent) {
    aabb_t result;
    result.mid = mid;
    result.extent = extent;
    return result;
  }

  static aabb_t initWithBounds(const v3_t &lo, const v3_t &hi) {
    aabb_t result;
    result.mid =
      v3_t::init((hi.x + lo.x) / 2.0,
               (hi.y + lo.y) / 2.0,
               (hi.z + lo.z) / 2.0);
    result.extent =
      v3_t::init(std::max(hi.x - result.mid.x, result.mid.x - lo.x),
               std::max(hi.y - result.mid.y, result.mid.y - lo.y),
               std::max(hi.z - result.mid.z, result.mid.z - lo.z));
    return result;
  }

  double xl() const { return mid.x - extent.x; }
  double xm() const { return mid.x;            }
  double xh() const { return mid.x + extent.x; }

  double yl() const { return mid.y - extent.y; }
  double ym() const { return mid.y;            }
  double yh() const { return mid.y + extent.y; }

  double zl() const { return mid.z - extent.z; }
  double zm() const { return mid.z;            }
  double zh() const { return mid.z + extent.z; }

  v3_t corner(unsigned c) const {
    return v3_t::init(mid.x + extent.x * (c & 1 ? +1.0 : -1.0),
                      mid.y + extent.y * (c & 2 ? +1.0 : -1.0),
                      mid.z + extent.z * (c & 4 ? +1.0 : -1.0));
  }


  bool intersects(const aabb_t &bbox) const {
    return abs(mid.x - bbox.mid.x) < extent.x + bbox.extent.x ||
           abs(mid.y - bbox.mid.y) < extent.y + bbox.extent.y ||
           abs(mid.z - bbox.mid.z) < extent.z + bbox.extent.z;
  }

  bool touchesFace(const aabb_t &bbox) const {
    double dx = abs(mid.x - bbox.mid.x) - extent.x + bbox.extent.x;
    double dy = abs(mid.y - bbox.mid.y) - extent.y + bbox.extent.y;
    double dz = abs(mid.z - bbox.mid.z) - extent.z + bbox.extent.z;
    return
      (dx == 0.0 && dy < 0.0 && dz < 0.0) ||
      (dy == 0.0 && dx < 0.0 && dz < 0.0) ||
      (dz == 0.0 && dx < 0.0 && dy < 0.0);
  }

  bool contains(const v3_t &p) const {
    return
      abs(p.x - mid.x) < extent.x &&
      abs(p.y - mid.y) < extent.y &&
      abs(p.z - mid.z) < extent.z;
  }

  bool intersects(const v3_t &p) const {
    return contains(p);
  }

  bool contains(const sphere_t &sp) const {
    return
      abs(sp.v.x - mid.x) < extent.x - sp.r &&
      abs(sp.v.y - mid.y) < extent.y - sp.r &&
      abs(sp.v.z - mid.z) < extent.z - sp.r;
  }

  bool intersects(const sphere_t &sp) const {
    double r = 0.0;
    for (unsigned i = 0; i < 3; ++i) {
      double t = fabs(sp.v.v[i] - mid.v[i]) - extent.v[i]; if (t > 0.0) r += t*t;
    }
    return r <= sp.r * sp.r;
  }
};



template<typename rng_t>
std::pair<v3_t, v3_t> randomLineSegmentInAABB(rng_t &rng, const aabb_t &aabb) {
  std::pair<v3_t, v3_t> result;
  result.first = randomOnAABB(rng, aabb);
  result.second = randomOnAABB(rng, aabb);
  return result;
}



template<typename rng_t>
v3_t randomOnAABB(rng_t &rng, const aabb_t &aabb) {
  boost::random::uniform_01<double> uniform;
  boost::random::uniform_int_distribution<> face(0, 5);
  boost::random::variate_generator<rng_t &, boost::uniform_01<double> > gen(rng, uniform);
  v3_t pos;
  int f = face(rng);

  if (f != 0 && f != 3) {
    pos.x = aabb.mid.x + (-1.0 + gen() * 2.0) * aabb.extent.x;
  } else if (f == 0) {
    pos.x = aabb.mid.x - aabb.extent.x;
  } else if (f == 3) {
    pos.x = aabb.mid.x + aabb.extent.x;
  }

  if (f != 1 && f != 4) {
    pos.y = aabb.mid.y + (-1.0 + gen() * 2.0) * aabb.extent.y;
  } else if (f == 1) {
    pos.y = aabb.mid.y - aabb.extent.y;
  } else if (f == 4) {
    pos.y = aabb.mid.y + aabb.extent.y;
  }

  if (f != 2 && f != 5) {
    pos.z = aabb.mid.z + (-1.0 + gen() * 2.0) * aabb.extent.z;
  } else if (f == 2) {
    pos.z = aabb.mid.z - aabb.extent.z;
  } else if (f == 5) {
    pos.z = aabb.mid.z + aabb.extent.z;
  }

  return pos;
}



template<typename rng_t>
v3_t randomInAABB(rng_t &rng, const aabb_t &aabb) {
  boost::uniform_01<double> uniform;
  boost::variate_generator<rng_t &, boost::uniform_01<double> > gen(rng, uniform);
  return v3_t::init(aabb.mid.x + (gen() - .5) * aabb.extent.x,
                    aabb.mid.y + (gen() - .5) * aabb.extent.y,
                    aabb.mid.z + (gen() - .5) * aabb.extent.z);
}



template<typename child_t>
class octree_t {
public:
  aabb_t bbox;
  child_t *children[8];

  ~octree_t() {
    clear_children();
  };

  octree_t(const aabb_t &_bbox) : bbox(_bbox) {
    std::fill(children, children+8, (child_t *)NULL);
  }

  void clear_children() {
    if (!is_leaf()) {
      for (size_t i = 0; i < 8; ++i) {
        delete children[i];
        children[i] = NULL;
      }
    }
  }

  bool is_leaf() const {
    return children[0] == NULL;
  }

  v3_t corner(unsigned c) const {
    return bbox.corner(c);
  }

  void split(bool (*filter)(int lev, child_t *), int curr_lev = 0) {
    if (!filter(curr_lev, this)) return;
    if (is_leaf()) split();
    for (size_t i = 0; i < 8; ++i) children[i]->split(filter, curr_lev + 1);
  }

  void split(int n_levels) {
    if (n_levels <= 0) return;
    if (is_leaf()) split();
    for (size_t i = 0; i < 8; ++i) children[i]->split(n_levels-1);
  }

  void split() {
    v3_t lo = v3_t::sub(bbox.mid, bbox.extent);
    v3_t hi = v3_t::add(bbox.mid, bbox.extent);
    v3_t md = bbox.mid;

    children[0] = new child_t(aabb_t::initWithBounds(v3_t::init(lo.x, lo.y, lo.z), v3_t::init(md.x, md.y, md.z)));
    children[1] = new child_t(aabb_t::initWithBounds(v3_t::init(md.x, lo.y, lo.z), v3_t::init(hi.x, md.y, md.z)));
    children[2] = new child_t(aabb_t::initWithBounds(v3_t::init(lo.x, md.y, lo.z), v3_t::init(md.x, hi.y, md.z)));
    children[3] = new child_t(aabb_t::initWithBounds(v3_t::init(md.x, md.y, lo.z), v3_t::init(hi.x, hi.y, md.z)));

    children[4] = new child_t(aabb_t::initWithBounds(v3_t::init(lo.x, lo.y, md.z), v3_t::init(md.x, md.y, hi.z)));
    children[5] = new child_t(aabb_t::initWithBounds(v3_t::init(md.x, lo.y, md.z), v3_t::init(hi.x, md.y, hi.z)));
    children[6] = new child_t(aabb_t::initWithBounds(v3_t::init(lo.x, md.y, md.z), v3_t::init(md.x, hi.y, hi.z)));
    children[7] = new child_t(aabb_t::initWithBounds(v3_t::init(md.x, md.y, md.z), v3_t::init(hi.x, hi.y, hi.z)));
  }

  struct visitor_t {
    virtual void pre(child_t *node, int depth) {}
    virtual void post(child_t *node, int depth) {}
  };

  struct const_visitor_t {
    virtual void pre(const child_t *node, int depth) {}
    virtual void post(const child_t *node, int depth) {}
  };

  void visit(visitor_t &visitor, int depth = 0) {
    visitor.pre(static_cast<child_t *>(this), depth);
    if (!is_leaf()) {
      for (size_t i = 0; i < 8; ++i) {
        children[i]->visit(visitor, depth + 1);
      }
    }
    visitor.post(static_cast<child_t *>(this), depth);
  }

  void visit(const_visitor_t &visitor, int depth = 0) const {
    visitor.pre(static_cast<const child_t *>(this), depth);
    if (!is_leaf()) {
      for (size_t i = 0; i < 8; ++i) {
        children[i]->visit(visitor, depth + 1);
      }
    }
    visitor.post(static_cast<const child_t *>(this), depth);
  }

  template<typename visitor_t>
  void visit_preorder(visitor_t &visitor, int depth = 0) {
    visitor(static_cast<child_t *>(this), depth);
    if (!is_leaf()) {
      for (size_t i = 0; i < 8; ++i) {
        children[i]->visit_preorder(visitor, depth + 1);
      }
    }
  }

  template<typename visitor_t>
  void visit_postorder(visitor_t &visitor, int depth = 0) {
    if (!is_leaf()) {
      for (size_t i = 0; i < 8; ++i) {
        children[i]->visit_postorder(visitor, depth + 1);
      }
    }
    visitor(static_cast<child_t *>(this), depth);
  }

  template<typename visitor_t>
  void visit_preorder(visitor_t &visitor, int depth = 0) const {
    visitor(static_cast<const child_t *>(this), depth);
    if (!is_leaf()) {
      for (size_t i = 0; i < 8; ++i) {
        children[i]->visit_preorder(visitor, depth + 1);
      }
    }
  }

  template<typename visitor_t>
  void visit_postorder(visitor_t &visitor, int depth = 0) const {
    if (!is_leaf()) {
      for (size_t i = 0; i < 8; ++i) {
        children[i]->visit_postorder(visitor, depth + 1);
      }
    }
    visitor(static_cast<const child_t *>(this), depth);
  }
};



template<typename node_t>
class octree_adjacent_cells_t {
  inline bool is_leaf(node_t *node) const {
    return node == NULL || node->is_leaf();
  }

  inline node_t *ch(node_t *node, size_t i) const {
    return is_leaf(node) ? node : node->children[i];
  }

  template<typename iter_t>
  void emit(node_t *q0,
            iter_t &out) const {
    std::vector<node_t *> r;
    r.resize(1);
    r[0] = q0;
    *out++ = r;
  }

  template<typename iter_t>
  void emit(node_t *q0, node_t *q1,
            iter_t &out) const {
    std::vector<node_t *> r;
    r.resize(2);
    r[0] = q0; r[1] = q1;
    *out++ = r;
  }

  template<typename iter_t>
  void emit(node_t *q0, node_t *q1, node_t *q2, node_t *q3,
            iter_t &out) const {
    std::vector<node_t *> r;
    r.resize(4);
    r[0] = q0; r[1] = q1; r[2] = q2; r[3] = q3;
    *out++ = r;
  }

  template<typename iter_t>
  void emit(node_t *q0, node_t *q1, node_t *q2, node_t *q3,
            node_t *q4, node_t *q5, node_t *q6, node_t *q7,
            iter_t &out) const {
    std::vector<node_t *> r;
    r.resize(8);
    r[0] = q0; r[1] = q1; r[2] = q2; r[3] = q3;
    r[4] = q4; r[5] = q5; r[6] = q6; r[7] = q7;
    *out++ = r;
  }

  // 1 arg - volume
  template<int lev, typename iter_t>
  void dual_3(node_t *q0, iter_t &out) const {
    if (is_leaf(q0)) {
      if (lev == 3) emit(q0, out);
      return;
    }

    dual_3<lev>(ch(q0, 0), out);
    dual_3<lev>(ch(q0, 1), out);
    dual_3<lev>(ch(q0, 2), out);
    dual_3<lev>(ch(q0, 3), out);
    dual_3<lev>(ch(q0, 4), out);
    dual_3<lev>(ch(q0, 5), out);
    dual_3<lev>(ch(q0, 6), out);
    dual_3<lev>(ch(q0, 7), out);

    if (lev <= 2) {
      dual_2x<lev>(ch(q0, 0), ch(q0, 1), out);
      dual_2x<lev>(ch(q0, 2), ch(q0, 3), out);
      dual_2x<lev>(ch(q0, 4), ch(q0, 5), out);
      dual_2x<lev>(ch(q0, 6), ch(q0, 7), out);

      dual_2y<lev>(ch(q0, 0), ch(q0, 2), out);
      dual_2y<lev>(ch(q0, 1), ch(q0, 3), out);
      dual_2y<lev>(ch(q0, 4), ch(q0, 6), out);
      dual_2y<lev>(ch(q0, 5), ch(q0, 7), out);

      dual_2z<lev>(ch(q0, 0), ch(q0, 4), out);
      dual_2z<lev>(ch(q0, 1), ch(q0, 5), out);
      dual_2z<lev>(ch(q0, 2), ch(q0, 6), out);
      dual_2z<lev>(ch(q0, 3), ch(q0, 7), out);
    }

    if (lev <= 1) {
      dual_1x<lev>(ch(q0, 0), ch(q0, 2), ch(q0, 4), ch(q0, 6), out);
      dual_1x<lev>(ch(q0, 1), ch(q0, 3), ch(q0, 5), ch(q0, 7), out);

      dual_1y<lev>(ch(q0, 0), ch(q0, 1), ch(q0, 4), ch(q0, 5), out);
      dual_1y<lev>(ch(q0, 2), ch(q0, 3), ch(q0, 6), ch(q0, 7), out);

      dual_1z<lev>(ch(q0, 0), ch(q0, 1), ch(q0, 2), ch(q0, 3), out);
      dual_1z<lev>(ch(q0, 4), ch(q0, 5), ch(q0, 6), ch(q0, 7), out);
    }

    if (lev == 0) {
      dual_0<lev>(ch(q0, 0), ch(q0, 1), ch(q0, 2), ch(q0, 3), ch(q0, 4), ch(q0, 5), ch(q0, 6), ch(q0, 7), out);
    }
  }

  // 2 args - face, meeting at a plane x=k
  template<int lev, typename iter_t>
  void dual_2x(node_t *q0, node_t *q1, iter_t &out) const {
    if (is_leaf(q0) && is_leaf(q1)) {
      if (lev == 2) emit(q0, q1, out);
      return;
    }

    dual_2x<lev>(ch(q0, 1), ch(q1, 0), out);
    dual_2x<lev>(ch(q0, 3), ch(q1, 2), out);
    dual_2x<lev>(ch(q0, 5), ch(q1, 4), out);
    dual_2x<lev>(ch(q0, 7), ch(q1, 6), out);

    if (lev <= 1) {
      dual_1y<lev>(ch(q0, 1), ch(q1, 0), ch(q0, 5), ch(q1, 4), out);
      dual_1y<lev>(ch(q0, 3), ch(q1, 2), ch(q0, 7), ch(q1, 6), out);

      dual_1z<lev>(ch(q0, 1), ch(q1, 0), ch(q0, 3), ch(q1, 2), out);
      dual_1z<lev>(ch(q0, 5), ch(q1, 4), ch(q0, 7), ch(q1, 6), out);
    }

    if (lev == 0) {
      dual_0<lev>(ch(q0, 1), ch(q1, 0), ch(q0, 3), ch(q1, 2), ch(q0, 5), ch(q1, 4), ch(q0, 7), ch(q1, 6), out);
    }
  }

  // 2 args - face, meeting at a plane y=k
  template<int lev, typename iter_t>
  void dual_2y(node_t *q0, node_t *q1, iter_t &out) const {
    if (is_leaf(q0) && is_leaf(q1)) {
      if (lev == 2) emit(q0, q1, out);
      return;
    }

    dual_2y<lev>(ch(q0, 2), ch(q1, 0), out);
    dual_2y<lev>(ch(q0, 3), ch(q1, 1), out);
    dual_2y<lev>(ch(q0, 6), ch(q1, 4), out);
    dual_2y<lev>(ch(q0, 7), ch(q1, 5), out);

    if (lev <= 1) {
      dual_1x<lev>(ch(q0, 2), ch(q1, 0), ch(q0, 6), ch(q1, 4), out);
      dual_1x<lev>(ch(q0, 3), ch(q1, 1), ch(q0, 7), ch(q1, 5), out);

      dual_1z<lev>(ch(q0, 2), ch(q0, 3), ch(q1, 0), ch(q1, 1), out);
      dual_1z<lev>(ch(q0, 6), ch(q0, 7), ch(q1, 4), ch(q1, 5), out);
    }

    if (lev == 0) {
      dual_0<lev>(ch(q0, 2), ch(q0, 3), ch(q1, 0), ch(q1, 1), ch(q0, 6), ch(q0, 7), ch(q1, 4), ch(q1, 5), out);
    }
  }

  // 2 args - face, meeting at a plane z=k
  template<int lev, typename iter_t>
  void dual_2z(node_t *q0, node_t *q1, iter_t &out) const {
    if (is_leaf(q0) && is_leaf(q1)) {
      if (lev == 2) emit(q0, q1, out);
      return;
    }

    dual_2z<lev>(ch(q0, 4), ch(q1, 0), out);
    dual_2z<lev>(ch(q0, 5), ch(q1, 1), out);
    dual_2z<lev>(ch(q0, 6), ch(q1, 2), out);
    dual_2z<lev>(ch(q0, 7), ch(q1, 3), out);

    if (lev <= 1) {
      dual_1x<lev>(ch(q0, 4), ch(q0, 6), ch(q1, 0), ch(q1, 2), out);
      dual_1x<lev>(ch(q0, 5), ch(q0, 7), ch(q1, 1), ch(q1, 3), out);

      dual_1y<lev>(ch(q0, 4), ch(q0, 5), ch(q1, 0), ch(q1, 1), out);
      dual_1y<lev>(ch(q0, 6), ch(q0, 7), ch(q1, 2), ch(q1, 3), out);
    }

    if (lev == 0) {
      dual_0<lev>(ch(q0, 4), ch(q0, 5), ch(q0, 6), ch(q0, 7), ch(q1, 0), ch(q1, 1), ch(q1, 2), ch(q1, 3), out);
    }
  }

  // 4 args - edge, meeting at a line y=i, z=j
  template<int lev, typename iter_t>
  void dual_1x(node_t *q0, node_t *q1, node_t *q2, node_t *q3, iter_t &out) const {
    if (is_leaf(q0) && is_leaf(q1) && is_leaf(q2) && is_leaf(q3)) {
      if (lev == 1) emit(q0, q1, q2, q3, out);
      return;
    }

    dual_1x<lev>(ch(q0, 6), ch(q1, 4), ch(q2, 2), ch(q3, 0), out);
    dual_1x<lev>(ch(q0, 7), ch(q1, 5), ch(q2, 3), ch(q3, 1), out);

    if (lev == 0) {
      dual_0<lev>(ch(q0, 6), ch(q0, 7), ch(q1, 4), ch(q1, 5), ch(q2, 2), ch(q2, 3), ch(q3, 0), ch(q3, 1), out);
    }
  }

  // 4 args - edge, meeting at a line x=i, z=j
  template<int lev, typename iter_t>
  void dual_1y(node_t *q0, node_t *q1, node_t *q2, node_t *q3, iter_t &out) const {
    if (is_leaf(q0) && is_leaf(q1) && is_leaf(q2) && is_leaf(q3)) {
      if (lev == 1) emit(q0, q1, q2, q3, out);
      return;
    }

    dual_1y<lev>(ch(q0, 5), ch(q1, 4), ch(q2, 1), ch(q3, 0), out);
    dual_1y<lev>(ch(q0, 7), ch(q1, 6), ch(q2, 3), ch(q3, 2), out);

    if (lev == 0) {
      dual_0<lev>(ch(q0, 5), ch(q1, 4), ch(q0, 7), ch(q1, 6), ch(q2, 1), ch(q3, 0), ch(q2, 3), ch(q3, 2), out);
    }
  }

  // 4 args - edge, meeting at a line x=i, y=j
  template<int lev, typename iter_t>
  void dual_1z(node_t *q0, node_t *q1, node_t *q2, node_t *q3, iter_t &out) const {
    if (is_leaf(q0) && is_leaf(q1) && is_leaf(q2) && is_leaf(q3)) {
      if (lev == 1) emit(q0, q1, q2, q3, out);
      return;
    }

    dual_1z<lev>(ch(q0, 3), ch(q1, 2), ch(q2, 1), ch(q3, 0), out);
    dual_1z<lev>(ch(q0, 7), ch(q1, 6), ch(q2, 5), ch(q3, 4), out);

    if (lev == 0) {
      dual_0<lev>(ch(q0, 3), ch(q1, 2), ch(q2, 1), ch(q3, 0), ch(q0, 7), ch(q1, 6), ch(q2, 5), ch(q3, 4), out);
    }
  }

  // 8 args - vertex
  template<int lev, typename iter_t>
  void dual_0(node_t *q0, node_t *q1, node_t *q2, node_t *q3,
              node_t *q4, node_t *q5, node_t *q6, node_t *q7,
              iter_t &out) const {
    while (!is_leaf(q0)) q0 = q0->children[7];
    while (!is_leaf(q1)) q1 = q1->children[6];
    while (!is_leaf(q2)) q2 = q2->children[5];
    while (!is_leaf(q3)) q3 = q3->children[4];
    while (!is_leaf(q4)) q4 = q4->children[3];
    while (!is_leaf(q5)) q5 = q5->children[2];
    while (!is_leaf(q6)) q6 = q6->children[1];
    while (!is_leaf(q7)) q7 = q7->children[0];

    emit(q0, q1, q2, q3, q4, q5, q6, q7, out);
  }

public:
  template<int lev, typename iter_t>
  void generate(node_t *root, iter_t out) const {
    if (lev <= 3) {
      dual_3<lev>(root, out);
    }

    if (lev <= 2) {
      dual_2x<lev>(root, NULL, out);
      dual_2x<lev>(NULL, root, out);
      dual_2y<lev>(root, NULL, out);
      dual_2y<lev>(NULL, root, out);
      dual_2z<lev>(root, NULL, out);
      dual_2z<lev>(NULL, root, out);
    }

    if (lev <= 1) {
      dual_1x<lev>(root, NULL, NULL, NULL, out);
      dual_1x<lev>(NULL, root, NULL, NULL, out);
      dual_1x<lev>(NULL, NULL, root, NULL, out);
      dual_1x<lev>(NULL, NULL, NULL, root, out);
      dual_1y<lev>(root, NULL, NULL, NULL, out);
      dual_1y<lev>(NULL, root, NULL, NULL, out);
      dual_1y<lev>(NULL, NULL, root, NULL, out);
      dual_1y<lev>(NULL, NULL, NULL, root, out);
      dual_1z<lev>(root, NULL, NULL, NULL, out);
      dual_1z<lev>(NULL, root, NULL, NULL, out);
      dual_1z<lev>(NULL, NULL, root, NULL, out);
      dual_1z<lev>(NULL, NULL, NULL, root, out);
    }

    if (lev <= 0) {
      dual_0<lev>(root, NULL, NULL, NULL, NULL, NULL, NULL, NULL, out);
      dual_0<lev>(NULL, root, NULL, NULL, NULL, NULL, NULL, NULL, out);
      dual_0<lev>(NULL, NULL, root, NULL, NULL, NULL, NULL, NULL, out);
      dual_0<lev>(NULL, NULL, NULL, root, NULL, NULL, NULL, NULL, out);
      dual_0<lev>(NULL, NULL, NULL, NULL, root, NULL, NULL, NULL, out);
      dual_0<lev>(NULL, NULL, NULL, NULL, NULL, root, NULL, NULL, out);
      dual_0<lev>(NULL, NULL, NULL, NULL, NULL, NULL, root, NULL, out);
      dual_0<lev>(NULL, NULL, NULL, NULL, NULL, NULL, NULL, root, out);
    }
  }
};
