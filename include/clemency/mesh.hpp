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

#include <tr1/unordered_set>
#include <tr1/unordered_map>

#include <set>

#include <clemency/geom.hpp>
#include <clemency/djset.hpp>



struct hash_size_t_pair {
  size_t operator()(const std::pair<size_t, size_t> &a) const {
    return (a.first * 131) ^ (a.second * 541);
  }
};



template<typename tri_t>
class tri_base_t {
private:
  tri_base_t &operator=(const tri_base_t &);
  tri_base_t(const tri_base_t &);

public:
  tri_t *next, *prev;
  size_t tag;

  union {
    struct { size_t a, b, c; };
    size_t v[3];
  };

  union {
    struct { tri_t *va, *vb, *vc; };
    tri_t *vptr[3];
  };

  tri_base_t(size_t _a, size_t _b, size_t _c) :
    tag(0),
    next(static_cast<tri_t *>(this)), prev(static_cast<tri_t *>(this)),
    a(_a), b(_b), c(_c),
    va(NULL), vb(NULL), vc(NULL) {
  }

  size_t vidx(size_t v) {
    if (a == v) return 0;
    if (b == v) return 1;
    if (c == v) return 2;
    return 3;
  }

  template<typename vlookup_t>
  v3d_t centroid(const vlookup_t &v) const {
    return (v(a) + v(b) + v(c)) / 3.0;
  }

  template<typename vlookup_t>
  aabb3d_t aabb(const vlookup_t &v) const {
    return aabb3d_t::initWithPoints(v(a), v(b), v(c));
  }

  template<typename vlookup_t>
  double area(const vlookup_t &v) const {
    return 0.5 * v3d_t::cross(v(b) - v(a), v(c) - v(a)).length();
    // double lab = (v(a) - v(b)).length();
    // double lbc = (v(b) - v(c)).length();
    // double lca = (v(c) - v(a)).length();
    // double s = (lab + lbc + lca) / 2.0;
    // return sqrt(s * (s - lab) * (s - lbc) * (s - lca));
  }

  template<typename vlookup_t>
  v3d_t normal(const vlookup_t &v) const {
    const v3d_t &pa = v(a), &pb = v(b), &pc = v(c);
    if (pa == pb || pa == pc || pb == pc) {
      std::cerr << "warning: calculating normal of degenerate triangle" << std::endl;
      return v3d_t::zero();
    }

    return v3d_t::cross(pb - pa, pc - pa).normalize();
  }
};



template<typename vert_t, typename tri_t>
class vert_base_t {
public:
  tri_t *face;
  v3d_t pos;

  vert_base_t() : face(NULL) {
  }

  vert_base_t(const v3d_t &_pos) : face(NULL), pos(_pos) {
  }

  void copy(const vert_base_t &v) {
    pos = v.pos;
  }
};



class tri_t : public tri_base_t<tri_t> {
  typedef tri_base_t<tri_t> super_t;

public:
  tri_t(size_t _a, size_t _b, size_t _c) : super_t(_a, _b, _c) {
  }
};



class vert_t : public vert_base_t<vert_t, tri_t> {
  typedef vert_base_t<vert_t, tri_t> super_t;

public:
  vert_t() : super_t() {
  }

  vert_t(const v3d_t &pos) : super_t(pos) {
  }

  void copy(const vert_t &v) {
    super_t::copy(v);
  }
};



struct vert_hash_t {
  size_t operator()(const vert_t &a) const {
    std::tr1::hash<double> h;
    size_t r = 0;
    r *= 131; r ^= h(a.pos.x);
    r *= 131; r ^= h(a.pos.y);
    r *= 131; r ^= h(a.pos.z);
    return r;
  }
};



struct vert_eq_t {
  bool operator()(const vert_t &a, const vert_t &b) const {
    return a.pos == b.pos;
  }
};



template<typename _vert_t, typename _tri_t>
struct triangle_mesh_t {
  typedef _vert_t vert_t;
  typedef _tri_t tri_t;

  struct face_iter {
    typedef std::bidirectional_iterator_tag iterator_category;
    typedef tri_t                           value_type;
    typedef ptrdiff_t                       difference_type;
    typedef tri_t *                         pointer;
    typedef tri_t &                         reference;
    tri_t *node;
    face_iter() : node() {}
    explicit face_iter(tri_t *_node) : node(_node) {}
    reference operator*() const { return *node; }
    pointer operator->() const { return node; }
    face_iter &operator++() { node = node->next; return *this; }
    face_iter operator++(int) { face_iter tmp = *this; node = node->next; return tmp; }
    face_iter &operator--() { node = node->prev; return *this; }
    face_iter operator--(int) { face_iter tmp = *this; node = node->prev; return tmp; }
    bool operator==(const face_iter &other) const { return node == other.node; }
    bool operator!=(const face_iter &other) const { return node != other.node; }
  };

  struct const_face_iter {
    typedef std::bidirectional_iterator_tag iterator_category;
    typedef tri_t                           value_type;
    typedef ptrdiff_t                       difference_type;
    typedef const tri_t *                   pointer;
    typedef const tri_t &                   reference;
    const tri_t *node;
    const_face_iter() : node() {}
    explicit const_face_iter(const tri_t *_node) : node(_node) {}
    reference operator*() const { return *node; }
    pointer operator->() const { return node; }
    const_face_iter &operator++() { node = node->next; return *this; }
    const_face_iter operator++(int) { const_face_iter tmp = *this; node = node->next; return tmp; }
    const_face_iter &operator--() { node = node->prev; return *this; }
    const_face_iter operator--(int) { const_face_iter tmp = *this; node = node->prev; return tmp; }
    bool operator==(const const_face_iter &other) const { return node == other.node; }
    bool operator!=(const const_face_iter &other) const { return node != other.node; }
  };

  size_t curr_tag;
  tri_t face_list;
  size_t face_count;
  std::vector<vert_t> vertices;

  triangle_mesh_t() : curr_tag(0), face_list(0,0,0), face_count(0), vertices() {
  }

  ~triangle_mesh_t() {
    clear_tris();
  }

  face_iter fbegin() { return face_iter(face_list.next); }
  face_iter fend() { return face_iter(&face_list); }

  const_face_iter fbegin() const { return const_face_iter(face_list.next); }
  const_face_iter fend() const { return const_face_iter(&face_list); }

  struct vlookup_t {
    const std::vector<vert_t> &vertices;
    vlookup_t(const std::vector<vert_t> &_vertices) : vertices(_vertices) {
    }
    const v3d_t &operator()(const size_t &v) const {
      return vertices[v].pos;
    }
  };

  vlookup_t vertex_getter() const { return vlookup_t(vertices); }

  aabb3d_t aabb(const tri_t &t) const {
    return t.aabb(vlookup_t(vertices));
  }

  v3d_t normal(const tri_t &t) const {
    return t.normal(vlookup_t(vertices));
  }

  v3d_t centroid(const tri_t &t) {
    return t.centroid(vlookup_t(vertices));
  }

  struct aabb_calc_t {
    const triangle_mesh_t *mesh;
    aabb_calc_t(const triangle_mesh_t *_mesh) : mesh(_mesh) {
    }
    aabb3d_t operator()(const tri_t &t) const { return mesh->aabb(t); }
    aabb3d_t operator()(const tri_t *t) const { return mesh->aabb(*t); }
  };

  aabb_calc_t aabb_getter() const { return aabb_calc_t(this); }

  void vnext(tri_t **&x, size_t a) {
    tri_t *t = *x;
    size_t j = t->vidx(a);
    assert(j != 3);
    x = &(t->vptr[j]);
  }

  void vnext(tri_t *&t, size_t a) {
    size_t j = t->vidx(a);
    assert(j != 3);
    t = t->vptr[j];
  }

  template<typename iter_t>
  void find_triangles(size_t a, size_t b, iter_t iter) {
    tri_t *ta = vertices[a].face;
    tri_t *tb = vertices[b].face;
    while (ta && tb) {
      if (ta == tb) {
        *iter = ta;
        ++iter;
        vnext(ta, a);
        vnext(tb ,b);
      } else if (ta < tb) {
        vnext(ta, a);
      } else {
        vnext(tb, b);
      }
    }
  }

  tri_t *find_tri(size_t a, size_t b) {
    tri_t *t = vertices[a].face;
    while (t != NULL) {
      size_t i = t->vidx(a);
      if (t->v[(i+1)%3] == b) return t;
      t = t->vptr[i];
    }
    return NULL;
  }

  void insert_tri_vert(tri_t *t, size_t i) {
    size_t a = t->v[i];

    tri_t **x = &(vertices[a].face);

    while (*x != NULL && *x < t) vnext(x, a);

    t->vptr[i] = *x;
    *x = t;
  }

  void remove_tri_vert(tri_t *t, size_t i) {
    size_t a = t->v[i];

    tri_t **x = &(vertices[a].face);

    while (*x != NULL && *x < t) {
      vnext(x, a);
    }
    assert (*x == t);

    *x = t->vptr[i];
    t->vptr[i] = NULL;
  }

  size_t add_vertex(const v3d_t &pos) {
    vertices.push_back(vert_t(pos));
    return vertices.size() - 1;
  }

  tri_t *add_tri(size_t a, size_t b, size_t c) {
    if (a == b || a == c || b == c) {
      return NULL;
    }

    tri_t *t = new tri_t(a, b, c);

    t->prev = face_list.prev;
    t->next = &face_list;

    face_list.prev->next = t;
    face_list.prev = t;

    for (size_t i = 0; i < 3; ++i) insert_tri_vert(t, i);

    ++face_count;

    return t;
  }

  tri_t *remove_tri(tri_t *t) {
    t->next->prev = t->prev;
    t->prev->next = t->next;

    t->next = t->prev = t;

    // std::cerr << "remove tri " << t << " ; " << t->a << "," << t->b << "," << t->c << std::endl;
    // std::cerr << "  a "; for (tri_t *t2 = vertices[t->a].face; t2 != NULL; vnext(t2, t->a)) std::cerr << " " << t2; std::cerr << std::endl;
    // std::cerr << "  b "; for (tri_t *t2 = vertices[t->b].face; t2 != NULL; vnext(t2, t->b)) std::cerr << " " << t2; std::cerr << std::endl;
    // std::cerr << "  c "; for (tri_t *t2 = vertices[t->c].face; t2 != NULL; vnext(t2, t->c)) std::cerr << " " << t2; std::cerr << std::endl;

    for (size_t i = 0; i < 3; ++i) remove_tri_vert(t, i);

    --face_count;

    return t;
  }

  void clear_tris() {
    tri_t *t_next;

    for (tri_t *t = face_list.next; t != &face_list; t = t_next) {
      t_next = t->next;
      delete t;
    }
    face_list.prev = face_list.next = &face_list;

    face_count = 0;
    for (size_t i = 0; i < vertices.size(); ++i) {
      vertices[i].face = NULL;
    }
  }

  void rewrite_vertex(tri_t *t, size_t a, size_t b) {
    size_t i = t->vidx(a);
    assert(i != 3);
    remove_tri_vert(t, i);
    t->v[i] = b;
    insert_tri_vert(t, i);
  }

  void trace(size_t a) {
    tri_t *t = vertices[a].face;
    while (t) { vnext(t, a); }
  }

  void collapse_edge(size_t a, size_t b) {
    // std::cerr << "collapse_edge(" << a << "," << b << ")" << std::endl;

    // {
    //   std::cerr << "in triangles incident to " << a << " = " << vertex_face_count(a) << std::endl;
    //   tri_t *t = vertices[a].face;
    //   while (t) {
    //     std::cerr << "    -> " << t << " " << t->a << ", " << t->b << ", " << t->c << std::endl;
    //     vnext(t, a);
    //   }
    //   std::cerr << std::endl;
    // }
    // {
    //   std::cerr << "in triangles incident to " << b << " = " << vertex_face_count(b) << std::endl;
    //   tri_t *t = vertices[b].face;
    //   while (t) {
    //     std::cerr << "    -> " << t << " " << t->a << ", " << t->b << ", " << t->c << std::endl;
    //     vnext(t, b);
    //   }
    //   std::cerr << std::endl;
    // }

    // remove triangles that have an edge ab or ba.
    // rewrite triangles with a vertex b to use a.
    if (!vertices[a].face || !vertices[b].face) return;

    tri_t **t = &(vertices[a].face);
    tri_t *u = vertices[b].face;

    size_t n_deleted = 0;
    while (u != NULL) {
      if (*t == u) {
        tri_t *tmp = u;
        size_t aidx = tmp->vidx(a);
        size_t bidx = tmp->vidx(b);
        size_t cidx = (aidx+1)%3; if (cidx == bidx) cidx = (bidx+1)%3;
        size_t c = tmp->v[cidx];
        assert(c != a && c != b);
        *t = tmp->vptr[aidx];        // unlink from a list
        u = tmp->vptr[bidx];         // advance u
        remove_tri_vert(tmp, cidx);  // unlink from c list
        --face_count;                // decrease face count
        tmp->next->prev = tmp->prev; // unlink from face list
        tmp->prev->next = tmp->next;
        n_deleted++;
        delete tmp;                  // free
      } else if (*t != NULL && *t < u) {
        vnext(t, a);                 // advance t
      } else {
        // insert the head of u into t.
        tri_t *tmp = u;
        size_t bidx = tmp->vidx(b);
        u = tmp->vptr[bidx];
        tmp->v[bidx] = a;            // remap b to a
        tmp->vptr[bidx] = *t;        // link tmp into t
        *t = tmp;
        t = &(tmp->vptr[bidx]);      // advance t
      }

      vertices[b].face = NULL;
    }
  }

  size_t vertex_face_count(size_t a) {
    size_t n = 0;
    tri_t *t = vertices[a].face;
    while (t) {
      vnext(t, a);
      ++n;
    }
    return n;
  }

  template<typename iter_t>
  void tag_faces(iter_t beg, iter_t end) {
    ++curr_tag;
    while (beg != end) {
      (*beg)->tag = curr_tag;
      ++beg;
    }
  }

  void tag_faces(tri_t *beg, tri_t *end) {
    ++curr_tag;
    while (beg != end) {
      beg->tag = curr_tag;
      beg = beg->next;
    }
  }

  inline bool is_tagged(tri_t *t, size_t tag) {
    return t != NULL && t->tag == tag;
  }

  bool all_faces_tagged(size_t vert, size_t tag) {
    tri_t *t = vertices[vert].face;
    while (t) {
      if (!is_tagged(t, tag)) return false;
      vnext(t, vert);
    }
    return true;
  }



  template<typename iter_t>
  void get_face_vertices(tri_t *t, iter_t iter) {
    *iter++ = t->a;
    *iter++ = t->b;
    *iter++ = t->c;
  }

  template<typename tri_iter_t, typename iter_t>
  void get_face_vertices(tri_iter_t beg, tri_iter_t end, iter_t iter) {
    for (;beg != end; ++beg) {
      get_face_vertices(*beg, iter);
    }
  }

  template<typename tri_iter_t, typename iter_t>
  void get_face_vertex_set(tri_iter_t beg, tri_iter_t end, iter_t iter) {
    std::set<size_t> tmp;
    get_face_vertices(beg, end, std::inserter(tmp, tmp.end()));
    std::copy(tmp.begin(), tmp.end(), iter);
  }



  template<typename iter_t>
  void get_face_edges(tri_t *t, iter_t iter) {
    *iter++ = std::make_pair(t->a, t->b);
    *iter++ = std::make_pair(t->b, t->c);
    *iter++ = std::make_pair(t->c, t->a);
  }

  template<typename tri_iter_t, typename iter_t>
  void get_face_edges(tri_iter_t beg, tri_iter_t end, iter_t iter) {
    for (;beg != end; ++beg) {
      get_face_edges(*beg, iter);
    }
  }

  template<typename tri_iter_t, typename iter_t>
  void get_face_edge_set(tri_iter_t beg, tri_iter_t end, iter_t iter) {
    std::set<std::pair<size_t, size_t> > tmp;
    get_face_edges(beg, end, std::inserter(tmp, tmp.end()));
    std::copy(tmp.begin(), tmp.end(), iter);
  }



  template<typename iter_t>
  void get_vertex_triangles(size_t a, iter_t iter) {
    for (tri_t *t = vertices[a].face; t; vnext(t, a)) *iter++ = t;
  }

  template<typename vert_iter_t, typename iter_t>
  void get_vertex_triangles(vert_iter_t beg, vert_iter_t end, iter_t iter) {
    for (;beg != end; ++beg) {
      get_vertex_triangles(*beg, iter);
    }
  }

  template<typename vert_iter_t, typename iter_t>
  void get_vertex_triangle_set(vert_iter_t beg, vert_iter_t end, iter_t iter) {
    std::set<tri_t *> tmp;
    for (;beg != end; ++beg) {
      get_vertex_triangles(*beg, std::inserter(tmp, tmp.end()));
    }
    std::copy(tmp.begin(), tmp.end(), iter);
  }



  template<typename iter_t>
  double _volume(iter_t begin, iter_t end, tri_t) const {
    double vol = 0.0;
    for (; begin != end; ++begin) {
      tri_t &t = *begin;
      vol += v3d_t::tetrahedron_volume(vertices[t.a].pos, vertices[t.b].pos, vertices[t.c].pos, v3d_t::zero());
    }
    return vol;
  }

  template<typename iter_t>
  double _volume(iter_t begin, iter_t end, tri_t *) const {
    double vol = 0.0;
    for (; begin != end; ++begin) {
      tri_t *t = *begin;
      vol += v3d_t::tetrahedron_volume(vertices[t->a].pos, vertices[t->b].pos, vertices[t->c].pos, v3d_t::zero());
    }
    return vol;
  }

  template<typename iter_t>
  double volume(iter_t begin, iter_t end) const {
    return _volume(begin, end, typename std::iterator_traits<iter_t>::value_type());
  }

  double volume() const {
    return volume(fbegin(), fend());
  }

  template<typename iter_t>
  v3d_t _centre_of_mass(iter_t begin, iter_t end, tri_t) const {
    v3d_t com = v3d_t::zero();
    double vol;
    for (; begin != end; ++begin) {
      tri_t &t = *begin;
      double tvol = v3d_t::tetrahedron_volume(vertices[t.a].pos, vertices[t.b].pos, vertices[t.c].pos, v3d_t::zero());
      com += tvol * (vertices[t.a].pos + vertices[t.b].pos + vertices[t.c].pos);
      vol += tvol;
    }
    return com / (vol * 4.0);
  }

  template<typename iter_t>
  v3d_t _centre_of_mass(iter_t begin, iter_t end, tri_t *) const {
    v3d_t com = v3d_t::zero();
    double vol;
    for (; begin != end; ++begin) {
      tri_t *t = *begin;
      double tvol = v3d_t::tetrahedron_volume(vertices[t->a].pos, vertices[t->b].pos, vertices[t->c].pos, v3d_t::zero());
      com += tvol * (vertices[t->a].pos + vertices[t->b].pos + vertices[t->c].pos);
      vol += tvol;
    }
    return com / (vol * 4.0);
  }

  template<typename iter_t>
  v3d_t centre_of_mass(iter_t begin, iter_t end) const {
    return _centre_of_mass(begin, end, typename std::iterator_traits<iter_t>::value_type());
  }

  v3d_t centre_of_mass() const {
    return centre_of_mass(fbegin(), fend());
  }

  bool is_closed_manifold() const;

  template<typename iter_t>
  static tri_t *_tri_ptr(iter_t iter, tri_t) { return &*iter; }

  template<typename iter_t>
  static tri_t *_tri_ptr(iter_t iter, tri_t *) { return *iter; }

  template<typename iter_t>
  static tri_t *tri_ptr(iter_t iter) {
    return _tri_ptr(iter, typename std::iterator_traits<iter_t>::value_type());
  }

  template<typename iter_t>
  triangle_mesh_t *submesh(iter_t begin, iter_t end) const {
    triangle_mesh_t *result = new triangle_mesh_t;

    std::tr1::unordered_map<size_t, size_t> vertex_map;
    std::vector<size_t> vertex_idx;

    for (iter_t i = begin; i != end; ++i) {
      tri_t *t = tri_ptr(i);
      for (size_t j = 0; j < 3; ++j) {
        std::tr1::unordered_map<size_t, size_t>::iterator iter = vertex_map.find(t->v[j]);
        if (iter == vertex_map.end()) {
          vertex_map[t->v[j]] = vertex_idx.size();
          vertex_idx.push_back(t->v[j]);
        }
      }
    }

    result->vertices.reserve(vertex_idx.size());

    for (size_t i = 0; i < vertex_idx.size(); ++i) {
      result->add_vertex(vertices[vertex_idx[i]].pos);
    }

    for (iter_t i = begin; i != end; ++i) {
      tri_t *t = tri_ptr(i);
      result->add_tri(vertex_map[t->a], vertex_map[t->b], vertex_map[t->c]);
    }

    return result;
  }
};



template<typename mesh_t, typename vert_hash_t, typename vert_eq_t>
mesh_t *merge_meshes(mesh_t *a, mesh_t *b) {
  typedef std::tr1::unordered_map<typename mesh_t::vert_t, size_t, vert_hash_t, vert_eq_t> vert_map_t;
  vert_map_t vert_map;

  std::vector<size_t> a_remap, b_remap;

  a_remap.reserve(a->vertices.size());
  b_remap.reserve(b->vertices.size());

  for (size_t i = 0; i < a->vertices.size(); ++i) {
    typename vert_map_t::iterator vi = vert_map.find(a->vertices[i]);
    if (vi != vert_map.end()) {
      a_remap.push_back((*vi).second);
    } else {
      a_remap.push_back(vert_map.size());
      vert_map[a->vertices[i]] = a_remap.back();
    }
  }

  for (size_t i = 0; i < b->vertices.size(); ++i) {
    typename vert_map_t::iterator vi = vert_map.find(b->vertices[i]);
    if (vi != vert_map.end()) {
      b_remap.push_back((*vi).second);
    } else {
      b_remap.push_back(vert_map.size());
      vert_map[b->vertices[i]] = b_remap.back();
    }
  }

  mesh_t *result = new mesh_t;

  result->vertices.resize(vert_map.size());

  for (size_t i = 0; i < a_remap.size(); ++i) {
    result->vertices[a_remap[i]].copy(a->vertices[i]);
  }

  for (size_t i = 0; i < b_remap.size(); ++i) {
    result->vertices[b_remap[i]].copy(b->vertices[i]);
  }

  for (typename mesh_t::tri_t *t = a->face_list.next; t != &a->face_list; t = t->next) {
    result->add_tri(a_remap[t->a], a_remap[t->b], a_remap[t->c]);
  }

  for (typename mesh_t::tri_t *t = b->face_list.next; t != &b->face_list; t = t->next) {
    result->add_tri(b_remap[t->a], b_remap[t->b], b_remap[t->c]);
  }

  return result;
}

#include <clemency/volume_integrals.hpp>

template<typename _vert_t, typename _tri_t>
bool triangle_mesh_t<_vert_t, _tri_t>::is_closed_manifold() const {
    typedef std::tr1::unordered_map<std::pair<size_t, size_t>, std::vector<size_t>, hash_size_t_pair> edge_map_t;

    std::vector<tri_t *> faces;
    faces.reserve(face_count);
    for (tri_t *t = face_list.next; t != &face_list; t = t->next) {
      faces.push_back(t);
    }

    edge_map_t edge_map;

    for (size_t i = 0; i < faces.size(); ++i) {
      tri_t *t = faces[i];
      edge_map[std::make_pair(t->a, t->b)].push_back(i);
      edge_map[std::make_pair(t->b, t->c)].push_back(i);
      edge_map[std::make_pair(t->c, t->a)].push_back(i);
    }

    djset connected_components(faces.size());

    std::set<tri_t *> border_faces;
    std::vector<std::set<tri_t *> > face_groups;

    for (typename edge_map_t::iterator i = edge_map.begin(); i != edge_map.end(); ++i) {
      typename edge_map_t::iterator j = edge_map.find(std::make_pair((*i).first.second, (*i).first.first));
      if (j == edge_map.end() || (*i).second.size() != 1 || (*j).second.size() != 1) {
        for (size_t k = 0; k < (*i).second.size(); ++k) {
          border_faces.insert(faces[(*i).second[k]]);
        }
        continue;
      }
      connected_components.merge_sets((*i).second[0], (*j).second[0]);
    }
    std::cerr << "connected component count: " << connected_components.count() << std::endl;
    connected_components.collate(faces.begin(), face_groups);
    for (size_t i = 0; i < face_groups.size(); ++i) {
      std::set<tri_t *> bf;
      std::set_intersection(face_groups[i].begin(), face_groups[i].end(),
                            border_faces.begin(), border_faces.end(),
                            std::inserter(bf, bf.end()));
      std::cerr << "  component " << i << "(" << face_groups[i].size() << " faces) ";
      if (bf.size()) {
        std::cerr << "not closed" << std::endl;
      } else {
        volint_t<triangle_mesh_t<_vert_t, _tri_t> > vol(this, face_groups[i].begin(), face_groups[i].end());
        std::cerr << "closed"
                  << " volume=" << volume(face_groups[i].begin(), face_groups[i].end())
                  << " com=" << centre_of_mass(face_groups[i].begin(), face_groups[i].end())
                  << " volume(2)=" << vol.volume()
                  << " com(2)=" << vol.centre_of_mass()
                  << std::endl;
      }
    }
    return true;
  }
