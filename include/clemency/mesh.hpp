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
  v3_t centroid(const vlookup_t &v) {
    return (v(a) + v(b) + v(c)) / 3.0;
  }

  template<typename vlookup_t>
  double area(const vlookup_t &v) {
    return 0.5 * v3_t::cross(v(b) - v(a), v(c) - v(a)).length();
    // double lab = (v(a) - v(b)).length();
    // double lbc = (v(b) - v(c)).length();
    // double lca = (v(c) - v(a)).length();
    // double s = (lab + lbc + lca) / 2.0;
    // return sqrt(s * (s - lab) * (s - lbc) * (s - lca));
  }

  template<typename vlookup_t>
  v3_t normal(const vlookup_t &v) {
    v3_t pa, pb, pc;
    pa = v(a); pb = v(b); pc = v(c);
    if (pa == pb || pa == pc || pb == pc) {
      std::cerr << "warning: calculating normal of degenerate triangle" << std::endl;
      return v3_t::zero();
    }

    return v3_t::cross(pb - pa, pc - pa).normalize();
  }
};



template<typename vert_t, typename tri_t>
class vert_base_t {
public:
  tri_t *face;
  v3_t pos;

  vert_base_t() : face(NULL) {
  }

  vert_base_t(const v3_t &_pos) : face(NULL), pos(_pos) {
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

  vert_t(const v3_t &pos) : super_t(pos) {
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

  size_t curr_tag;
  _tri_t face_list;
  size_t face_count;
  std::vector<_vert_t> vertices;

  triangle_mesh_t() : curr_tag(0), face_list(0,0,0), face_count(0), vertices() {
  }

  ~triangle_mesh_t() {
    clear_tris();
  }

  v3_t centroid(const _tri_t &t) {
    return (vertices[t.a].pos +
            vertices[t.b].pos +
            vertices[t.c].pos) / 3.0;
  }

  v3_t normal(const _tri_t &t) {
    if (vertices[t.a].pos == vertices[t.b].pos ||
        vertices[t.a].pos == vertices[t.c].pos ||
        vertices[t.b].pos == vertices[t.c].pos) {
      std::cerr << "warning: calculating normal of degenerate triangle" << std::endl;
      return v3_t::zero();
    }

    return v3_t::cross(vertices[t.b].pos - vertices[t.a].pos,
                       vertices[t.c].pos - vertices[t.a].pos).normalize();
  }

  void vnext(_tri_t **&x, size_t a) {
    _tri_t *t = *x;
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
    _tri_t *ta = vertices[a].face;
    _tri_t *tb = vertices[b].face;
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

  _tri_t *find_tri(size_t a, size_t b) {
    _tri_t *t = vertices[a].face;
    while (t != NULL) {
      size_t i = t->vidx(a);
      if (t->v[(i+1)%3] == b) return t;
      t = t->vptr[i];
    }
    return NULL;
  }

  void insert_tri_vert(_tri_t *t, size_t i) {
    size_t a = t->v[i];

    _tri_t **x = &(vertices[a].face);

    while (*x != NULL && *x < t) vnext(x, a);

    t->vptr[i] = *x;
    *x = t;
  }

  void remove_tri_vert(_tri_t *t, size_t i) {
    size_t a = t->v[i];

    _tri_t **x = &(vertices[a].face);

    while (*x != NULL && *x < t) {
      vnext(x, a);
    }
    assert (*x == t);

    *x = t->vptr[i];
    t->vptr[i] = NULL;
  }

  size_t add_vertex(const v3_t &pos) {
    vertices.push_back(vert_t(pos));
    return vertices.size() - 1;
  }

  _tri_t *add_tri(size_t a, size_t b, size_t c) {
    if (a == b || a == c || b == c) {
      return NULL;
    }

    _tri_t *t = new _tri_t(a, b, c);

    t->prev = face_list.prev;
    t->next = &face_list;

    face_list.prev->next = t;
    face_list.prev = t;

    for (size_t i = 0; i < 3; ++i) insert_tri_vert(t, i);

    ++face_count;

    return t;
  }

  _tri_t *remove_tri(_tri_t *t) {
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
    _tri_t *t_next;

    for (_tri_t *t = face_list.next; t != &face_list; t = t_next) {
      t_next = t->next;
      delete t;
    }
    face_list.prev = face_list.next = &face_list;

    face_count = 0;
    for (size_t i = 0; i < vertices.size(); ++i) {
      vertices[i].face = NULL;
    }
  }

  void rewrite_vertex(_tri_t *t, size_t a, size_t b) {
    size_t i = t->vidx(a);
    assert(i != 3);
    remove_tri_vert(t, i);
    t->v[i] = b;
    insert_tri_vert(t, i);
  }

  void trace(size_t a) {
    _tri_t *t = vertices[a].face;
    while (t) { vnext(t, a); }
  }

  void collapse_edge(size_t a, size_t b) {
    // std::cerr << "collapse_edge(" << a << "," << b << ")" << std::endl;

    // {
    //   std::cerr << "in triangles incident to " << a << " = " << vertex_face_count(a) << std::endl;
    //   _tri_t *t = vertices[a].face;
    //   while (t) {
    //     std::cerr << "    -> " << t << " " << t->a << ", " << t->b << ", " << t->c << std::endl;
    //     vnext(t, a);
    //   }
    //   std::cerr << std::endl;
    // }
    // {
    //   std::cerr << "in triangles incident to " << b << " = " << vertex_face_count(b) << std::endl;
    //   _tri_t *t = vertices[b].face;
    //   while (t) {
    //     std::cerr << "    -> " << t << " " << t->a << ", " << t->b << ", " << t->c << std::endl;
    //     vnext(t, b);
    //   }
    //   std::cerr << std::endl;
    // }

    // remove triangles that have an edge ab or ba.
    // rewrite triangles with a vertex b to use a.
    if (!vertices[a].face || !vertices[b].face) return;

    _tri_t **t = &(vertices[a].face);
    _tri_t *u = vertices[b].face;

    size_t n_deleted = 0;
    while (u != NULL) {
      if (*t == u) {
        _tri_t *tmp = u;
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
        _tri_t *tmp = u;
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
    _tri_t *t = vertices[a].face;
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

  void tag_faces(_tri_t *beg, _tri_t *end) {
    ++curr_tag;
    while (beg != end) {
      beg->tag = curr_tag;
      beg = beg->next;
    }
  }

  inline bool is_tagged(_tri_t *t, size_t tag) {
    return t != NULL && t->tag == tag;
  }

  bool all_faces_tagged(size_t vert, size_t tag) {
    _tri_t *t = vertices[vert].face;
    while (t) {
      if (!is_tagged(t, tag)) return false;
      vnext(t, vert);
    }
    return true;
  }



  template<typename iter_t>
  void get_face_vertices(_tri_t *t, iter_t iter) {
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
  void get_face_edges(_tri_t *t, iter_t iter) {
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
    for (_tri_t *t = vertices[a].face; t; vnext(t, a)) *iter++ = t;
  }

  template<typename vert_iter_t, typename iter_t>
  void get_vertex_triangles(vert_iter_t beg, vert_iter_t end, iter_t iter) {
    for (;beg != end; ++beg) {
      get_vertex_triangles(*beg, iter);
    }
  }

  template<typename vert_iter_t, typename iter_t>
  void get_vertex_triangle_set(vert_iter_t beg, vert_iter_t end, iter_t iter) {
    std::set<_tri_t *> tmp;
    for (;beg != end; ++beg) {
      get_vertex_triangles(*beg, std::inserter(tmp, tmp.end()));
    }
    std::copy(tmp.begin(), tmp.end(), iter);
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
