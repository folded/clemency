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



#include "v3.hpp"
#include "aabb3.hpp"



template<typename child_t, typename num_t>
class octree_t {
public:
  aabb3_t<num_t> bbox;
  child_t *children[8];

  ~octree_t() {
    clear_children();
  };

  octree_t(const aabb3_t<num_t> &_bbox) : bbox(_bbox) {
    std::fill(children, children + 8, (child_t *)NULL);
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

  v3_t<num_t> corner(unsigned c) const {
    return bbox.corner(c);
  }

  void split(octree_t *node, const v3_t<num_t> &expand = v3_t<num_t>::zero()) {
    if (children[0]) throw std::runtime_error("node already split");

    v3_t<num_t> lo = bbox.mid - bbox.extent;
    v3_t<num_t> hi = bbox.mid + bbox.extent;
    v3_t<num_t> md = bbox.mid;

    children[0] = new child_t(aabb3_t<num_t>::initWithPoints(v3_t<num_t>::init(lo.x, lo.y, lo.z) - expand, v3_t<num_t>::init(md.x, md.y, md.z) + expand));
    children[1] = new child_t(aabb3_t<num_t>::initWithPoints(v3_t<num_t>::init(md.x, lo.y, lo.z) - expand, v3_t<num_t>::init(hi.x, md.y, md.z) + expand));
    children[2] = new child_t(aabb3_t<num_t>::initWithPoints(v3_t<num_t>::init(lo.x, md.y, lo.z) - expand, v3_t<num_t>::init(md.x, hi.y, md.z) + expand));
    children[3] = new child_t(aabb3_t<num_t>::initWithPoints(v3_t<num_t>::init(md.x, md.y, lo.z) - expand, v3_t<num_t>::init(hi.x, hi.y, md.z) + expand));

    children[4] = new child_t(aabb3_t<num_t>::initWithPoints(v3_t<num_t>::init(lo.x, lo.y, md.z) - expand, v3_t<num_t>::init(md.x, md.y, hi.z) + expand));
    children[5] = new child_t(aabb3_t<num_t>::initWithPoints(v3_t<num_t>::init(md.x, lo.y, md.z) - expand, v3_t<num_t>::init(hi.x, md.y, hi.z) + expand));
    children[6] = new child_t(aabb3_t<num_t>::initWithPoints(v3_t<num_t>::init(lo.x, md.y, md.z) - expand, v3_t<num_t>::init(md.x, hi.y, hi.z) + expand));
    children[7] = new child_t(aabb3_t<num_t>::initWithPoints(v3_t<num_t>::init(md.x, md.y, md.z) - expand, v3_t<num_t>::init(hi.x, hi.y, hi.z) + expand));
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
