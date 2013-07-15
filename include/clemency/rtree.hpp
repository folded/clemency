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

#include <clemency/geom.hpp>
#include <clemency/mesh.hpp>

#include <iostream>

#include <cmath>
#include <limits>

template<typename _data_t>
struct rtree_node_t {
  typedef _data_t data_t;

  aabb3d_t bbox;
  rtree_node_t *child;
  rtree_node_t *sibling;
  std::vector<data_t *> data;

  struct data_aabb_t {
    data_t *data;
    aabb3d_t bbox;

    data_aabb_t() { }
    data_aabb_t(data_t *_data, aabb3d_t _bbox) : data(_data), bbox(_bbox) { }
    const aabb3d_t &aabb() const { return bbox; }
    operator data_t *() const { return data; }
  };

  const aabb3d_t &aabb() const {
    return bbox;
  }

  struct get_aabb3d_t {
    const aabb3d_t &operator()(const rtree_node_t *n) const { return n->aabb(); }
    const aabb3d_t &operator()(const data_aabb_t &d) const { return d.aabb(); }
  };

  template<typename iter_t, typename aabb_calc_t>
  aabb3d_t fit_aabb(iter_t begin, iter_t end, const aabb_calc_t &get_aabb) {
    aabb3d_t result = get_aabb(*begin);
    while (++begin != end) {
      result.expand(get_aabb(*begin));
    }
    return result;
  }

  // Fill an rtree node from a set of (data, aabb) pairs.
  template<typename iter_t>
  void _fill(iter_t begin, iter_t end, data_aabb_t) {
    data.clear();
    std::copy(begin, end, std::back_inserter(data));
    bbox = fit_aabb(begin, end, get_aabb3d_t());
  }

  // Fill an rtree node from a set of child nodes.
  template<typename iter_t>
  void _fill(iter_t begin, iter_t end, rtree_node_t *) {
    iter_t i = begin;
    rtree_node_t *curr = child = *i;
    while (++i != end) {
      curr->sibling = *i;
      curr = curr->sibling;
    }
    bbox = fit_aabb(begin, end, get_aabb3d_t());
  }

  template<typename iter_t>
  rtree_node_t(iter_t begin, iter_t end) : bbox(), child(NULL), sibling(NULL), data() {
    _fill(begin, end, typename std::iterator_traits<iter_t>::value_type());
  }

  ~rtree_node_t() {
    if (child) {
      rtree_node_t *next = child;
      while (next) {
        rtree_node_t *curr = next;
        next = next->sibling;
        delete curr;
      }
    }
  }

  // Search the rtree for objects that intersect obj (generally an aabb).
  // The aabb class must provide a method intersects(obj_t).
  template<typename obj_t, typename out_iter_t>
  void search(const obj_t &obj, out_iter_t out) const {
    if (!bbox.intersects(obj)) return;
    if (child) {
      for (rtree_node_t *node = child; node; node = node->sibling) {
        node->search(obj, out);
      }
    } else {
      std::copy(data.begin(), data.end(), out);
    }
  }

  // update the bounding box extents of nodes that intersect obj (generally an aabb).
  // The aabb class must provide a method intersects(obj_t).
  template<typename obj_t, typename aabb_calc_t>
  void update_extents(const obj_t &obj, const aabb_calc_t &get_aabb) {
    if (!bbox.intersects(obj)) return;

    if (child) {
      rtree_node_t *node = child;
      node->update_extents(obj);
      bbox = node->bbox;
      for (node = node->sibling; node; node = node->sibling) {
        node->update_extents(obj);
        bbox.expand(node->bbox);
      }
    } else {
      bbox = fit_aabb(data.begin(), data.end(), get_aabb);
    }
  }

  template<typename aabb_calc_t>
  bool _remove(data_t *val, const aabb3d_t &val_aabb, const aabb_calc_t &get_aabb) {
    if (!bbox.intersects(val_aabb)) return false;

    if (child) {
      rtree_node_t *node = child;

      bool removed = node->_remove(val, val_aabb, get_aabb);
      bbox = node->bbox;

      for (node = node->sibling; node; node = node->sibling) {
        if (!removed) {
          removed = node->_remove(val, val_aabb, get_aabb);
        }
        bbox.expand(node->bbox);
      }

      return removed;
    } else {
      typename std::vector<data_t *>::iterator i = std::remove(data.begin(), data.end(), val);

      if (i == data.end()) {
        return false;
      }

      data.erase(i, data.end());
      bbox = fit_aabb(data.begin(), data.end(), get_aabb);
      return true;
    }
  }

  // remove an object from the rtree.
  template<typename aabb_calc_t>
  bool remove(data_t *val, const aabb_calc_t &get_aabb) {
    return _remove(val, get_aabb(val), get_aabb);
  }

  // functor for ordering nodes by increasing aabb midpoint, along a specified axis.
  struct aabb_cmp_mid {
    size_t dim;
    aabb_cmp_mid(size_t _dim) : dim(_dim) { }

    bool operator()(const data_aabb_t &a, const data_aabb_t &b) const {
      return a.aabb()._m(dim) < b.aabb()._m(dim);
    }
    bool operator()(const rtree_node_t *a, const rtree_node_t *b) const {
      return a->aabb()._m(dim) < b->aabb()._m(dim);
    }
  };

  // functor for ordering nodes by increasing aabb minimum, along a specified axis.
  struct aabb_cmp_min {
    size_t dim;
    aabb_cmp_min(size_t _dim) : dim(_dim) { }

    bool operator()(const data_aabb_t &a, const data_aabb_t &b) const {
      return a.aabb()._l(dim) < b.aabb()._l(dim);
    }
    bool operator()(const rtree_node_t *a, const rtree_node_t *b) const {
      return a->aabb()._l(dim) < b->aabb()._l(dim);
    }
  };

  // functor for ordering nodes by increasing aabb maximum, along a specified axis.
  struct aabb_cmp_max {
    size_t dim;
    aabb_cmp_max(size_t _dim) : dim(_dim) { }

    bool operator()(const data_aabb_t &a, const data_aabb_t &b) const {
      return a.aabb()._h(dim) < b.aabb()._h(dim);
    }
    bool operator()(const rtree_node_t *a, const rtree_node_t *b) const {
      return a->aabb()._h(dim) < b->aabb()._h(dim);
    }
  };

  // facade for projecting node bounding box onto an axis.
  struct aabb_extent {
    size_t dim;
    aabb_extent(size_t _dim) : dim(_dim) { }

    double min(const data_aabb_t &a) const { return a.aabb()._l(dim); }
    double max(const data_aabb_t &a) const { return a.aabb()._h(dim); }
    double len(const data_aabb_t &a) const { return max(a) - min(a); }

    double min(const rtree_node_t *a) const { return a->aabb()._l(dim); }
    double max(const rtree_node_t *a) const { return a->aabb()._h(dim); }
    double len(const rtree_node_t *a) const { return max(a) - min(a); }
  };

  template<typename iter_t>
  static void makeNodes(const iter_t begin,
                        const iter_t end,
                        size_t dim_num,
                        uint32_t dim_mask,
                        size_t child_size,
                        std::vector<rtree_node_t *> &out) {
    const size_t N = std::distance(begin, end);

    size_t dim = 3;
    double r_best = N+1;

    // find the sparsest remaining dimension to partition by.
    for (size_t i = 0; i < 3; ++i) {
      if (dim_mask & (1U << i)) continue;
      aabb_extent extent(i);
      double dmin, dmax, dsum;

      dmin = extent.min(*begin);
      dmax = extent.max(*begin);
      dsum = 0.0;
      for (iter_t j = begin; j != end; ++j) {
        dmin = std::min(dmin, extent.min(*j));
        dmax = std::max(dmax, extent.max(*j));
        dsum += extent.len(*j);
      }
      double r = dsum ? dsum / (dmax - dmin) : 0.0;
      if (r_best > r) {
        dim = i;
        r_best = r;
      }
    }

    assert(dim < 3);

    // dim = dim_num;

    const size_t P = (N + child_size - 1) / child_size;
    const size_t n_parts = (size_t)std::ceil(std::pow((double)P, 1.0 / (3 - dim_num)));

    std::sort(begin, end, aabb_cmp_mid(dim));

    if (dim_num == 2 || n_parts == 1) {
      for (size_t i = 0, s = 0, e = 0; i < P; ++i, s = e) {
        e = N * (i+1) / P;
        assert(e - s <= child_size);
        out.push_back(new rtree_node_t(begin + s, begin + e));
      }
    } else {
      for (size_t i = 0, s = 0, e = 0; i < n_parts; ++i, s = e) {
        e = N * (i+1) / n_parts;
        makeNodes(begin + s, begin + e, dim_num + 1, dim_mask | (1U << dim), child_size, out);
      }
    }
  }

  struct partition_info {
    double score;
    size_t partition_pos;

    partition_info() : score(std::numeric_limits<double>::max()), partition_pos(0) {
    }
    partition_info(double _score, size_t _partition_pos) :
      score(_score),
      partition_pos(_partition_pos) {
    }
  };

  static partition_info findPartition(typename std::vector<data_aabb_t>::iterator base,
                                      std::vector<size_t>::iterator begin,
                                      std::vector<size_t>::iterator end,
                                      size_t part_size) {
    assert(begin < end);

    partition_info best(std::numeric_limits<double>::max(), 0);
    const size_t N = (size_t)std::distance(begin, end);

    std::vector<double> rhs_vol(N, 0.0);

    aabb3d_t rhs = base[begin[N-1]].aabb;
    rhs_vol[N-1] = rhs.volume();
    for (size_t i = N - 1; i > 0; ) {
      rhs.expand(base[begin[--i]].aabb());
      rhs_vol[i] = rhs.volume();
    }

    aabb3d_t lhs = base[begin[0]].aabb;
    for (size_t i = 1; i < N; ++i) {
      lhs.expand(base[begin[i]].aabb());
      if (i % part_size == 0 || (N - i) % part_size == 0) {
        partition_info curr(lhs.volume() + rhs_vol[i], i);
        if (best.score > curr.score) best = curr;
      }
    }
    return best;
  }

  static void partition(typename std::vector<data_aabb_t>::iterator base,
                        std::vector<size_t>::iterator begin,
                        std::vector<size_t>::iterator end,
                        size_t part_size,
                        std::vector<size_t> &part_num,
                        size_t &part_next) {
    assert(begin < end);

    const size_t N = (size_t)std::distance(begin, end);

    partition_info best;
    partition_info curr;
    size_t part_curr = part_num[*begin];

    std::vector<size_t> tmp(begin, end);

    for (size_t dim = 0; dim < 3; ++dim) {
      std::sort(tmp.begin(), tmp.end(), make_index_sort(base, aabb_cmp_min(dim)));
      curr = findPartition(base, tmp.begin(), tmp.end(), part_size);
      if (best.score > curr.score) {
        best = curr;
        std::copy(tmp.begin(), tmp.end(), begin);
      }

      std::sort(tmp.begin(), tmp.end(), make_index_sort(base, aabb_cmp_mid(dim)));
      curr = findPartition(base, tmp.begin(), tmp.end(), part_size);
      if (best.score > curr.score) {
        best = curr;
        std::copy(tmp.begin(), tmp.end(), begin);
      }

      std::sort(tmp.begin(), tmp.end(), make_index_sort(base, aabb_cmp_max(dim)));
      curr = findPartition(base, tmp.begin(), tmp.end(), part_size);
      if (best.score > curr.score) {
        best = curr;
        std::copy(tmp.begin(), tmp.end(), begin);
      }
    }

    for (size_t j = 0; j < best.partition_pos; ++j) part_num[begin[(ssize_t)j]] = part_curr;
    for (size_t j = best.partition_pos; j < N; ++j) part_num[begin[(ssize_t)j]] = part_next;
    ++part_next;

    if (best.partition_pos > part_size) {
      partition(base, begin, begin + best.partition_pos, part_size, part_num, part_next);
    }
    if (N - best.partition_pos > part_size) {
      partition(base, begin + best.partition_pos, end, part_size, part_num, part_next);
    }
  }

  static size_t makePartitions(typename std::vector<data_aabb_t>::iterator begin,
                               typename std::vector<data_aabb_t>::iterator end,
                               size_t part_size,
                               std::vector<size_t> &part_num) {
    const size_t N = std::distance(begin, end);
    std::vector<size_t> idx;
    idx.reserve(N);
    for (size_t i = 0; i < N; ++i) { idx.push_back(i); }
    size_t part_next = 1;

    partition(begin, idx.begin(), idx.end(), part_size, part_num, part_next);
    return part_next;
  }

  static rtree_node_t *construct_STR(std::vector<data_aabb_t> &data,
                                     size_t leaf_size,
                                     size_t internal_size) {
    std::vector<rtree_node_t *> out;
    makeNodes(data.begin(), data.end(), 0, 0, leaf_size, out);

    while (out.size() > 1) {
      std::vector<rtree_node_t *> next;
      makeNodes(out.begin(), out.end(), 0, 0, internal_size, next);
      std::swap(out, next);
    }

    assert(out.size() == 1);
    return out[0];
  }

  template<typename iter_t, typename aabb_calc_t>
  static rtree_node_t *construct_STR(iter_t begin,
                                     iter_t end,
                                     const aabb_calc_t &get_aabb,
                                     size_t leaf_size,
                                     size_t internal_size) {
    std::vector<data_aabb_t> data;
    data.reserve(std::distance(begin, end));
    for (iter_t i = begin; i != end; ++i) {
      data.push_back(data_aabb_t(&*i, get_aabb(*i)));
    }
    return construct_STR(data, leaf_size, internal_size);
  }


  template<typename iter_t, typename aabb_calc_t>
  static rtree_node_t *construct_STR(iter_t begin1,
                                     iter_t end1,
                                     iter_t begin2,
                                     iter_t end2,
                                     const aabb_calc_t &get_aabb,
                                     size_t leaf_size,
                                     size_t internal_size) {
    std::vector<data_aabb_t> data;
    data.reserve(std::distance(begin1, end1) + std::distance(begin2, end2));
    for (iter_t i = begin1; i != end1; ++i) {
      data.push_back(data_aabb_t(&*i, get_aabb(*i)));
    }
    for (iter_t i = begin2; i != end2; ++i) {
      data.push_back(data_aabb_t(&*i, get_aabb(*i)));
    }
    return construct_STR(data, leaf_size, internal_size);
  }

  static rtree_node_t *construct_TGS(typename std::vector<data_aabb_t>::iterator begin,
                                     typename std::vector<data_aabb_t>::iterator end,
                                     size_t leaf_size,
                                     size_t internal_size) {
    size_t N = std::distance(begin, end);

    if (N <= leaf_size) {
      return new rtree_node_t(begin, end);
    } else {
      size_t P = (N + internal_size - 1) / internal_size;
      std::vector<size_t> part_num(N, 0);
      P = makePartitions(begin, end, P, part_num);

      size_t S = 0, E = 0;
      std::vector<rtree_node_t *> children;
      for (size_t i = 0; i < P; ++i) {
        size_t j = S, k = N;
        while (true) {
          while (true) {
            if (j == k) goto done;
            else if (part_num[j] == i) ++j;
            else break;
          }
          --k;
          while (true) {
            if (j == k) goto done;
            else if (part_num[k] != i) --k;
            else break;
          }
          std::swap(*(begin+j), *(begin+k));
          std::swap(part_num[j], part_num[k]);
          ++j;
        }
      done:
        E = j;
        children.push_back(construct_TGS(begin + S, begin + E, leaf_size, internal_size));
        S = E;
      }
      return new rtree_node_t(children.begin(), children.end());
    }
  }

  template<typename iter_t, typename aabb_calc_t>
  static rtree_node_t *construct_TGS(iter_t begin,
                                     iter_t end,
                                     const aabb_calc_t &get_aabb,
                                     size_t leaf_size,
                                     size_t internal_size) {
    std::vector<data_aabb_t> data;
    data.reserve(std::distance(begin, end));
    for (iter_t i = begin; i != end; ++i) {
      data.push_back(data_aabb_t(&*i, get_aabb(*i)));
    }
    return construct_TGS(data.begin(), data.end(), leaf_size, internal_size);
  }

  template<typename iter_t, typename aabb_calc_t>
  static rtree_node_t *construct_TGS(iter_t begin1,
                                     iter_t end1,
                                     iter_t begin2,
                                     iter_t end2,
                                     const aabb_calc_t &get_aabb,
                                     size_t leaf_size,
                                     size_t internal_size) {
    std::vector<data_aabb_t> data;
    data.reserve(std::distance(begin1, end1) + std::distance(begin2, end2));
    for (iter_t i = begin1; i != end1; ++i) {
      data.push_back(data_aabb_t(&*i, get_aabb(*i)));
    }
    for (iter_t i = begin2; i != end2; ++i) {
      data.push_back(data_aabb_t(&*i, get_aabb(*i)));
    }
    return construct_TGS(data.begin(), data.end(), leaf_size, internal_size);
  }
};
