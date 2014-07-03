#include <clemency/geom.hpp>

#include <deque>
#include <list>
#include <vector>
#include <tr1/unordered_map>

#define OFFSET .001
#define DEPTH 7

class dist_t {
public:
  virtual double dist(v3d_t) const =0;
  virtual v3d_t grad(v3d_t) const =0;
};

inline double smin_exp(double a, double b, const double smooth) {
  a *= -smooth;
  b *= -smooth;
  double m = std::max(a, b);
  double res = m + ::log(1.0 + ::exp(std::min(a, b) - m));
  return res / -smooth;
}

class scene_t : public dist_t {
public:
  virtual double dist(v3d_t v) const {
    return smin_exp(
      v3d_t::sub(v, v3d_t::init(-.5,0,.3)).length() - 1.0,
      v3d_t::sub(v, v3d_t::init(0,.9,0)).length() - 1.0,
      30.0
    );
  }
  virtual v3d_t grad(v3d_t v) const {
    const double DEL = 1.0/2048.0;
    const v3d_t basis0 = v3d_t::init(DEL, 0.0, 0.0);
    const v3d_t basis1 = v3d_t::init(0.0, DEL, 0.0);
    const v3d_t basis2 = v3d_t::init(0.0, 0.0, DEL);
    double gx1 = dist(v - basis0);
    double gx2 = dist(v + basis0);
    double gy1 = dist(v - basis1);
    double gy2 = dist(v + basis1);
    double gz1 = dist(v - basis2);
    double gz2 = dist(v + basis2);

    double  gradX,gradY, gradZ;

    gradX = gx2 - gx1;
    gradY = gy2 - gy1;
    gradZ = gz2 - gz1;

    return v3d_t::init(gradX, gradY, gradZ).normalize();
  }
};
dist_t *df = new scene_t;

class isosurface_octree_t : public octree_t<isosurface_octree_t, double> {
public:
  struct funcdata_t {
    v3d_t grad;
    double val;
  };

  struct dual_info_t {
    v3d_t c3_pt;    // volume
    v3d_t c2_pt[6]; // face
    v3d_t c1_pt[8]; // edge
    funcdata_t c3_dat;
    funcdata_t c2_dat[6];
    funcdata_t c1_dat[8];
  };

  typedef octree_t<isosurface_octree_t, double> super;
  dual_info_t *dual_info;
  funcdata_t data[8];
  int depth;

  isosurface_octree_t(const aabb3d_t &_bbox) : super(_bbox), dual_info(NULL), depth(0) {
  }

  ~isosurface_octree_t() {
    if (dual_info) delete dual_info;
  }

  void calculate_dual_info(dist_t &df) {
    if (dual_info == NULL) dual_info = new dual_info_t;
  }
};



struct isosurface_splitter_t : public isosurface_octree_t::visitor_t {
  dist_t &df;

  isosurface_splitter_t(dist_t &_df) : df(_df) {
  }

  void split(isosurface_octree_t *root) {
    root->visit(*this);
  }

  void pre(isosurface_octree_t *node, int depth) {
    node->depth = depth;

    if (depth > DEPTH) return;
    for (size_t i = 0; i < 8; ++i) {
      node->data[i].val = df.dist(node->corner(i));
    }
    if (node->is_leaf()) {
      node->split();
    }
  }

  void post(isosurface_octree_t *node, int depth) {
    if (!node->is_leaf()) {
      bool all_empty = true;
      for (size_t c = 0; all_empty && c < 8; ++c) {
        size_t n_in = 0;
        for (size_t i = 0; i < 8; ++i) {
          if (node->children[c]->data[i].val <= OFFSET) {
            n_in++;
          }
        }
        if (n_in > 0 && n_in < 8) all_empty = false;
      }
      if (all_empty) {
        node->clear_children();
        assert (node->is_leaf());
      }
    }
  }
};



template<typename node_t>
struct cell_t {
  node_t *cell;
  int pos_idx;
  enum cell_type_t {
    VERTEX,
    EDGE,
    FACE,
    VOL
  } cell_type;

  cell_t(node_t *_cell, cell_type_t _cell_type, int _pos_idx) :
    cell(_cell), cell_type(_cell_type), pos_idx(_pos_idx) {
  }
  cell_t() : cell(NULL), cell_type(VOL), pos_idx(0) {
  }

  struct hash {
    size_t operator()(const cell_t &pos) const {
      return (size_t)pos.cell ^ pos.cell_type ^ pos.pos_idx * 4;
    }
    size_t operator()(const std::pair<cell_t, cell_t> &pos_pair) const {
      size_t a = operator()(pos_pair.first);
      size_t b = operator()(pos_pair.second);
      return a ^ (b >> 16 ^ b << 16);
    }
  };

  struct eq {
    bool operator()(const cell_t &a, const cell_t &b) const {
      return a.cell == b.cell && a.cell_type == b.cell_type && a.pos_idx == b.pos_idx;
    }
    bool operator()(const std::pair<cell_t, cell_t> &a,
                    const std::pair<cell_t, cell_t> &b) const {
      return operator()(a.first, b.first) && operator()(a.second, b.second);
    }
  };

  struct lt {
    bool operator()(const cell_t &a, const cell_t &b) {
      if (a.cell == NULL) return false;
      if (b.cell == NULL) return true;

      if (a.cell->depth > b.cell->depth) return true;
      if (a.cell->depth < b.cell->depth) return false;

      return a.cell->bbox.mid < b.cell->bbox.mid;
    }
  };

  template<typename iter_t>
  static cell_t min(iter_t a,
             iter_t b) {
    return *std::min_element(a, b, lt());
  }

  static cell_t min(const cell_t &a,
             const cell_t &b) {
    return lt()(a, b) ? a : b;
  }

  static cell_t min(const cell_t &a,
                    const cell_t &b,
                    const cell_t &c,
                    const cell_t &d) {
    cell_t r = a;
    if (lt()(b, r)) r = b;
    if (lt()(c, r)) r = c;
    if (lt()(d, r)) r = d;
    return r;
  }

  static cell_t min(const cell_t &a,
                    const cell_t &b,
                    const cell_t &c,
                    const cell_t &d,
                    const cell_t &e,
                    const cell_t &f,
                    const cell_t &g,
                    const cell_t &h) {
    cell_t r = a;
    if (lt()(b, r)) r = b;
    if (lt()(c, r)) r = c;
    if (lt()(d, r)) r = d;
    if (lt()(e, r)) r = e;
    if (lt()(f, r)) r = f;
    if (lt()(g, r)) r = g;
    if (lt()(h, r)) r = h;
    return r;
  }


  static void make_edge(node_t *a,
                        node_t *b,
                        node_t *c,
                        node_t *d,
                        int axis,
                        int dir,
                        cell_t cell0,
                        std::back_insert_iterator<std::list<std::vector<cell_t> > > out) {
    if ((a == b && c == d) || (a == c && b == d)) {
      return;
    }

    std::vector<cell_t> tet;
    tet.resize(4);

    cell_t cell1, cell2ab, cell2ac, cell2bd, cell2cd, cell3;

    cell1 = min(cell_t(a, cell_t::EDGE, axis * 4 + 3),
                cell_t(b, cell_t::EDGE, axis * 4 + 2),
                cell_t(c, cell_t::EDGE, axis * 4 + 1),
                cell_t(d, cell_t::EDGE, axis * 4 + 0));

    int a1, a2;

    switch (axis) {
    case 0: a1 = 1; a2 = 2; break;
    case 1: a1 = 2; a2 = 0; break;
    case 2: a1 = 0; a2 = 1; break;
    }

    cell2ab = min(cell_t(a, cell_t::FACE, a1*2 + 1),
                  cell_t(b, cell_t::FACE, a1*2));
    cell2ac = min(cell_t(a, cell_t::FACE, a2*2 + 1),
                  cell_t(c, cell_t::FACE, a2*2));
    cell2bd = min(cell_t(b, cell_t::FACE, a2*2 + 1),
                  cell_t(d, cell_t::FACE, a2*2));
    cell2cd = min(cell_t(c, cell_t::FACE, a1*2 + 1),
                  cell_t(d, cell_t::FACE, a1*2));

    tet[0] = cell0;

    int A, B;
    if (dir == 0) { A = 1; B = 2; } else { A = 2; B = 1; }

    if (a) {
      tet[3] = cell_t(a, cell_t::VOL, 0);
      if (a != b) { tet[B] = cell1; tet[A] = cell2ab; *out++ = tet; }
      if (a != c) { tet[A] = cell1; tet[B] = cell2ac; *out++ = tet; }
    }

    if (b) {
      tet[3] = cell_t(b, cell_t::VOL, 0);
      if (a != b) { tet[A] = cell1; tet[B] = cell2ab; *out++ = tet; }
      if (b != d) { tet[B] = cell1; tet[A] = cell2bd; *out++ = tet; }
    }

    if (c) {
      tet[3] = cell_t(c, cell_t::VOL, 0);
      if (c != d) { tet[A] = cell1; tet[B] = cell2cd; *out++ = tet; }
      if (c != a) { tet[B] = cell1; tet[A] = cell2ac; *out++ = tet; }
    }

    if (d) {
      tet[3] = cell_t(d, cell_t::VOL, 0);
      if (d != c) { tet[B] = cell1; tet[A] = cell2cd; *out++ = tet; }
      if (d != b) { tet[A] = cell1; tet[B] = cell2bd; *out++ = tet; }
    }
  }

  static void generate_tetrahedra(const std::vector<node_t *> &dc,
                                  std::back_insert_iterator<std::list<std::vector<cell_t> > > out) {
    cell_t cell0;

    cell0 =
      min(cell_t(dc[0], cell_t::VERTEX, 7),
          cell_t(dc[1], cell_t::VERTEX, 6),
          cell_t(dc[2], cell_t::VERTEX, 5),
          cell_t(dc[3], cell_t::VERTEX, 4),
          cell_t(dc[4], cell_t::VERTEX, 3),
          cell_t(dc[5], cell_t::VERTEX, 2),
          cell_t(dc[6], cell_t::VERTEX, 1),
          cell_t(dc[7], cell_t::VERTEX, 0));

    make_edge(dc[0], dc[2], dc[4], dc[6], 0, 0, cell0, out); // -x (0,2,4,6)
    make_edge(dc[1], dc[3], dc[5], dc[7], 0, 1, cell0, out); // +x (1,3,5,7)

    make_edge(dc[0], dc[4], dc[1], dc[5], 1, 0, cell0, out); // -y (0,1,4,5)
    make_edge(dc[2], dc[6], dc[3], dc[7], 1, 1, cell0, out); // +y (2,3,6,7)

    make_edge(dc[0], dc[1], dc[2], dc[3], 2, 0, cell0, out); // -z (0,1,2,3)
    make_edge(dc[4], dc[5], dc[6], dc[7], 2, 1, cell0, out); // +z (4,5,6,7)
  }
};



template<typename node_t>
struct result_t {
  typedef cell_t<node_t> cell_type;

  const dist_t &df;

  struct position_t {
    v3d_t pos;
    enum orient_t { IN = -1, ON = 0, OUT = +1 } orient;
    double dist;

    position_t(v3d_t _pos, orient_t _orient, double _dist) :
        pos(_pos), orient(_orient), dist(_dist) {
    }
    position_t() {
    }
  };

  typedef std::tr1::unordered_map<cell_type, position_t, typename cell_type::hash, typename cell_type::eq> tvert_info_map_t;
  tvert_info_map_t tvert_info;

  typedef std::tr1::unordered_map<std::pair<cell_type, cell_type>, size_t, typename cell_type::hash, typename cell_type::eq> int_idx_map_t;
  int_idx_map_t int_idx;
  std::vector<v3d_t> int_pos;

  struct tri_t { int a, b, c; };
  std::vector<tri_t> tris;

  result_t(const dist_t &_df) : df(_df) {
  }

  void dump(std::ostream &out) {
    out <<
      "ply\n" <<
      "format ascii 1.0\n" <<
      "element vertex " << int_pos.size() << "\n" <<
      "property double x\n" <<
      "property double y\n" <<
      "property double z\n" <<
      "element face " << tris.size() << "\n" <<
      "property list int int vertex_indices\n" <<
      "end_header\n";
    for (size_t i = 0; i < int_pos.size(); ++i) {
      out << int_pos[i].x << " " << int_pos[i].y << " " << int_pos[i].z << std::endl;
    }
    for (size_t i = 0; i < tris.size(); ++i) {
      out << "3 " << tris[i].a << " " << tris[i].b <<  " " << tris[i].c << std::endl;
    }
  }

  void binsearch(v3d_t a, v3d_t b, v3d_t &r) {
    double tgt = OFFSET;
    v3d_t test;

    double a_dist = df.dist(a);
    if (a_dist == tgt) { r = a; return; }

    double b_dist = df.dist(b);
    if (b_dist == tgt) { r = b; return; }

    assert(a_dist > tgt);
    assert(b_dist < tgt);

    for (size_t i = 0; i < 10; ++i) {
      double frac = (tgt - b_dist) / (a_dist - b_dist);
      v3d_t c = b + frac * (a - b);
      double c_dist = df.dist(c);
      if (c_dist < tgt) {
        b = c;
        b_dist = c_dist;
      } else if (c_dist == tgt) {
        r = c;
        return;
      } else {
        a = c;
        a_dist = c_dist;
      }
      r = c;
    }
  }

  int get_intersection(cell_type a, cell_type b) {
    position_t pa = get_tvert_info(a);
    position_t pb = get_tvert_info(b);
    assert(pa.orient == position_t::OUT);
    assert(pb.orient != position_t::OUT);
    if (pb.orient == position_t::ON) {
      typename int_idx_map_t::iterator i = int_idx.find(std::make_pair(b, b));
      if (i != int_idx.end()) return i->second;

      int r = int_idx[std::make_pair(b, b)] = int_pos.size();

      int_pos.push_back(pb.pos);

      return r;
    } else {
      typename int_idx_map_t::iterator i = int_idx.find(std::make_pair(a, b));
      if (i != int_idx.end()) return i->second;

      int r = int_idx[std::make_pair(a, b)] = int_pos.size();
      v3d_t intersection;

      binsearch(pa.pos, pb.pos, intersection);
      int_pos.push_back(intersection);

      return r;
    }
  }

  position_t choose_position(cell_type cell) {
    v3d_t lo = cell.cell->corner(0);
    v3d_t hi = cell.cell->corner(7);

    position_t result;
    v3d_t r;

    switch (cell.cell_type) {
    case cell_type::VERTEX: {
      result.pos = cell.cell->corner(cell.pos_idx);
      result.dist = df.dist(result.pos);
      result.orient = result.dist <= OFFSET ? position_t::IN : position_t::OUT;
      break;
    }

    case cell_type::EDGE: {
      r = cell.cell->bbox.mid;
      int a = cell.pos_idx >> 2;
      int b = (a+1)%3;
      int c = (a+2)%3;

      r.v[b] += cell.pos_idx & 1 ? +cell.cell->bbox.extent.v[b] : -cell.cell->bbox.extent.v[b];
      r.v[c] += cell.pos_idx & 2 ? +cell.cell->bbox.extent.v[c] : -cell.cell->bbox.extent.v[c];

      result.pos = r;
      result.dist = df.dist(result.pos);
      result.orient = result.dist <= OFFSET ? position_t::IN : position_t::OUT;
      break;
    }

    case cell_type::FACE: {
      r = cell.cell->bbox.mid;
      int a = cell.pos_idx >> 1;
      int b = (a+1)%3;
      int c = (a+2)%3;

      r.v[a] += cell.pos_idx & 1 ? +cell.cell->bbox.extent.v[a] : -cell.cell->bbox.extent.v[a];

      result.pos = r;
      result.dist = df.dist(result.pos);
      result.orient = result.dist <= OFFSET ? position_t::IN : position_t::OUT;
      break;
    }

    case cell_type::VOL:
      r = cell.cell->bbox.mid;

      result.pos = r;
      result.dist = df.dist(result.pos);
      result.orient = result.dist <= OFFSET ? position_t::IN : position_t::OUT;
      break;
    }

    return result;
  }

  position_t &get_tvert_info(cell_type cell) {
    typename tvert_info_map_t::iterator i = tvert_info.find(cell);
    if (i != tvert_info.end()) return i->second;
    return tvert_info[cell] = choose_position(cell);
  }

  void add_tri(int a, int b, int c) {
    if (a != b && a != c && b != c) {
      tri_t t;
      t.a = a; t.b = b; t.c = c;
      tris.push_back(t);
    }
  }

  void generate_one(const std::vector<cell_type> &tet) {
    int bits = 0;
    int a, b, c, d;

    position_t t[4];
    for (size_t i = 0; i < 4; ++i) {
      t[i] = get_tvert_info(tet[i]);
      if (t[i].orient != position_t::OUT) bits |= 1 << i;
    }

    switch (bits) {
    case 0:  // 0000
      break;

    case 1:  // 0001
      a = get_intersection(tet[3], tet[0]);
      b = get_intersection(tet[2], tet[0]);
      c = get_intersection(tet[1], tet[0]);
      add_tri(a,b,c);
      break;

    case 2:  // 0010
      a = get_intersection(tet[3], tet[1]);
      b = get_intersection(tet[2], tet[1]);
      c = get_intersection(tet[0], tet[1]);
      add_tri(a,c,b);
      break;

    case 3:  // 0011
      a = get_intersection(tet[3], tet[1]);
      b = get_intersection(tet[3], tet[0]);
      c = get_intersection(tet[2], tet[1]);
      d = get_intersection(tet[2], tet[0]);
      add_tri(a,b,c);
      add_tri(c,b,d);
      break;

    case 4:  // 0100
      a = get_intersection(tet[3], tet[2]);
      b = get_intersection(tet[1], tet[2]);
      c = get_intersection(tet[0], tet[2]);
      add_tri(a,b,c);
      break;

    case 5:  // 0101
      a = get_intersection(tet[3], tet[2]);
      b = get_intersection(tet[3], tet[0]);
      c = get_intersection(tet[1], tet[2]);
      d = get_intersection(tet[1], tet[0]);
      add_tri(a,c,b);
      add_tri(b,c,d);
      break;

    case 6:  // 0110
      a = get_intersection(tet[3], tet[2]);
      b = get_intersection(tet[3], tet[1]);
      c = get_intersection(tet[0], tet[2]);
      d = get_intersection(tet[0], tet[1]);
      add_tri(a,b,c);
      add_tri(c,b,d);
      break;

    case 7:  // 0111
      a = get_intersection(tet[3], tet[2]);
      b = get_intersection(tet[3], tet[1]);
      c = get_intersection(tet[3], tet[0]);
      add_tri(a,b,c);
      break;

    case 8:  // 1000
      a = get_intersection(tet[2], tet[3]);
      b = get_intersection(tet[1], tet[3]);
      c = get_intersection(tet[0], tet[3]);
      add_tri(a,c,b);
      break;

    case 9:  // 1001
      a = get_intersection(tet[2], tet[3]);
      b = get_intersection(tet[2], tet[0]);
      c = get_intersection(tet[1], tet[3]);
      d = get_intersection(tet[1], tet[0]);
      add_tri(a,b,c);
      add_tri(c,b,d);
      break;

    case 10: // 1010
      a = get_intersection(tet[2], tet[3]);
      b = get_intersection(tet[2], tet[1]);
      c = get_intersection(tet[0], tet[3]);
      d = get_intersection(tet[0], tet[1]);
      add_tri(a,c,b);
      add_tri(b,c,d);
      break;

    case 11: // 1011
      a = get_intersection(tet[2], tet[3]);
      b = get_intersection(tet[2], tet[1]);
      c = get_intersection(tet[2], tet[0]);
      add_tri(a,c,b);
      break;

    case 12: // 1100
      a = get_intersection(tet[1], tet[3]);
      b = get_intersection(tet[1], tet[2]);
      c = get_intersection(tet[0], tet[3]);
      d = get_intersection(tet[0], tet[2]);
      add_tri(a,b,c);
      add_tri(c,b,d);
      break;

    case 13: // 1101
      a = get_intersection(tet[1], tet[3]);
      b = get_intersection(tet[1], tet[2]);
      c = get_intersection(tet[1], tet[0]);
      add_tri(a,b,c);
      break;

    case 14: // 1110
      a = get_intersection(tet[0], tet[3]);
      b = get_intersection(tet[0], tet[2]);
      c = get_intersection(tet[0], tet[1]);
      add_tri(a,c,b);
      break;

    case 15: // 1111
      break;
    }
  }

  void generate(const std::list<std::vector<cell_type> > &tetrahedra) {
    for (typename std::list<std::vector<cell_type> >::const_iterator i = tetrahedra.begin(); i != tetrahedra.end(); ++i) {
      generate_one(*i);
    }
  }
};

struct leaf_count_t {
  size_t count;
  leaf_count_t() : count(0) {
  }

  template<typename node_t>
  void operator()(const node_t *node, int depth) {
    if (node->is_leaf()) ++count;
  }
};



int main(int argc, char **argv) {
  isosurface_octree_t pt(aabb3d_t::initWithPoints(v3d_t::init(-2.0, -2.0, -2.0), v3d_t::init(+2.0, +2.0, +2.0)));

  isosurface_splitter_t splitter(*df);
  splitter.split(&pt);

  leaf_count_t counter;
  pt.visit_preorder(counter);
  std::cerr << "leaf count after split: " << counter.count << std::endl;

  std::cerr << "generating dual..." << std::endl;
  std::list<std::vector<isosurface_octree_t *> > dual_cells;
  octree_adjacent_cells_t<isosurface_octree_t> dual;
  dual.generate<0>(&pt, std::back_inserter(dual_cells));

  std::cerr << "n(dual cells): " << dual_cells.size() << std::endl;

  std::cerr << "generating tetrahedra..." << std::endl;
  std::list<std::vector<cell_t<isosurface_octree_t> > > tetrahedra;
  int q = 0;
  for (std::list<std::vector<isosurface_octree_t *> >::iterator i = dual_cells.begin(); i != dual_cells.end(); ++i) {
    cell_t<isosurface_octree_t>::generate_tetrahedra(*i, std::back_inserter(tetrahedra));
  }
  std::cerr << "n(tetrahedra): " << tetrahedra.size() << std::endl;

  result_t<isosurface_octree_t> result(*df);
  result.generate(tetrahedra);
  result.dump(std::cout);
}
