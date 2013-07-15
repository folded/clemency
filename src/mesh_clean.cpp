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

#include <clemency/config.h>

#include <utility>

#include <clemency/mesh.hpp>
#include <clemency/mesh_io.hpp>
#include <clemency/rtree.hpp>
#include <clemency/volume_integrals.hpp>

#include <fstream>
#include <deque>
#include <list>
#include <limits>
#include <vector>
#include <ext/algorithm>
#include <set>
#include <map>

#include <tr1/unordered_set>
#include <tr1/unordered_map>

#include <boost/heap/binomial_heap.hpp>
#include <boost/heap/fibonacci_heap.hpp>

#include <boost/program_options.hpp>
#include <boost/scoped_array.hpp>

namespace opt = boost::program_options;

typedef triangle_mesh_t<vert_t, tri_t> mesh_t;

typedef rtree_node_t<mesh_t::tri_t> rtree_t;



struct closest_point_info_t {
  double current_dist2;
  double orientation;
  mesh_t::tri_t *closest_tri;
  v3d_t closest_point;

  void _find(const rtree_t *rt,
             const mesh_t::vlookup_t &vlookup,
             const v3d_t &pt) {

    std::vector<std::pair<double, rtree_t *> > sorted_children;

    for (size_t i = 0; i < rt->data.size(); ++i) {
      mesh_t::tri_t *t = rt->data[i];
      tri3d_t t_pt(vlookup(t->a), vlookup(t->b), vlookup(t->c));
      v3d_t c = t_pt.closest_point(pt);

      double d2 = (c - pt).lengthsq();
      if (d2 < current_dist2) {
        double o = t_pt.orient(pt);
        orientation = o;
        closest_tri = t;
        closest_point = c;
        current_dist2 = d2;
      }
    }

    for (rtree_t *ch = rt->child; ch != NULL; ch = ch->sibling) {
      double bbox_dist2 = ch->bbox.distancesq(pt);
      if (bbox_dist2 > current_dist2) continue;
      sorted_children.push_back(std::make_pair(bbox_dist2, ch));
    }

    std::sort(sorted_children.begin(), sorted_children.end());

    for (size_t i = 0; i < sorted_children.size(); ++i) {
      if (sorted_children[i].first > current_dist2) break;
      _find(sorted_children[i].second, vlookup, pt);
    }
  }

  void find(const rtree_t *rt, const mesh_t::vlookup_t &vlookup, const v3d_t &pt) {
    closest_tri = NULL;
    current_dist2 = std::numeric_limits<double>::max();
    orientation = 0.0;
    _find(rt, vlookup, pt);
  }
};



void process(mesh_t *mesh) {
  rtree_t *rtree = rtree_t::construct_STR(mesh->fbegin(), mesh->fend(), mesh->aabb_getter(), 4, 4);
  std::cerr << "mesh bounds: " << rtree->bbox << std::endl;

  const double STEP = 0.05;
  const double RAD = 0.25;
  size_t tag = ++mesh->curr_tag;
  for (double x = rtree->bbox.xl(); x <= rtree->bbox.xh(); x += STEP) {
    std::cerr << "x=" << x << std::endl;
    for (double y = rtree->bbox.yl(); y <= rtree->bbox.yh(); y += STEP) {
      for (double z = rtree->bbox.zl(); z <= rtree->bbox.zh(); z += STEP) {
        closest_point_info_t info;
        info.find(rtree, mesh->vertex_getter(), v3d_t::init(x, y, z));
        double d = sqrt(info.current_dist2);
        if (d > RAD) {
          info.closest_tri->tag = tag;
          double n_steps = floor((d-RAD) / STEP);
          if (n_steps > 1) z += STEP * (n_steps-1);
        }
      }
    }
  }
  std::vector<mesh_t::tri_t *> tagged_faces;
  for (mesh_t::face_iter iter = mesh->fbegin(); iter != mesh->fend(); ++iter) {
    if (iter->tag == tag) {
      tagged_faces.push_back(&*iter);
    }
  }
  mesh_t *submesh = mesh->submesh(tagged_faces.begin(), tagged_faces.end());

  gloop::ply::PlyWriter file(true, false);
  io::write_mesh(std::cout, file, submesh);
}



int main(int argc, char **argv) {
  std::string src_file;

  opt::options_description generic;
  generic.add_options()
    ("help,h",                                                       "Print help and exit")
    ("version,v",                                                    "Print version string")
    ;

  opt::options_description hidden;
  hidden.add_options()
    ("mesh", opt::value<std::string>(&src_file), "Mesh to load")
    ;

  opt::positional_options_description pos;
  pos.add("mesh", -1);

  opt::options_description cmdline;
  cmdline.add(generic).add(hidden);

  opt::options_description help("Allowed options");
  help.add(generic);

  opt::variables_map vm;
  try {
    opt::store(opt::command_line_parser(argc, argv).options(cmdline).positional(pos).run(), vm);
    opt::notify(vm);
  } catch(opt::error &e) {
    std::cerr << e.what() << std::endl;
    std::cerr << help << "\n";
    exit(1);
  }

  if (vm.count("help")) {
    std::cout << help << "\n";
    exit(0);
  }

  if (vm.count("version")) {
    std::cout << VERSION << "\n";
    exit(0);
  }

  try {
    mesh_t *mesh;

    gloop::ply::PlyReader file;
    if (src_file == "-") {
      io::read_mesh(std::cin, file, mesh);
    } else {
      std::ifstream inf(src_file.c_str());
      if (!inf.good()) {
        std::cerr << "could not open " << src_file << " for reading" << std::endl;
        exit(1);
      }
      io::read_mesh(inf, file, mesh);
    }

    if (!mesh) {
      std::cerr << "could not read mesh from " << src_file << std::endl;
      exit(1);
    }

    process(mesh);

    delete mesh;

  } catch(std::runtime_error e) {
    std::cerr << "oops. something went wrong.\n" << e.what() << "\n";
  }
}
