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
#include <algorithm>
#include <set>
#include <map>

#include <unordered_set>
#include <unordered_map>

#include <boost/heap/binomial_heap.hpp>
#include <boost/heap/fibonacci_heap.hpp>

#include <boost/program_options.hpp>
#include <boost/scoped_array.hpp>

namespace opt = boost::program_options;

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
    typedef triangle_mesh_t<vert_t, tri_t> mesh_t;
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

    std::cout << "closed manifold: "
              << (mesh->is_closed_manifold() ? "true" : "false")
              << std::endl;

    delete mesh;

  } catch(std::runtime_error e) {
    std::cerr << "oops. something went wrong.\n" << e.what() << "\n";
  }
}
