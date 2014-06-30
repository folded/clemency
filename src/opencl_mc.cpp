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

#include <clemency/cl.hpp>
#include <clemency/geom.hpp>
#include <clemency/mesh.hpp>
#include <clemency/mesh_io.hpp>
#include <clemency/quadric.hpp>

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

#include "mc_table.hpp"

#include <boost/program_options.hpp>
#include <boost/scoped_array.hpp>

namespace opt = boost::program_options;

template<typename mesh_t>
struct marching_cubes_t {
  struct edgeint_t {
    uint16_t a[3], b[3];

    edgeint_t(uint16_t _a[3], uint16_t _b[3]) {
      std::copy(_a, _a+3, a);
      std::copy(_b, _b+3, b);
    }

    edgeint_t(const edgeint_t &o) {
      std::copy(o.a, o.a+3, a);
      std::copy(o.b, o.b+3, b);
    }

    edgeint_t &operator=(const edgeint_t &o) {
      if (this != &o) {
        std::copy(o.a, o.a+3, a);
        std::copy(o.b, o.b+3, b);
      }
      return *this;
    }

    bool operator==(const edgeint_t &o) const {
      return
        a[0] == o.a[0] && a[1] == o.a[1] && a[2] == o.a[2] &&
        b[0] == o.b[0] && b[1] == o.b[1] && b[2] == o.b[2];
    }
  };

  struct edgeint_hash_t {
    size_t operator()(const edgeint_t &e) const {
      size_t r = 0;
      r *= 131; r ^= e.a[0]; r *= 131; r ^= e.a[1]; r *= 131; r ^= e.a[2];
      r *= 131; r ^= e.b[0]; r *= 131; r ^= e.b[1]; r *= 131; r ^= e.b[2];
      return r;
    }
  };

  std::vector<double> px;
  std::vector<double> py;
  std::vector<double> pz;

  typedef std::tr1::unordered_map<edgeint_t, size_t, edgeint_hash_t> edgeint_map_t;
  edgeint_map_t int_idx;

  mesh_t *mesh;

  marching_cubes_t() :
      px(), py(), pz(),
      int_idx(),
      mesh(NULL) {
  }

  ~marching_cubes_t() {
    if (mesh) delete mesh;
  }

  v3d_t position_of_vertex(uint16_t x, uint16_t y, uint16_t z) {
    return v3d_t::init(px[x], py[y], pz[z]);
  }

  edgeint_t decode_edge(uint16_t x, uint16_t y, uint16_t z, uint8_t code) {
    uint16_t a[3];
    uint16_t b[3];

    a[0] = x + ((code >> 4) & 1);
    a[1] = y + ((code >> 5) & 1);
    a[2] = z + ((code >> 6) & 1);
    b[0] = x + ((code >> 0) & 1);
    b[1] = y + ((code >> 1) & 1);
    b[2] = z + ((code >> 2) & 1);

    return edgeint_t(a, b);
  }

  size_t get_intersection(const edgeint_t &e) {
    typename edgeint_map_t::iterator i = int_idx.find(e);

    if (i != int_idx.end()) return (*i).second;
    size_t r = int_idx[e] = mesh->vertices.size();

    mesh->vertices.push_back(typename mesh_t::vert_t(v3d_t::zero()));

    return r;
  }

  void init_mesh() {
    if (mesh) delete mesh;
    mesh = new mesh_t;
  }

  mesh_t *take_mesh() {
    mesh_t *temp = mesh;
    mesh = NULL;
    return temp;
  }

  void grid_eval(cl::program_t &prog, size_t block0, size_t block1, size_t refine) {

    const size_t Nx = px.size();
    const size_t Ny = py.size();
    const size_t Nz = pz.size();

    cl::ctx_t ctx(prog.context());
    std::cerr << "ctx="<<ctx<< std::endl;
    std::vector<cl_device_id> devices = ctx.devices();
    std::cerr << "n_devices=" << devices.size() << std::endl;
    std::cerr << "ctx.device(0)="<<ctx.device(0)<< std::endl;
    cl::queue_t queue(ctx, ctx.device(0));
    std::cerr << "made queue" << std::endl;
    init_mesh();

    {
      cl::kernel_t kern = prog.kernel("grid");

      boost::scoped_array<uint8_t> out(new uint8_t[block0*block0*block0]);
      boost::scoped_array<float> gx(new float[block0]);
      boost::scoped_array<float> gy(new float[block0]);
      boost::scoped_array<float> gz(new float[block0]);

      cl::mem_t gx_mem = ctx.create_buffer(CL_MEM_READ_ONLY, block0 * sizeof(float), NULL);
      cl::mem_t gy_mem = ctx.create_buffer(CL_MEM_READ_ONLY, block0 * sizeof(float), NULL);
      cl::mem_t gz_mem = ctx.create_buffer(CL_MEM_READ_ONLY, block0 * sizeof(float), NULL);

      cl::mem_t out_mem = ctx.create_buffer(CL_MEM_READ_WRITE, block0*block0*block0, NULL);

      const size_t Py = block0;
      const size_t Pz = block0*block0;

      kern.arg(0, gx_mem);
      kern.arg(1, gy_mem);
      kern.arg(2, gz_mem);
      kern.arg(3, out_mem);
      kern.arg(4, Py);
      kern.arg(5, Pz);

      for (size_t _z = 0; _z < Nz; _z += block0-1) {
        const size_t Wz = std::min(block0, Nz - _z);
        for (size_t i = 0; i < Wz; ++i) gz[i] = (float)pz[_z+i];
        queue.sync(cl::cmd_write_buffer(gz_mem, CL_FALSE, 0, Wz * sizeof(float), gz.get()));

        for (size_t _y = 0; _y < Ny; _y += block0-1) {
          const size_t Wy = std::min(block0, Ny - _y);
          for (size_t i = 0; i < Wy; ++i) gy[i] = (float)py[_y+i];
          queue.sync(cl::cmd_write_buffer(gy_mem, CL_FALSE, 0, Wy * sizeof(float), gy.get()));

          for (size_t _x = 0; _x < Nx; _x += block0-1) {
            std::cerr << "block: " << _x << "," << _y << "," << _z << std::endl;
            const size_t Wx = std::min(block0, Nx - _x);
            for (size_t i = 0; i < Wx; ++i) gx[i] = (float)px[_x+i];
            queue.sync(cl::cmd_write_buffer(gx_mem, CL_FALSE, 0, Wx * sizeof(float), gx.get()));

            queue.sync(cl::cmd_task(kern).wrk_sz(Wx, Wy, Wz));

            queue.sync(cl::cmd_read_buffer(out_mem, CL_FALSE, 0, block0*block0*block0, out.get()));

            for (size_t z = 0; z < Wz-1; ++z) {
              for (size_t y = 0; y < Wy-1; ++y) {
                for (size_t x = 0; x < Wx-1; ++x) {
                  size_t bits = 0;

                  if (out[(x+0) + (y+0)*Py + (z+0)*Pz]) bits |= 0x01;
                  if (out[(x+1) + (y+0)*Py + (z+0)*Pz]) bits |= 0x02;
                  if (out[(x+0) + (y+1)*Py + (z+0)*Pz]) bits |= 0x04;
                  if (out[(x+1) + (y+1)*Py + (z+0)*Pz]) bits |= 0x08;
                  if (out[(x+0) + (y+0)*Py + (z+1)*Pz]) bits |= 0x10;
                  if (out[(x+1) + (y+0)*Py + (z+1)*Pz]) bits |= 0x20;
                  if (out[(x+0) + (y+1)*Py + (z+1)*Pz]) bits |= 0x40;
                  if (out[(x+1) + (y+1)*Py + (z+1)*Pz]) bits |= 0x80;

                  const uint8_t *row = mc_table[bits];

                  size_t n_tri = *row++;
                  if (n_tri == 0) continue;

                  {
                    for (size_t t = 0; t < n_tri; ++t) {
                      uint8_t ea, eb, ec;
                      size_t va, vb, vc;
                      ea = *row++;
                      eb = *row++;
                      ec = *row++;
                      
                      va = get_intersection(decode_edge(_x+x, _y+y, _z+z, ea));
                      vb = get_intersection(decode_edge(_x+x, _y+y, _z+z, eb));
                      vc = get_intersection(decode_edge(_x+x, _y+y, _z+z, ec));

                      mesh->add_tri(va, vb, vc);
                    }
                  }
                }
              }
            }
          }
        }
      }
    }

    std::vector<std::pair<edgeint_t, size_t> > vert_calc;
    vert_calc.reserve(int_idx.size());
    for (typename edgeint_map_t::iterator i = int_idx.begin(); i != int_idx.end(); ++i) {
      vert_calc.push_back(std::make_pair((*i).first, (*i).second));
    }

    {
      boost::scoped_array<cl_float3> a(new cl_float3[block1]);
      boost::scoped_array<cl_float3> b(new cl_float3[block1]);
      boost::scoped_array<cl_float3> c(new cl_float3[block1]);

      cl::mem_t a_mem = ctx.create_buffer(CL_MEM_READ_ONLY,  block1 * sizeof(cl_float3), NULL);
      cl::mem_t b_mem = ctx.create_buffer(CL_MEM_READ_ONLY,  block1 * sizeof(cl_float3), NULL);
      cl::mem_t c_mem = ctx.create_buffer(CL_MEM_READ_WRITE, block1 * sizeof(cl_float3), NULL);

      cl::kernel_t kern = prog.kernel("chop");

      kern.arg(0, a_mem);
      kern.arg(1, b_mem);
      kern.arg(2, c_mem);
      kern.arg(3, refine);

      for (size_t i = 0; i < vert_calc.size(); i += block1) {
        const size_t N = std::min(vert_calc.size() - i, block1);

        for (size_t j = 0; j < N; ++j) {
          edgeint_t &e = vert_calc[i+j].first;

          v3d_t a_pos =  position_of_vertex(e.a[0], e.a[1], e.a[2]);
          v3d_t b_pos =  position_of_vertex(e.b[0], e.b[1], e.b[2]);

          a[j].x = (float)a_pos.x; a[j].y = (float)a_pos.y; a[j].z = (float)a_pos.z;
          b[j].x = (float)b_pos.x; b[j].y = (float)b_pos.y; b[j].z = (float)b_pos.z;
        }

        queue.sync(cl::cmd_write_buffer(a_mem, CL_FALSE, 0, N * sizeof(cl_float3), a.get()));
        queue.sync(cl::cmd_write_buffer(b_mem, CL_FALSE, 0, N * sizeof(cl_float3), b.get()));

        queue.sync(cl::cmd_task(kern).wrk_sz(N));

        queue.sync(cl::cmd_read_buffer(c_mem, CL_FALSE, 0, N * sizeof(cl_float3), c.get()));

        for (size_t j = 0; j < N; ++j) {
          mesh->vertices[vert_calc[i+j].second].pos = v3d_t::init(c[j].x, c[j].y, c[j].z);
        }
      }
    }
  }

  void grid_eval(cl::program_t &prog, const aabb3d_t &aabb, const size_t nx, const size_t ny, const size_t nz, size_t block0, size_t block1, size_t refine) {
    assert(nx < 65535);
    assert(ny < 65535);
    assert(nz < 65535);

    interp(aabb.xl(), aabb.xh(), nx, px);
    interp(aabb.yl(), aabb.yh(), ny, py);
    interp(aabb.zl(), aabb.zh(), nz, pz);

    grid_eval(prog, block0, block1, refine);
  }

  void interp(double lo, double hi, size_t n, std::vector<double> &out) {
    out.resize(n+1);

    for (size_t i = 0; i <= n; ++i) {
      double f = (double)i / (double)(n);
      out[i] = (1.0 - f) * lo + f * hi;
    }
  }
};



void validate(boost::any &v,
              const std::vector<std::string> &values,
              aabb3d_t *,
              int)  {
  opt::validators::check_first_occurrence(v);
  const std::string &s = opt::validators::get_single_string(values);

  double xl, yl, zl, xh, yh, zh;
  if (sscanf(s.c_str(), "%lg,%lg,%lg:%lg,%lg,%lg", &xl, &yl, &zl, &xh, &yh, &zh) != 6) {
    throw opt::validation_error(opt::validation_error::invalid_option_value);
  }
  v = boost::any(aabb3d_t::initWithPoints(v3d_t::init(xl,yl,zl), v3d_t::init(xh,yh,zh)));
}

int main(int argc, char **argv) {
  size_t grid, grid_x, grid_y, grid_z;
  size_t block0, block1;
  size_t refine;
  std::string src_file;
  std::string dst_file;
  aabb3d_t aabb;
  aabb3d_t aabb0 = aabb3d_t::initWithPoints(v3d_t::init(-1.0, -1.0, -1.0), v3d_t::init(+1.0, +1.0, +1.0));

  opt::options_description generic;
  generic.add_options()
    ("help,h",                                                            "Print help and exit")
    ("version,v",                                                         "Print version string")
    ("cpu,c",                                                             "Execute on CPU")
    ("list,l",                                                            "List available OpenCL devices")
    ("cl-opts,C", opt::value<std::vector<std::string> >(),                "Set OpenCL compiler option")
    ("grid,g",    opt::value<size_t>(&grid)->default_value(256),          "Size of marching cubes grid")
    ("gx",        opt::value<size_t>(&grid_x)->default_value(0),          "X dimension of mc grid")
    ("gy",        opt::value<size_t>(&grid_y)->default_value(0),          "Y dimension of mc grid")
    ("gz",        opt::value<size_t>(&grid_z)->default_value(0),          "Z dimension of mc grid")
    ("aabb,a",    opt::value<aabb3d_t>(&aabb)->default_value(aabb0),      "AABB of marching cubes grid")
    ("block0",    opt::value<size_t>(&block0)->default_value(256),        "Block size for grid evaluation")
    ("block1",    opt::value<size_t>(&block1)->default_value(1024),       "Block size for intersection refinement")
    ("refine,r",  opt::value<size_t>(&refine)->default_value(10),         "Number of refinement steps")
    ("output,o",  opt::value<std::string>(&dst_file)->default_value("-"), "Ouput file (- for stdout)")
    ;

  opt::options_description hidden;
  hidden.add_options()
    ("kernel", opt::value<std::string>(&src_file), "OpenCL kernel to load")
    ;

  opt::positional_options_description pos;
  pos.add("kernel", -1);

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
    cl::opencl_t cl;
 
    cl.init();

    if (vm.count("list")) {
      cl.list(std::cout);
      exit(0);
    }

    if (!vm.count("kernel")) {
      std::cerr << "no kernel provided" << std::endl;
      exit(1);
    }

    cl::ctx_t ctx;

    if (!vm.count("cpu")) {
      ctx = cl::ctx_t(CL_DEVICE_TYPE_GPU);
    } else {
      ctx = cl::ctx_t(CL_DEVICE_TYPE_CPU);
    }

    std::ifstream inf(src_file.c_str());

    if (!inf.good()) {
      std::cerr << "could not read kernel from " << src_file << std::endl;
      exit(1);
    }

    std::string src;

    while (!inf.eof()) {
      char buf[1024];
      inf.read(buf, 1024);
      src.append(buf, inf.gcount());
    }

    cl::program_t prog = ctx.create_program_from_source(src);

    std::string cl_opts;

    if (vm.count("cl-opts")) {
      std::vector<std::string> cl_opt_strings;
      cl_opt_strings = vm["cl-opts"].as<std::vector<std::string> >();
      for (size_t i = 0; i < cl_opt_strings.size(); ++i) {
        if (i) cl_opts += " ";
        cl_opts += cl_opt_strings[i];
      }
    }

    if (cl_opts.size()) {
      std::cerr << "compiling kernel with options: " << cl_opts << std::endl;
      prog.build(cl_opts);
    } else {
      prog.build();
    }

    if (prog.build_status(ctx.device(0)) != CL_BUILD_SUCCESS) {
      std::cerr << "kernel compilation failed\n\n";
      std::cerr << prog.build_log(ctx.device(0)) << std::endl;
      exit(1);
    }

    typedef triangle_mesh_t<vert_t, tri_t> mesh_t;

    if (!grid_x) grid_x = grid;
    if (!grid_y) grid_y = grid;
    if (!grid_z) grid_z = grid;

    marching_cubes_t<mesh_t> mc;
    mc.grid_eval(prog, aabb, grid_x, grid_y, grid_z, block0, block1, refine);

    mesh_t *mesh = mc.take_mesh();

    if (dst_file == "-") {
      gloop::ply::PlyWriter file(true, false);
      io::write_mesh(std::cout, file, mesh);
    } else {
      std::ofstream outf(dst_file.c_str());
      if (!outf.good()) {
        std::cerr << "could not open " << dst_file << " for writing" << std::endl;
        exit(1);
      }
      gloop::ply::PlyWriter file(false, false);
      io::write_mesh(outf, file, mesh);
    }

    delete mesh;

  } catch(std::runtime_error e) {
    std::cerr << "oops. something went wrong.\n" << e.what() << "\n";
  }
}
