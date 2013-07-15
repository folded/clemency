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

#include <gloop/gloop.hpp>
#include <gloop/gloop-model.hpp>


namespace io {
  namespace detail {
    namespace read {
      template<typename mesh_t>
      struct vertex : public gloop::stream::null_reader {
        mesh_t *mesh;

        vertex(mesh_t *_mesh) : mesh(_mesh) {
        }
        virtual void next() {
          mesh->add_vertex(v3d_t::zero());
        }
        virtual void length(int l) {
          if (l > 0) mesh->vertices.reserve(mesh->vertices.size() + l);
        }
        virtual void end() {
        }
        typename mesh_t::vert_t &curr() const {
          return mesh->vertices.back();
        }

        template<int idx>
        struct component : public gloop::stream::reader<double> {
          const vertex *i;
          component(const vertex *_i) : i(_i) { }
          virtual void value(double val) {
            i->curr().pos.v[idx] = val;
          }
        };
      };



      template<typename mesh_t>
      struct face_idx : public gloop::stream::reader<int> {
        mesh_t *mesh;
        mutable std::vector<int> vidx;

        face_idx(mesh_t *_mesh) : mesh(_mesh), vidx() { }

        virtual void length(int l) { vidx.clear(); vidx.reserve(l); }
        virtual void value(int val) { vidx.push_back(val); }
        virtual void end() {
          if (vidx.size() != 3) {
            throw std::runtime_error("Only triangular meshes are supported");
          }
          mesh->add_tri(vidx[0], vidx[1], vidx[2]);
        }
      };



      template<typename mesh_t>
      struct tristrip_idx : public gloop::stream::reader<int> {
        mesh_t *mesh;
        mutable int a, b, c;
        mutable bool clk;
        
        tristrip_idx(mesh_t *_mesh) : mesh(_mesh), a(-1), b(-1), c(-1), clk(true) { }

        virtual void value(int val) {
          a = b; b = c; c = val;
          if (a == -1 || b == -1 || c == -1) {
            clk = true;
          } else {
            if (clk) {
              mesh->add_tri(a, b, c);
            } else {
              mesh->add_tri(c, b, a);
            }
            clk = !clk;
          }
        }

        virtual void length(int len) {
        }
      };
    }



    namespace write {
      template<typename container_t>
      struct vertex : public gloop::stream::null_writer {
        const container_t &cnt;
        int i;

        vertex(const container_t &_cnt) : cnt(_cnt), i(-1) { }

        virtual void next() { ++i; }
        virtual int length() { return cnt.size(); }
        virtual const v3d_t &curr() const {  return cnt[i].pos; }

        template<int idx>
        struct component : public gloop::stream::writer<double> {
          vertex &r;

          component(vertex &_r) : r(_r) { }

          virtual double value() { return r.curr().v[idx]; }
        };
      };


  
      template<typename mesh_t>
      struct face : public gloop::stream::null_writer {
        const mesh_t *mesh;
        const typename mesh_t::tri_t *curr_face;

        face(const mesh_t *_mesh) : mesh(_mesh), curr_face(&_mesh->face_list) { }
        virtual void next() { curr_face = curr_face->next; }
        virtual int length() { return mesh->face_count; }
        const typename mesh_t::tri_t *curr() const { return curr_face; }

        struct idx : public gloop::stream::writer<size_t> {
          face &r;
          gloop::stream::Type data_type;
          const typename mesh_t::tri_t *f;
          int i;

          idx(face &_r, gloop::stream::Type _data_type) :
            r(_r), data_type(_data_type), f(NULL), i(0) {
          }
          virtual void begin() {
            f = r.curr();
            i = 0;
          }
          virtual int length() { return 3; }
          virtual bool isList() { return true; }
          virtual gloop::stream::Type dataType() { return data_type; }
          virtual int maxLength() { return 3; }
          virtual size_t value() {
            return f->v[i++];
          }
        };
      };
    }
  }
 


  template<typename mesh_t, typename filetype_t>
  bool read_mesh(std::istream &in, filetype_t &f, mesh_t *&mesh) {
    mesh = new mesh_t;

    typedef detail::read::vertex<mesh_t> vertex_t;

    vertex_t *vi = new vertex_t(mesh);
    f.addReader("polyhedron.vertex", vi);
    f.addReader("polyhedron.vertex.x", new typename vertex_t::template component<0>(vi));
    f.addReader("polyhedron.vertex.y", new typename vertex_t::template component<1>(vi));
    f.addReader("polyhedron.vertex.z", new typename vertex_t::template component<2>(vi));

    f.addReader("polyhedron.face.vertex_indices", new detail::read::face_idx<mesh_t>(mesh));

    f.addReader("polyhedron.tristrips.vertex_indices", new detail::read::tristrip_idx<mesh_t>(mesh));

    if (!f.read(in)) {
      delete mesh;
      return false;
    }

    return true;
  }



  template<typename mesh_t, typename filetype_t>
  bool write_mesh(std::ostream &out, filetype_t &f, const mesh_t *mesh) {
    f.newBlock("polyhedron");

    typedef detail::write::vertex<std::vector<typename mesh_t::vert_t> > vertex_t;
    typedef detail::write::face<mesh_t> face_t;

    vertex_t *vi = new vertex_t(mesh->vertices);
    f.addWriter("polyhedron.vertex", vi);
    f.addWriter("polyhedron.vertex.x", new typename vertex_t::template component<0>(*vi));
    f.addWriter("polyhedron.vertex.y", new typename vertex_t::template component<1>(*vi));
    f.addWriter("polyhedron.vertex.z", new typename vertex_t::template component<2>(*vi));

    face_t *fi = new face_t(mesh);
    f.addWriter("polyhedron.face", fi);
    f.addWriter("polyhedron.face.vertex_indices",
                new typename face_t::idx(*fi, gloop::stream::smallest_type(mesh->vertices.size())));

    return f.write(out);
  }
}
