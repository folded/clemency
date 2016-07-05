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

#if defined(__APPLE__)
#  include <OpenCL/opencl.h>
#else
#  include <CL/cl.h>
#endif

#include <stdexcept>
#include <vector>
#include <iostream>
#include <sstream>



namespace cl {
  static const char *errstr(cl_int err) {
    switch (err) {
    case CL_SUCCESS:                            return "Success!";
    case CL_DEVICE_NOT_FOUND:                   return "Device not found.";
    case CL_DEVICE_NOT_AVAILABLE:               return "Device not available";
    case CL_COMPILER_NOT_AVAILABLE:             return "Compiler not available";
    case CL_MEM_OBJECT_ALLOCATION_FAILURE:      return "Memory object allocation failure";
    case CL_OUT_OF_RESOURCES:                   return "Out of resources";
    case CL_OUT_OF_HOST_MEMORY:                 return "Out of host memory";
    case CL_PROFILING_INFO_NOT_AVAILABLE:       return "Profiling information not available";
    case CL_MEM_COPY_OVERLAP:                   return "Memory copy overlap";
    case CL_IMAGE_FORMAT_MISMATCH:              return "Image format mismatch";
    case CL_IMAGE_FORMAT_NOT_SUPPORTED:         return "Image format not supported";
    case CL_BUILD_PROGRAM_FAILURE:              return "Program build failure";
    case CL_MAP_FAILURE:                        return "Map failure";
    case CL_INVALID_VALUE:                      return "Invalid value";
    case CL_INVALID_DEVICE_TYPE:                return "Invalid device type";
    case CL_INVALID_PLATFORM:                   return "Invalid platform";
    case CL_INVALID_DEVICE:                     return "Invalid device";
    case CL_INVALID_CONTEXT:                    return "Invalid context";
    case CL_INVALID_QUEUE_PROPERTIES:           return "Invalid queue properties";
    case CL_INVALID_COMMAND_QUEUE:              return "Invalid command queue";
    case CL_INVALID_HOST_PTR:                   return "Invalid host pointer";
    case CL_INVALID_MEM_OBJECT:                 return "Invalid memory object";
    case CL_INVALID_IMAGE_FORMAT_DESCRIPTOR:    return "Invalid image format descriptor";
    case CL_INVALID_IMAGE_SIZE:                 return "Invalid image size";
    case CL_INVALID_SAMPLER:                    return "Invalid sampler";
    case CL_INVALID_BINARY:                     return "Invalid binary";
    case CL_INVALID_BUILD_OPTIONS:              return "Invalid build options";
    case CL_INVALID_PROGRAM:                    return "Invalid program";
    case CL_INVALID_PROGRAM_EXECUTABLE:         return "Invalid program executable";
    case CL_INVALID_KERNEL_NAME:                return "Invalid kernel name";
    case CL_INVALID_KERNEL_DEFINITION:          return "Invalid kernel definition";
    case CL_INVALID_KERNEL:                     return "Invalid kernel";
    case CL_INVALID_ARG_INDEX:                  return "Invalid argument index";
    case CL_INVALID_ARG_VALUE:                  return "Invalid argument value";
    case CL_INVALID_ARG_SIZE:                   return "Invalid argument size";
    case CL_INVALID_KERNEL_ARGS:                return "Invalid kernel arguments";
    case CL_INVALID_WORK_DIMENSION:             return "Invalid work dimension";
    case CL_INVALID_WORK_GROUP_SIZE:            return "Invalid work group size";
    case CL_INVALID_WORK_ITEM_SIZE:             return "Invalid work item size";
    case CL_INVALID_GLOBAL_OFFSET:              return "Invalid global offset";
    case CL_INVALID_EVENT_WAIT_LIST:            return "Invalid event wait list";
    case CL_INVALID_EVENT:                      return "Invalid event";
    case CL_INVALID_OPERATION:                  return "Invalid operation";
    case CL_INVALID_GL_OBJECT:                  return "Invalid OpenGL object";
    case CL_INVALID_BUFFER_SIZE:                return "Invalid buffer size";
    case CL_INVALID_MIP_LEVEL:                  return "Invalid mip-map level";
    default:                                    return "Unknown";
    }
  }

  static inline void __check(cl_int status, const std::string &message = "") {
    if (status != CL_SUCCESS) {
      std::ostringstream msg;
      msg << "OpenCL error " << status << " (" << errstr(status) << ")";
      if (message.size()) msg << " " << message;
      throw std::runtime_error(msg.str());
    }
  }
#define _str(x) #x
#define _xstr(x) _str(x)
#define _check(status) do { __check(status, std::string(__PRETTY_FUNCTION__) + " " + __FILE__ + ":" + _xstr(__LINE__)); } while(0)



  template<typename id_type_t>
  struct info_traits_t {
  };

  template<>
  struct info_traits_t<cl_device_id> {
    typedef cl_device_info info_type;
    typedef cl_int (*info_func_t)(cl_device_id, cl_device_info, size_t, void *, size_t *);

    static inline info_func_t func() { return clGetDeviceInfo; }
  };

  template<>
  struct info_traits_t<cl_platform_id> {
    typedef cl_platform_info info_type;
    typedef cl_int (*info_func_t)(cl_platform_id, cl_platform_info, size_t, void *, size_t *);

    static inline info_func_t func() { return clGetPlatformInfo; }
  };

  template<>
  struct info_traits_t<cl_context> {
    typedef cl_context_info info_type;
    typedef cl_int (*info_func_t)(cl_context, cl_context_info, size_t, void *, size_t *);

    static inline info_func_t func() { return clGetContextInfo; }
  };

  template<>
  struct info_traits_t<cl_command_queue> {
    typedef cl_command_queue_info info_type;
    typedef cl_int (*info_func_t)(cl_command_queue, cl_command_queue_info, size_t, void *, size_t *);

    static inline info_func_t func() { return clGetCommandQueueInfo; }
  };

  template<>
  struct info_traits_t<cl_program> {
    typedef cl_program_info info_type;
    typedef cl_int (*info_func_t)(cl_program, cl_program_info, size_t, void *, size_t *);

    static inline info_func_t func() { return clGetProgramInfo; }
  };

  template<>
  struct info_traits_t<cl_kernel> {
    typedef cl_kernel_info info_type;
    typedef cl_int (*info_func_t)(cl_kernel, cl_kernel_info, size_t, void *, size_t *);

    static inline info_func_t func() { return clGetKernelInfo; }
  };

  template<>
  struct info_traits_t<cl_event> {
    typedef cl_event_info info_type;
    typedef cl_int (*info_func_t)(cl_event, cl_event_info, size_t, void *, size_t *);

    static inline info_func_t func() { return clGetEventInfo; }
  };

  template<>
  struct info_traits_t<cl_mem> {
    typedef cl_mem_info info_type;
    typedef cl_int (*info_func_t)(cl_mem, cl_mem_info, size_t, void *, size_t *);

    static inline info_func_t func() { return clGetMemObjectInfo; }
  };

  template<>
  struct info_traits_t<std::pair<cl_program, cl_device_id> > {
    typedef cl_program_build_info info_type;
    typedef cl_int (*info_func_t)(std::pair<cl_program, cl_device_id>, cl_program_build_info, size_t, void *, size_t *);

    static inline cl_int _clGetProgramBuildInfo(std::pair<cl_program, cl_device_id> id, cl_program_build_info info, size_t param_sz, void *param_val, size_t *param_sz_ret) {
      return clGetProgramBuildInfo(id.first, id.second, info, param_sz, param_val, param_sz_ret);
    }

    static inline info_func_t func() { return _clGetProgramBuildInfo; }
  };



  template<typename id_type_t>
  struct info_t {
    static inline cl_int _call(id_type_t id,
                               typename info_traits_t<id_type_t>::info_type info,
                               size_t param_sz,
                               void *param_val,
                               size_t *param_sz_ret) {
      return info_traits_t<id_type_t>::func()(id, info, param_sz, param_val, param_sz_ret);
    }

    static inline std::string getstr(id_type_t id,
                                     typename info_traits_t<id_type_t>::info_type info) {
      size_t len;
      _check(_call(id, info, 0, NULL, &len));
      char *buf = new char[len+1];
      buf[len] = 0;
      _check(_call(id, info, len, buf, NULL));
      std::string ret(buf);
      delete [] buf;
      return ret;
    };

    template<typename T>
    static inline void getvec(id_type_t id,
                              typename info_traits_t<id_type_t>::info_type info,
                              T *buf,
                              size_t n_buf) {
      _check(_call(id, info, sizeof(T) * n_buf, (void *)buf, NULL));
    }

    template<typename T>
    static inline T get(id_type_t id,
                        typename info_traits_t<id_type_t>::info_type info) {
      T ret;
      getvec(id, info, &ret, 1);
      return ret;
    }
  };



  template<typename derived>
  struct cmd_rect {
    size_t v_src_origin[3];
    size_t v_dst_origin[3];
    size_t v_region[3];
    size_t src_row_pitch;
    size_t src_slice_pitch;
    size_t dst_row_pitch;
    size_t dst_slice_pitch;

    cmd_rect() {
      v_src_origin[0] = v_src_origin[1] = v_src_origin[2] = 0;
      v_dst_origin[0] = v_dst_origin[1] = v_dst_origin[2] = 0;
      v_region[0] = v_region[1] = v_region[2] = 0;
      src_row_pitch = 0;
      src_slice_pitch = 0;
      dst_row_pitch = 0;
      dst_slice_pitch = 0;
    }

    derived &rect(size_t x, size_t y, size_t z = 0) { v_region[0] = x; v_region[1] = y; v_region[2] = z; return *this; }
    derived &src(size_t x, size_t y, size_t z = 0) { v_src_origin[0] = x; v_src_origin[1] = y; v_src_origin[2] = z; return *this; }
    derived &dst(size_t x, size_t y, size_t z = 0) { v_dst_origin[0] = x; v_dst_origin[1] = y; v_dst_origin[2] = z; return *this; }
    derived &src_pitch(size_t x, size_t y = 0) { src_row_pitch = x; src_slice_pitch = y; return *this; }
    derived &dst_pitch(size_t x, size_t y = 0) { dst_row_pitch = x; dst_slice_pitch = y; return *this; }
  };



  struct cmd_read_rect : public cmd_rect<cmd_read_rect> {
    typedef cmd_rect<cmd_read_rect> super;
    cl_mem buffer;
    cl_bool blocking;
    void *ptr;

    cmd_read_rect(cl_mem _buffer, cl_bool _blocking, void *_ptr) : 
      super(), buffer(_buffer), blocking(_blocking), ptr(_ptr) {
    }

    cl_int operator()(cl_command_queue queue, size_t n_wait, const cl_event *wait, cl_event *event) const {
      return clEnqueueReadBufferRect(queue,
                                     buffer,
                                     blocking,
                                     v_src_origin, v_dst_origin, v_region,
                                     src_row_pitch, src_slice_pitch,
                                     dst_row_pitch, dst_slice_pitch,
                                     ptr, n_wait, wait, event);
    }
  };



  struct cmd_write_rect : public cmd_rect<cmd_write_rect> {
    typedef cmd_rect<cmd_write_rect> super;
    cl_mem buffer;
    cl_bool blocking;
    const void *ptr;

    cmd_write_rect(cl_mem _buffer, cl_bool _blocking, const void *_ptr) :
      super(), buffer(_buffer), blocking(_blocking), ptr(_ptr) {
    }

    cl_int operator()(cl_command_queue queue, size_t n_wait, const cl_event *wait, cl_event *event) const {
      return clEnqueueWriteBufferRect(queue,
                                      buffer,
                                      blocking,
                                      v_dst_origin, v_src_origin, v_region,
                                      dst_row_pitch, dst_slice_pitch,
                                      src_row_pitch, src_slice_pitch,
                                      ptr, n_wait, wait, event);
    }
  };



  struct cmd_copy_rect : public cmd_rect<cmd_copy_rect> {
    typedef cmd_rect<cmd_copy_rect> super;
    cl_mem src_buffer;
    cl_mem dst_buffer;

    cmd_copy_rect(cl_mem _src_buffer, cl_mem _dst_buffer) :
      super(), src_buffer(_src_buffer), dst_buffer(_dst_buffer) {
    }

    cl_int operator()(cl_command_queue queue, size_t n_wait, const cl_event *wait, cl_event *event) const {
      return clEnqueueCopyBufferRect(queue,
                                     src_buffer, dst_buffer,
                                     v_src_origin, v_dst_origin, v_region,
                                     src_row_pitch, src_slice_pitch,
                                     dst_row_pitch, dst_slice_pitch,
                                     n_wait, wait, event);
    }
  };



  struct cmd_read_buffer {
    cl_mem buffer;
    cl_bool blocking;
    size_t offset;
    size_t size;
    void *ptr;

    cmd_read_buffer(cl_mem _buffer,
                    cl_bool _blocking,
                    size_t _offset,
                    size_t _size,
                    void *_ptr) :
      buffer(_buffer), blocking(_blocking), offset(_offset), size(_size), ptr(_ptr) {
    }

    cl_int operator()(cl_command_queue queue, size_t n_wait, const cl_event *wait, cl_event *event) const {
      return clEnqueueReadBuffer(queue, buffer, blocking, offset, size, ptr, n_wait, wait, event);
    }
  };



  struct cmd_write_buffer {
    cl_mem buffer;
    cl_bool blocking;
    size_t offset;
    size_t size;
    const void *ptr;

    cmd_write_buffer(cl_mem _buffer,
                     cl_bool _blocking,
                     size_t _offset,
                     size_t _size,
                     const void *_ptr) :
      buffer(_buffer), blocking(_blocking), offset(_offset), size(_size), ptr(_ptr) {
    }

    cl_int operator()(cl_command_queue queue, size_t n_wait, const cl_event *wait, cl_event *event) const {
      return clEnqueueWriteBuffer(queue, buffer, blocking, offset, size, ptr, n_wait, wait, event);
    }
  };



  struct cmd_copy_buffer {
    cl_mem src_buffer;
    cl_mem dst_buffer;
    size_t cb;
    size_t src_offset;
    size_t dst_offset;

    cmd_copy_buffer(cl_mem _src_buffer,
                    cl_mem _dst_buffer,
                    size_t _cb,
                    size_t _src_offset = 0,
                    size_t _dst_offset = 0) :
      src_buffer(_src_buffer), dst_buffer(_dst_buffer),
      cb(_cb),
      src_offset(_src_offset), dst_offset(_dst_offset) {
    }

    cl_int operator()(cl_command_queue queue, size_t n_wait, const cl_event *wait, cl_event *event) const {
      return clEnqueueCopyBuffer(queue, src_buffer, dst_buffer, src_offset, dst_offset, cb, n_wait, wait, event);
    }
  };



  struct cmd_task {
    cl_kernel kernel;
    std::vector<size_t> v_offset;
    std::vector<size_t> v_wrk_sz;
    std::vector<size_t> v_grp_sz;

    cmd_task(cl_kernel _kernel) : kernel(_kernel), v_offset(), v_wrk_sz(), v_grp_sz() {
    }

    cmd_task &offset(size_t x)                     { v_offset.resize(1); v_offset[0] = x; return *this; }
    cmd_task &offset(size_t x, size_t y)           { v_offset.resize(2); v_offset[0] = x; v_offset[1] = y; return *this; }
    cmd_task &offset(size_t x, size_t y, size_t z) { v_offset.resize(3); v_offset[0] = x; v_offset[1] = y; v_offset[2] = z;  return *this; }
    cmd_task &offset(const std::vector<size_t> &v) { v_offset = v; return *this; }

    cmd_task &wrk_sz(size_t x)                     { v_wrk_sz.resize(1); v_wrk_sz[0] = x; return *this; }
    cmd_task &wrk_sz(size_t x, size_t y)           { v_wrk_sz.resize(2); v_wrk_sz[0] = x; v_wrk_sz[1] = y; return *this; }
    cmd_task &wrk_sz(size_t x, size_t y, size_t z) { v_wrk_sz.resize(3); v_wrk_sz[0] = x; v_wrk_sz[1] = y; v_wrk_sz[2] = z; return *this; }
    cmd_task &wrk_sz(const std::vector<size_t> &v) { v_wrk_sz = v; return *this; }

    cmd_task &grp_sz(size_t x)                     { v_grp_sz.resize(1); v_grp_sz[0] = x; return *this; }
    cmd_task &grp_sz(size_t x, size_t y)           { v_grp_sz.resize(2); v_grp_sz[0] = x; v_grp_sz[1] = y; return *this; }
    cmd_task &grp_sz(size_t x, size_t y, size_t z) { v_grp_sz.resize(3); v_grp_sz[0] = x; v_grp_sz[1] = y; v_grp_sz[2] = z; return *this; }
    cmd_task &grp_sz(const std::vector<size_t> &v) { v_grp_sz = v; return *this; }

    cl_int operator()(cl_command_queue queue, size_t n_wait, const cl_event *wait, cl_event *event) const {
      if (!v_wrk_sz.size()) {
        return clEnqueueTask(queue, kernel, n_wait, wait, event);
      } else {
        size_t n_dim = v_wrk_sz.size();
        const size_t   *global_work_size = &v_wrk_sz[0];
        const size_t *global_work_offset = (v_offset.size() == n_dim) ? &v_offset[0] : NULL;
        const size_t    *local_work_size = (v_grp_sz.size() == n_dim) ? &v_grp_sz[0] : NULL;
        return clEnqueueNDRangeKernel(queue, kernel, n_dim, global_work_offset, global_work_size, local_work_size, n_wait, wait, event);
      }
    }
  };



  struct dev_t {
    typedef info_t<cl_device_id> info_type;

    cl_device_id id;

    dev_t(cl_device_id _id) : id(_id) {
    }

    cl_device_type               type() const { return info_type::get<cl_device_type>(id, CL_DEVICE_TYPE);                   }
    cl_uint                 vendor_id() const { return info_type::get<cl_uint>(id, CL_DEVICE_VENDOR_ID);                     }
    cl_uint         max_compute_units() const { return info_type::get<cl_uint>(id, CL_DEVICE_MAX_COMPUTE_UNITS);             }
    cl_uint          max_workitem_dim() const { return info_type::get<cl_uint>(id, CL_DEVICE_MAX_WORK_ITEM_DIMENSIONS);      }
    size_t         max_workgroup_size() const { return info_type::get<size_t>(id, CL_DEVICE_MAX_WORK_GROUP_SIZE);            }
    cl_uint         max_clockfreq_mhz() const { return info_type::get<cl_uint>(id, CL_DEVICE_MAX_CLOCK_FREQUENCY);           }

    cl_uint           global_mem_size() const { return info_type::get<cl_uint>(id, CL_DEVICE_GLOBAL_MEM_SIZE);               }
    cl_uint            local_mem_size() const { return info_type::get<cl_uint>(id, CL_DEVICE_LOCAL_MEM_SIZE);                }
    cl_uint                 max_alloc() const { return info_type::get<cl_uint>(id, CL_DEVICE_MAX_MEM_ALLOC_SIZE);            }

    cl_bool         has_image_support() const { return info_type::get<cl_bool>(id, CL_DEVICE_IMAGE_SUPPORT);                 }
    cl_bool          device_available() const { return info_type::get<cl_bool>(id, CL_DEVICE_AVAILABLE);                     }
    cl_bool        compiler_available() const { return info_type::get<cl_bool>(id, CL_DEVICE_COMPILER_AVAILABLE);            }

    size_t         max_parameter_size() const { return info_type::get<size_t>(id, CL_DEVICE_MAX_PARAMETER_SIZE);             }
    cl_uint         max_constant_args() const { return info_type::get<cl_uint>(id, CL_DEVICE_MAX_CONSTANT_ARGS);             }

    std::string                  name() const { return info_type::getstr(id, CL_DEVICE_NAME);                                }
    std::string                vendor() const { return info_type::getstr(id, CL_DEVICE_VENDOR);                              }
    std::string               version() const { return info_type::getstr(id, CL_DEVICE_VERSION);                             }
    std::string             c_version() const { return info_type::getstr(id, CL_DEVICE_OPENCL_C_VERSION);                    }
    std::string        driver_version() const { return info_type::getstr(id, CL_DRIVER_VERSION);                             }
    std::string               profile() const { return info_type::getstr(id, CL_DEVICE_PROFILE);                             }
    std::string            extensions() const { return info_type::getstr(id, CL_DEVICE_EXTENSIONS);                          }

    bool is_cpu() const { return type() & CL_DEVICE_TYPE_CPU; }
    bool is_gpu() const { return type() & CL_DEVICE_TYPE_GPU; }

    std::vector<cl_uint> max_workitem_sizes() const {
      std::vector<cl_uint> sizes;
      sizes.resize(max_workitem_dim());
      info_type::getvec<cl_uint>(id, CL_DEVICE_MAX_WORK_ITEM_SIZES, &sizes[0], sizes.size());
      return sizes;
    }

    operator cl_device_id() const { return id; }
  };



  struct platform_t {
    typedef info_t<cl_platform_id> info_type;

    cl_platform_id id;
    std::vector<dev_t *> cl_devices;

    platform_t(cl_platform_id _id) : id(_id) {
      cl_uint n_dev;
      clGetDeviceIDs(id, CL_DEVICE_TYPE_ALL, 0, NULL, &n_dev);
      cl_devices.reserve(n_dev);
      cl_device_id *dev;
      dev = new cl_device_id[n_dev];
      _check(clGetDeviceIDs(id, CL_DEVICE_TYPE_ALL, n_dev, dev, NULL));
      for (cl_uint i = 0; i < n_dev; ++i) cl_devices.push_back(new dev_t(dev[i]));
      delete [] dev;
    }

    ~platform_t() {
      for (size_t i = 0; i < cl_devices.size(); ++i) {
        delete cl_devices[i];
      }
    }

    std::string    profile() const { return info_type::getstr(id, CL_PLATFORM_PROFILE);    }
    std::string    version() const { return info_type::getstr(id, CL_PLATFORM_VERSION);    }
    std::string       name() const { return info_type::getstr(id, CL_PLATFORM_NAME);       }
    std::string     vendor() const { return info_type::getstr(id, CL_PLATFORM_VENDOR);     }
    std::string extensions() const { return info_type::getstr(id, CL_PLATFORM_EXTENSIONS); }

    operator cl_platform_id() const { return id; }
  };



  struct event_t {
    typedef info_t<cl_event> info_type;

    cl_event id;

    event_t &operator=(const event_t &event) {
      if (this != &event) {
        if (event.id) {
          clRetainEvent(event.id);
        }
        if (id) {
          clReleaseEvent(id);
        }
        id = event.id;
      }
      return *this;
    }

    event_t() : id(0) {
    }

    event_t(const event_t &event) : id(0) {
      *this = event;
    }

    event_t(cl_event _id) : id(_id) {
      if (id) {
        clRetainEvent(id);
      }
    }

    ~event_t() {
      if (id) {
        clReleaseEvent(id);
      }
    }

    cl_uint         refcount() const { return info_type::get<cl_uint>(id, CL_EVENT_REFERENCE_COUNT);         }
    cl_context       context() const { return info_type::get<cl_context>(id, CL_EVENT_CONTEXT);              }
    cl_command_queue   queue() const { return info_type::get<cl_command_queue>(id, CL_EVENT_COMMAND_QUEUE);  }
    cl_command_type     type() const { return info_type::get<cl_command_type>(id, CL_EVENT_COMMAND_TYPE);    }
    cl_int            status() const { return info_type::get<cl_int>(id, CL_EVENT_COMMAND_EXECUTION_STATUS); }

    void set_status(cl_int status) const {
      _check(clSetUserEventStatus(id, status));
    }

    void wait() const {
      cl_int result = clWaitForEvents(1, &id);
      if (result == CL_EXEC_STATUS_ERROR_FOR_EVENTS_IN_WAIT_LIST) {
        _check(status());
      }
      _check(result);
    }

    operator cl_event() const { return id; }
  };



  struct mem_t {
    typedef info_t<cl_mem> info_type;

    cl_mem id;

    mem_t &operator=(const mem_t &mem) {
      if (this != &mem) {
        if (mem.id) clRetainMemObject(mem.id);
        if (id) clReleaseMemObject(id);
        id = mem.id;
      }
      return *this;
    }

    mem_t() : id(0) {
    }

    mem_t(const mem_t &mem) : id(0) {
      *this = mem;
    }

    mem_t(cl_mem _id) : id(_id) {
      clRetainMemObject(id);
    }

    ~mem_t() {
      if (id) clReleaseMemObject(id);
    }

    cl_uint         refcount() const { return info_type::get<cl_uint>(id, CL_MEM_REFERENCE_COUNT);     }
    cl_mem_object_type  type() const { return info_type::get<cl_mem_object_type>(id, CL_MEM_TYPE);     }
    cl_mem_flags       flags() const { return info_type::get<cl_mem_flags>(id, CL_MEM_FLAGS);          }
    size_t              size() const { return info_type::get<size_t>(id, CL_MEM_SIZE);                 }
    void *          host_ptr() const { return info_type::get<void *>(id, CL_MEM_HOST_PTR);             }
    cl_uint         mapcount() const { return info_type::get<cl_uint>(id, CL_MEM_MAP_COUNT);           }
    cl_mem            parent() const { return info_type::get<cl_mem>(id, CL_MEM_ASSOCIATED_MEMOBJECT); }
    size_t            offset() const { return info_type::get<size_t>(id, CL_MEM_OFFSET);               }

    operator cl_mem() const { return id; }
  };



  struct kernel_t {
    typedef info_t<cl_kernel> info_type;

    cl_kernel id;

    kernel_t &operator=(const kernel_t &kern) {
      if (this != &kern) {
        if (kern.id) clRetainKernel(kern.id);
        if (id) clReleaseKernel(id);
        id = kern.id;
      }
      return *this;
    }

    kernel_t() : id(0) {
    }

    kernel_t(const kernel_t &kern) : id(0) {
      *this = kern;
    }

    kernel_t(cl_kernel _id) : id(_id) {
      clRetainKernel(id);
    }

    ~kernel_t() {
      if (id) clReleaseKernel(id);
    }

    cl_uint         refcount() const { return info_type::get<cl_uint>(id, CL_KERNEL_REFERENCE_COUNT); }
    std::string         name() const { return info_type::getstr(id, CL_KERNEL_FUNCTION_NAME);         }
    cl_uint         num_args() const { return info_type::get<cl_uint>(id, CL_KERNEL_NUM_ARGS);        }
    cl_context       context() const { return info_type::get<cl_context>(id, CL_KERNEL_CONTEXT);      }
    cl_program       program() const { return info_type::get<cl_program>(id, CL_KERNEL_PROGRAM);      }

    template<typename T>
    void arg(int arg_num, T val) const {
      _check(clSetKernelArg(id, arg_num, sizeof(T), (void *)&val));
    }

    operator cl_kernel() const { return id; }
  };



  struct program_t {
    typedef info_t<std::pair<cl_program, cl_device_id> > build_info_type;
    typedef info_t<cl_program> info_type;

    cl_program id;

    program_t &operator=(const program_t &prog) {
      if (this != &prog) {
        if (prog.id) clRetainProgram(prog.id);
        if (id) clReleaseProgram(id);
        id = prog.id;
      }
      return *this;
    }

    program_t() : id(0) {
    }

    program_t(const program_t &prog) : id(0) {
      *this = prog;
    }

    program_t(cl_program _id) : id(_id) {
      clRetainProgram(id);
    }

    ~program_t() {
      if (id) clReleaseProgram(id);
    }

    cl_uint         refcount() const { return info_type::get<cl_uint>(id, CL_CONTEXT_REFERENCE_COUNT); }
    cl_context       context() const { return info_type::get<cl_context>(id, CL_PROGRAM_CONTEXT);      }
    std::string       source() const { return info_type::getstr(id, CL_PROGRAM_SOURCE);                }
    cl_uint      num_devices() const { return info_type::get<cl_uint>(id, CL_CONTEXT_NUM_DEVICES);     }

    cl_build_status build_status(cl_device_id dev) const {
      return build_info_type::get<cl_build_status>(std::make_pair(id, dev), CL_PROGRAM_BUILD_STATUS);
    }

    std::string build_options(cl_device_id dev) const {
      return build_info_type::getstr(std::make_pair(id, dev), CL_PROGRAM_BUILD_OPTIONS);
    }

    std::string build_log(cl_device_id dev) const {
      return build_info_type::getstr(std::make_pair(id, dev), CL_PROGRAM_BUILD_LOG);
    }

    std::vector<cl_device_id> devices() const {
      std::vector<cl_device_id> dev;
      dev.resize(num_devices());
      info_type::getvec<cl_device_id>(id, CL_CONTEXT_DEVICES, &dev[0], dev.size());
      return dev;
    }

    cl_device_id device(int idx) const {
      return devices()[idx];
    }

    bool build(size_t n_devices, cl_device_id *devices, const std::string &opts = "") {
      cl_int err = clBuildProgram(id, n_devices, devices, opts.c_str(), NULL, NULL);
      switch (err) {
      case CL_SUCCESS:               return true;
      case CL_BUILD_PROGRAM_FAILURE: return false;
      default:                       _check(err); return false;
      }
    }

    bool build(const std::string &opts = "") {
      return build(0, NULL, opts);
    }

    bool build(cl_device_id dev, const std::string &opts = "") {
      return build(1, &dev, opts);
    }

    kernel_t kernel(const std::string &kern) const {
      cl_int err;
      kernel_t k(clCreateKernel(id, kern.c_str(), &err));
      _check(err);
      return k;
    }

    operator cl_program() const { return id; }
  };



  struct ctx_t {
    typedef info_t<cl_context> info_type;

    cl_context id;

    void notify(const char *errinfo,
                const void *private_info,
                size_t cb) {
      //throw std::runtime_error(errinfo);
    }

    static void _notify(const char *errinfo,
                        const void *private_info,
                        size_t cb,
                        void *user_data) {
      static_cast<ctx_t *>(user_data)->notify(errinfo, private_info, cb);
    }

    ctx_t &operator=(const ctx_t &ctx) {
      if (this != &ctx) {
        if (ctx.id) clRetainContext(ctx.id);
        if (id) clReleaseContext(id);
        id = ctx.id;
      }
      return *this;
    }

    ctx_t() : id(0) {
    }

    ctx_t(cl_device_type device_type, platform_t *plat = NULL) {
      cl_int err;
      if (plat) {
        cl_context_properties props[] = { CL_CONTEXT_PLATFORM, (cl_context_properties)plat->id, 0 };
        id = clCreateContextFromType(props, device_type, _notify, (void *)this, &err);
      } else {
        id = clCreateContextFromType(NULL, device_type, _notify, (void *)this, &err);
      }
      _check(err);
    }

    ctx_t(const ctx_t &ctx) : id(0) {
      *this = ctx;
    }

    ctx_t(cl_context _id) : id(_id) {
      clRetainContext(id);
    }

    ctx_t(dev_t *dev, platform_t *plat = NULL) {
      cl_int err;
      if (plat) {
        cl_context_properties props[] = { CL_CONTEXT_PLATFORM, (cl_context_properties)plat->id, 0 };
        id = clCreateContext(props, 1, &dev->id, _notify, (void *)this, &err);
      } else {
        id = clCreateContext(NULL, 1, &dev->id, _notify, (void *)this, &err);
      }
      _check(err);
    }

    ~ctx_t() {
      if (id) clReleaseContext(id);
    }

    cl_uint    refcount() const { return info_type::get<cl_uint>(id, CL_CONTEXT_REFERENCE_COUNT); }
    cl_uint num_devices() const { return info_type::get<cl_uint>(id, CL_CONTEXT_NUM_DEVICES);     }

    std::vector<cl_device_id> devices() const {
      std::vector<cl_device_id> dev;
      dev.resize(num_devices());
      info_type::getvec<cl_device_id>(id, CL_CONTEXT_DEVICES, &dev[0], dev.size());
      return dev;
    }

    cl_device_id device(int idx) const {
      return devices()[idx];
    }

    cl_mem create_buffer(cl_mem_flags flags,
                         size_t size,
                         void *host_ptr) {
      cl_int err;
      cl_mem ret = clCreateBuffer(id, flags, size, host_ptr, &err);

      _check(err);
      return ret;
    }

    program_t create_program_from_source(cl_uint n_ch, const char **ch) {
      cl_int err;
      cl_program ret = clCreateProgramWithSource(id, n_ch, ch, NULL, &err);
      _check(err);
      return program_t(ret);
    }

    program_t create_program_from_source(const std::string &str) {
      const char *ch[1] = { str.c_str() };
      return create_program_from_source(1, ch);
    }

    program_t create_program_from_source(const char *str) {
      const char *ch[1] = { str };
      return create_program_from_source(1, ch);
    }

    event_t create_event() const {
      cl_int err;
      event_t event(clCreateUserEvent(id, &err));
      _check(err);
      return event;
    }

    operator cl_context() const { return id; }
  };



  struct queue_t {
    typedef info_t<cl_command_queue> info_type;

    cl_command_queue id;

    queue_t &operator=(const queue_t &q) {
      if (this != &q) {
        if (q.id) clRetainCommandQueue(q.id);
        if (id) clReleaseCommandQueue(id);
        id = q.id;
      }
      return *this;
    }

    queue_t() : id(0) {
    }

    queue_t(const queue_t &q) : id(0) {
      *this = q;
    }

    queue_t(cl_command_queue _id) : id(_id) {
      clRetainCommandQueue(id);
    }

    queue_t(cl_context ctx, cl_device_id dev, cl_command_queue_properties prop = 0) {
      cl_int err;
      id = clCreateCommandQueue(ctx, dev, prop, &err);
      _check(err);
    }

    ~queue_t() {
      if (id) clReleaseCommandQueue(id);
    }

    cl_uint         refcount() const { return info_type::get<cl_uint>(id, CL_QUEUE_REFERENCE_COUNT); }
    cl_context       context() const { return info_type::get<cl_context>(id, CL_QUEUE_CONTEXT);      }
    cl_device_id      device() const { return info_type::get<cl_device_id>(id, CL_QUEUE_DEVICE);     }

    cl_command_queue_properties properties() const {
      return info_type::get<cl_command_queue_properties>(id, CL_QUEUE_PROPERTIES);
    }

    template<typename cmd_t>
    event_t async(const cmd_t &cmd) {
      event_t ev;
      _check(cmd(id, 0, NULL, &(ev.id)));
      return ev;
    }

    template<typename cmd_t, typename iter_t>
    event_t async(const cmd_t &cmd, iter_t wait_beg, iter_t wait_end) {
      ssize_t n_wait = std::distance(wait_beg, wait_end);
      cl_event wait[n_wait];
      for (size_t i =0; wait_beg != wait_end; ++i, ++wait_beg) {
        wait[i] = *wait_beg;
      }

      event_t ev;
      _check(cmd(id, n_wait, wait, &(ev.id)));
      return ev;
    }

    template<typename cmd_t>
    event_t async(const cmd_t &cmd, cl_event event) {
      return async(cmd, &event, &event+1);
    }

    template<typename cmd_t>
    void sync(const cmd_t &cmd) {
      async(cmd).wait();
    }

    template<typename cmd_t, typename iter_t>
    void sync(const cmd_t &cmd, iter_t wait_beg, iter_t wait_end) {
      async(cmd, wait_beg, wait_end).wait();
    }

    template<typename cmd_t>
    void sync(const cmd_t &cmd, cl_event event) {
      async(cmd, event).wait();
    }

    operator cl_command_queue() const { return id; }
  };



  struct opencl_t {
    std::vector<platform_t *> cl_platforms;

    opencl_t() {
    }

    ~opencl_t() {
      for (size_t i = 0; i < cl_platforms.size(); ++i) {
        delete cl_platforms[i];
      }
    }

    void init() {
      cl_uint n_plat;
      _check(clGetPlatformIDs(0, NULL, &n_plat));

      cl_platforms.reserve(n_plat);
      cl_platform_id *plat;
      plat = new cl_platform_id[n_plat];
      _check(clGetPlatformIDs(n_plat, plat, NULL));
      for (cl_uint i = 0; i < n_plat; ++i) cl_platforms.push_back(new platform_t(plat[i]));
      delete [] plat;
    }

    void list(std::ostream &out) const {
      for (size_t i = 0; i < cl_platforms.size(); ++i) {
        const platform_t &plat = *cl_platforms[i];
        out << plat.id << " " << plat.name()
            << " (" << plat.profile() << " ver: " << plat.version() << " vendor: " << plat.vendor() << ")" << std::endl;

        for (size_t j = 0; j < plat.cl_devices.size(); ++j) {
          const dev_t &dev = *plat.cl_devices[j];
          out << "\t" << dev.id << " " << dev.name()
              << " (" << dev.profile() << " ver: " << dev.version() << " vendor: " << dev.vendor() << ")" << std::endl;
        }
      }
    }

    std::vector<dev_t *> gpu_devices() { 
      std::vector<dev_t *> out;
      for (size_t i = 0; i < cl_platforms.size(); ++i) {
        platform_t &plat = *cl_platforms[i];
        for (size_t j = 0; j < plat.cl_devices.size(); ++j) {
          dev_t *dev = plat.cl_devices[j];
          if (dev->is_gpu()) out.push_back(dev);
        }
      }
      return out;
    }

    std::vector<dev_t *> cpu_devices() { 
      std::vector<dev_t *> out;
      for (size_t i = 0; i < cl_platforms.size(); ++i) {
        platform_t &plat = *cl_platforms[i];
        for (size_t j = 0; j < plat.cl_devices.size(); ++j) {
          dev_t *dev = plat.cl_devices[j];
          if (dev->is_cpu()) out.push_back(dev);
        }
      }
      return out;
    }
  };
}
