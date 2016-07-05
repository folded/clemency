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

#include <gtest/gtest.h>
#include <clemency/cl.hpp>

TEST(OpenCL, OpenCL) {
  cl::opencl_t opencl;

  opencl.init();

  cl::dev_t* gpu_dev = opencl.gpu_devices()[0];

  ASSERT_TRUE(gpu_dev != NULL);
  ASSERT_TRUE(gpu_dev->is_gpu());
  ASSERT_FALSE(gpu_dev->is_cpu());

  cl::ctx_t ctx(CL_DEVICE_TYPE_GPU);

  ASSERT_TRUE(ctx.id != 0);

  ASSERT_EQ(ctx.device(0), gpu_dev->id);
  ASSERT_EQ(ctx.refcount(), 1);

  cl::queue_t queue(ctx, *gpu_dev);

  cl::program_t prog = ctx.create_program_from_source(
      "__kernel void saxpy(const float alpha,\n"
      "                    __global const float* X,\n"
      "                    __global float* Y)\n"
      "{\n"
      "    Y[get_global_id(0)] += alpha * X[get_global_id(0)];\n"
      "}\n");

  ASSERT_EQ(prog.build_status(*gpu_dev), CL_BUILD_NONE);
  ASSERT_TRUE(prog.build());
  ASSERT_EQ(prog.build_status(*gpu_dev), CL_BUILD_SUCCESS);

  const size_t N = 256;
  float cpuX[N], cpuY[N];

  for (size_t i = 0; i < N; ++i) {
    cpuX[i] = float(i);
    cpuY[i] = 1.0;
  }

  cl_mem memX = ctx.create_buffer(CL_MEM_READ_ONLY | CL_MEM_USE_HOST_PTR,
                                  N * sizeof(float), cpuX);
  cl_mem memY = ctx.create_buffer(CL_MEM_READ_WRITE | CL_MEM_USE_HOST_PTR,
                                  N * sizeof(float), cpuY);

  ASSERT_NE(memX, (cl_mem)NULL);
  ASSERT_NE(memY, (cl_mem)NULL);

  {
    cl::event_t e1 = queue.async(
        cl::cmd_write_buffer(memX, CL_FALSE, 0, N * sizeof(float), cpuX));
    cl::event_t e2 = queue.async(
        cl::cmd_write_buffer(memY, CL_FALSE, 0, N * sizeof(float), cpuY), e1);
    e2.wait();
  }

  cl::kernel_t kern = prog.kernel("saxpy");
  ASSERT_EQ(kern.name(), "saxpy");

  kern.arg(0, 1.5f);
  kern.arg(1, memX);
  kern.arg(2, memY);

  queue.sync(cl::cmd_task(kern).wrk_sz(N));

  {
    cl::event_t e1 = queue.async(
        cl::cmd_read_buffer(memY, CL_FALSE, 0, N * sizeof(float), cpuY));
    e1.wait();
  }

  for (size_t i = 0; i < N; ++i) {
    ASSERT_EQ(cpuY[i], 1.0f + 1.5f * cpuX[i]);
  }

  clReleaseMemObject(memX);
  clReleaseMemObject(memY);
}
