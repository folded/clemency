clemency
========

An OpenCL marching cubes implementation.

Copyright (C) 2013, Tobias Sargeant <tobias.sargeant@gmail.com>

All rights reserved.

Usage
=====

	clmc [options] kernel
	
	Allowed options:
	
	  -h [ --help ]                       Print help and exit
	  -v [ --version ]                    Print version string
	  -c [ --cpu ]                        Execute on CPU
	  -l [ --list ]                       List available OpenCL devices
	  -C [ --cl-opts ] arg                Set OpenCL compiler option
	  -g [ --grid ] arg (=256)            Size of marching cubes grid
	  --gx arg (=0)                       X dimension of mc grid
	  --gy arg (=0)                       Y dimension of mc grid
	  --gz arg (=0)                       Z dimension of mc grid
	  -a [ --aabb ] arg (=-1,-1,-1:1,1,1) aabb of marching cubes grid
	  --block0 arg (=256)                 Block size for grid evaluation
	  --block1 arg (=1024)                Block size for intersection refinement
	  -r [ --refine ] arg (=10)           Number of refinement steps
	  -o [ --output ] arg (=-)            Ouput file (- for stdout)


Output is produced in stanford .ply format.

New Kernels
===========

New kernels should define a function 'inside' with prototype:

	bool inside(float3 pos);

that returns true if pos is inside the marching cubes isosurface.

This should be followed by:

	#include <core.cl>

to pull in the two OpenCL kernels (that make use of inside()) required by clmc.

Executing a new kernel might require providing an include path
(using an argument like -C "-I /path/to/opencl") to enable the
OpenCL compiler to find core.cl.

Similarly, you can use -C "-D PARAMETER=value" to set kernel parameters,
if needed.

Sample kernels are provided in the opencl directory.
