cmake_minimum_required(VERSION 2.6)

project(gloop)

set(gloop_VERSION_MAJOR 1)
set(gloop_VERSION_MINOR 0)
set(gloop_VERSION_PATCH 0)

set(GLOOP_VERSION ${gloop_VERSION_MAJOR}.${gloop_VERSION_MINOR}.${gloop_VERSION_PATCH})

include_directories("${gloop_SOURCE_DIR}/include")

find_package(OpenGL)
find_package(GLEW)

if (GLEW_FOUND)
  if (OPENGL_FOUND)
    include_directories(${OPENGL_INCLUDE_DIR})
    include_directories(${GLEW_INCLUDE_DIR})

    add_library(gloop STATIC src/shader.cpp src/surface.cpp src/texparam.cpp src/vbo.cpp src/fbo.cpp)
    target_link_libraries(gloop ${OPENGL_LIBRARIES} ${GLEW_LIBRARY})
  endif()
endif()

add_library(gloop_math STATIC src/matrix.cpp src/quaternion.cpp)

add_library(gloop_image STATIC src/radiance.cpp)

add_library(gloop_model STATIC  src/model/stream.cpp src/model/ply_format.cpp src/model/obj_format.cpp src/model/vtk_format.cpp)

set(GLOOP_MODEL_LIBRARY gloop_model)
