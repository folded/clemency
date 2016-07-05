include_directories("${clemency_SOURCE_DIR}/include")

if(OPENCL_FOUND)
  include_directories(${OPENCL_INCLUDE_DIRS})
endif()

add_executable       (mesh_info mesh_info.cpp)
target_link_libraries(mesh_info Boost::program_options clemency gloop_model)

add_executable       (mesh_clean mesh_clean.cpp)
target_link_libraries(mesh_clean Boost::program_options clemency gloop_model)

add_executable       (flip flip.cpp)
target_link_libraries(flip Boost::program_options clemency gloop_model)

add_executable       (isosurface isosurface.cpp)
target_link_libraries(isosurface Boost::program_options clemency gloop_model)

if(OPENCL_FOUND)
  add_executable       (clmc opencl_mc.cpp mc_table.cpp)
  target_link_libraries(clmc Boost::program_options ${OPENCL_LIBRARIES} clemency gloop_model)
endif()

if(OPENCL_FOUND)
  add_executable       (clmc-dif opencl_dif.cpp mc_table.cpp)
  target_link_libraries(clmc-dif Boost::program_options ${OPENCL_LIBRARIES} clemency gloop_model)
endif()
