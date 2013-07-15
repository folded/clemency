include_directories("${clemency_SOURCE_DIR}/include")

if(OPENMP_FOUND)
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()

if(OPENCL_FOUND)
  include_directories(${OPENCL_INCLUDE_DIRS})
endif()

add_executable       (mesh_info mesh_info.cpp)
target_link_libraries(mesh_info ${Boost_PROGRAM_OPTIONS_LIBRARY} clemency gloop_model)

add_executable       (mesh_clean mesh_clean.cpp)
target_link_libraries(mesh_clean ${Boost_PROGRAM_OPTIONS_LIBRARY} clemency gloop_model)

add_executable       (flip flip.cpp)
target_link_libraries(flip ${Boost_PROGRAM_OPTIONS_LIBRARY} clemency gloop_model)

if(OPENCL_FOUND)
  add_executable       (clmc opencl_mc.cpp mc_table.cpp)
  target_link_libraries(clmc ${Boost_PROGRAM_OPTIONS_LIBRARY} ${OPENCL_LIBRARIES} clemency gloop_model)
endif()

