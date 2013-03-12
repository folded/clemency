include_directories("${clemency_SOURCE_DIR}/include")

if(OPENMP_FOUND)
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()

if(OPENCL_FOUND)
  include_directories(${OPENCL_INCLUDE_DIRS})
endif()

if(OPENCL_FOUND)
  add_executable       (clmc opencl_mc.cpp mc_table.cpp)
  target_link_libraries(clmc ${Boost_PROGRAM_OPTIONS_LIBRARY} ${OPENCL_LIBRARIES})
endif()

