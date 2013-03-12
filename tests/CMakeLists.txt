
if(BUILD_TESTING)
  include(CTest)

  enable_testing()

  find_package(GTest REQUIRED)

  include_directories(${GTEST_INCLUDE_DIRS})

  if(OPENCL_FOUND)
    include_directories(${OPENCL_INCLUDE_DIRS})
  endif()

  include_directories("${clemency_SOURCE_DIR}/include")

  add_executable(test_m3inv test_m3inv.cpp)
  target_link_libraries(test_m3inv ${GTEST_BOTH_LIBRARIES})

  add_executable(test_bezier test_bezier.cpp)
  target_link_libraries(test_bezier ${GTEST_BOTH_LIBRARIES})

  gtest_add_tests(test_m3inv "" test_m3inv.cpp)
  gtest_add_tests(test_bezier "" test_bezier.cpp)

  if(OPENCL_FOUND)
    add_executable       (test_opencl test_opencl.cpp)
    target_link_libraries(test_opencl ${OPENCL_LIBRARIES} ${GTEST_BOTH_LIBRARIES})
  endif()

endif(BUILD_TESTING)
