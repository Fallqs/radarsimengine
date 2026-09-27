# Phase-0 smoke + Phase-1 golden tests. Registered with CTest so the upstream
# build scripts (--test=on) exercise them via `ctest`.
add_executable(rsim_smoke test_smoke.cpp)
target_link_libraries(rsim_smoke PRIVATE radarsimcpp)
add_test(NAME rsim_smoke COMMAND rsim_smoke)

add_executable(rsim_point test_point.cpp)
target_link_libraries(rsim_point PRIVATE radarsimcpp)
add_test(NAME rsim_point COMMAND rsim_point)

add_executable(rsim_geometry test_geometry.cpp)
target_link_libraries(rsim_geometry PRIVATE radarsimcpp)
add_test(NAME rsim_geometry COMMAND rsim_geometry)

# keep assert() live in Release builds (the smoke test is assert-based)
if(CMAKE_CXX_COMPILER_ID MATCHES "MSVC")
    set(RSIM_TEST_UNDEBUG /UNDEBUG)
else()
    set(RSIM_TEST_UNDEBUG -UNDEBUG)
endif()
target_compile_options(rsim_smoke PRIVATE ${RSIM_TEST_UNDEBUG})
target_compile_options(rsim_point PRIVATE ${RSIM_TEST_UNDEBUG})
target_compile_options(rsim_geometry PRIVATE ${RSIM_TEST_UNDEBUG})

# the test executables need radarsimcpp.dll next to them on Windows
if(WIN32)
    set_target_properties(rsim_smoke rsim_point rsim_geometry PROPERTIES
        RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})
endif()
