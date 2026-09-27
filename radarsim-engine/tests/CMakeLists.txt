# Phase-0 smoke + Phase-1 golden tests. Registered with CTest so the upstream
# build scripts (--test=on) exercise them via `ctest`.
add_executable(rsim_smoke test_smoke.cpp)
target_link_libraries(rsim_smoke PRIVATE radarsimcpp)
add_test(NAME rsim_smoke COMMAND rsim_smoke)

add_executable(rsim_point test_point.cpp)
target_link_libraries(rsim_point PRIVATE radarsimcpp)
add_test(NAME rsim_point COMMAND rsim_point)

# the test executables need radarsimcpp.dll next to them on Windows
if(WIN32)
    set_target_properties(rsim_smoke rsim_point PROPERTIES
        RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})
endif()
