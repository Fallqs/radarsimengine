# Phase-0 smoke tests. Registered with CTest so the upstream build scripts
# (--test=on) exercise them via `ctest`.
add_executable(rsim_smoke test_smoke.cpp)
target_link_libraries(rsim_smoke PRIVATE radarsimcpp)

add_test(NAME rsim_smoke COMMAND rsim_smoke)
