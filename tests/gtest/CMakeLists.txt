# Second level CMakeLists.txt: register individual tests

include(GoogleTest)

# Add tests: specify add_executable and target_link_libraries

add_executable(dmultinom_test
  dmultinom_upd-unit.cpp
)

target_include_directories(dmultinom_test
  PUBLIC
    ${CMAKE_SOURCE_DIR}/../
)

target_link_libraries(dmultinom_test
  gtest_main
)

gtest_discover_tests(dmultinom_test)
