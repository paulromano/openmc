set(TEST_NAMES
  test_distribution
  test_file_utils
  test_tally
  test_interpolate
  test_math
  test_mcpl_stat_sum
  test_mesh
  test_photon
  test_ray
  test_region
  test_tensor
  test_geometry
  test_recoil
  # Add additional unit test files here
)

foreach(test ${TEST_NAMES})
  add_executable(${test} ${test}.cpp)
  target_link_libraries(${test} Catch2::Catch2WithMain libopenmc)
  # Some tests replay a checked-in fixture. The binary runs from the build
  # tree, so it is told where the source-tree data lives rather than guessing.
  target_compile_definitions(${test} PRIVATE
    OPENMC_UNIT_TEST_DATA_DIR="${CMAKE_CURRENT_SOURCE_DIR}/data")
  add_test(NAME ${test} COMMAND ${test} WORKING_DIRECTORY ${UNIT_TEST_BIN_OUTPUT_DIR})
endforeach()
