include(FetchContent)

# Download and build GoogleTest
FetchContent_Declare(
  googletest
  URL https://github.com/google/googletest/archive/release-1.12.1.zip
  DOWNLOAD_EXTRACT_TIMESTAMP TRUE
)


# Needed for GoogleTest on some platforms (especially Windows)
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)

FetchContent_MakeAvailable(googletest)

# Define the test executable
add_executable(test_autocalibration
  test_autocalibration.cpp
  ../libmv.cpp
)

# Include directories for headers
target_include_directories(test_autocalibration PRIVATE
  ../
  ${CMAKE_CURRENT_SOURCE_DIR}
  /Data/gautt/bundle_adjustment/robust_povar/external/eigen  # <-- Replace with your actual Eigen path
)

# Link with GoogleTest and pthread
target_link_libraries(test_autocalibration
  gtest_main
  pthread
)

# Register tests with CTest
include(GoogleTest)
gtest_discover_tests(test_autocalibration)
