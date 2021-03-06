# CMake configuration file for Benchmark subdirectory.
# Top level docs for 3.1.3 at: https://cmake.org/cmake/help/v3.1/
# Commands herein described at: https://cmake.org/cmake/help/v3.1/manual/cmake-commands.7.html

# Hides options.
mark_as_advanced(FORCE LIBRT)
mark_as_advanced(FORCE BENCHMARK_BUILD_32_BITS)
mark_as_advanced(FORCE BENCHMARK_ENABLE_EXCEPTIONS)
mark_as_advanced(FORCE BENCHMARK_ENABLE_INSTALL)
mark_as_advanced(FORCE BENCHMARK_ENABLE_TESTING)
mark_as_advanced(FORCE BENCHMARK_ENABLE_LTO)
mark_as_advanced(FORCE BENCHMARK_USE_LIBCXX)

set(BENCHMARK_ENABLE_TESTING OFF CACHE BOOL "Enable testing of the benchmark library.")

# Add subdirectory to build.
# For cmd details, see: https://cmake.org/cmake/help/v3.1/command/add_subdirectory.html
# Adds benchmark here...
add_subdirectory(benchmark)

# Add include directories to build...
# For cmd details, see: https://cmake.org/cmake/help/v3.1/command/include_directories.html
include_directories(${PlayRho_SOURCE_DIR})

file(GLOB Benchmark_SRCS *.cpp)

# Add an executable to the project using specified source files.
# See details at: https://cmake.org/cmake/help/v3.1/command/add_executable.html
add_executable(Benchmark ${Benchmark_SRCS})

# Link a target to given libraries.
# See details at: https://cmake.org/cmake/help/v3.1/command/target_link_libraries.html
target_link_libraries(Benchmark PlayRho benchmark)
