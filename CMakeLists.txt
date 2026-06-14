cmake_minimum_required(VERSION 3.24)

project(rtlblogv4_gnss_sdr VERSION 0.1.0 LANGUAGES CXX)

option(ENABLE_UNIT_TESTS "Build unit test binaries" ON)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# Debug builds use -Og instead of the default -O0. -O0 is too slow for high-rate captures (a 25 MHz file
# decimated to ~2 MHz runs ~0.5x real time at -O0 but ~real time at -Og), and -Og keeps full debug info
# and steppability (unlike -O3). To force plain -O0, edit this line (or override CMAKE_CXX_FLAGS_DEBUG
# in a fresh build dir before this runs).
set(CMAKE_CXX_FLAGS_DEBUG "-Og -g")

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${PROJECT_SOURCE_DIR}/bin")
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${PROJECT_SOURCE_DIR}/bin")
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${PROJECT_SOURCE_DIR}/bin")

if(ENABLE_UNIT_TESTS)
  include(CTest)
  find_package(Catch2 CONFIG REQUIRED)
endif()

add_subdirectory(sdr)