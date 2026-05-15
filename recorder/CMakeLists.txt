cmake_minimum_required(VERSION 3.20)
project(Rtlblogv4_recorder C CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

function(SET_DEBUG_POSTFIX target)
    set_target_properties(${target} PROPERTIES DEBUG_POSTFIX "_d")
endfunction()

find_package(cxxopts CONFIG REQUIRED)
find_package(Eigen3 CONFIG REQUIRED)
find_package(fmt CONFIG REQUIRED)

add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/../libs/lib_utils ${CMAKE_BINARY_DIR}/lib_utils)

file(GLOB SRC_FILES CONFIGURE_DEPENDS
    "*.cpp"
)

file(GLOB HDR_FILES CONFIGURE_DEPENDS
    "*.h"
)

add_executable( 
	Rtlblogv4_recorder
	${SRC_FILES}
	${HDR_FILES}
)

target_include_directories(Rtlblogv4_recorder 
	PRIVATE 
		${CMAKE_CURRENT_BINARY_DIR}
)

SET_DEBUG_POSTFIX( Rtlblogv4_recorder )

target_link_libraries(Rtlblogv4_recorder
	PRIVATE
        lib_utils
        Threads::Threads
        Eigen3::Eigen
        cxxopts::cxxopts
        fmt::fmt
)