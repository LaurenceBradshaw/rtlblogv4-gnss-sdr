find_package(fmt REQUIRED)

set(COMMON_SOURCE_DIR ${CMAKE_CURRENT_LIST_DIR})

file(GLOB_RECURSE LIB_UTILS_HEADERS CONFIGURE_DEPENDS
    "${COMMON_SOURCE_DIR}/*.h"
)

file(GLOB_RECURSE LIB_UTILS_SOURCES CONFIGURE_DEPENDS
    "${COMMON_SOURCE_DIR}/*.cpp"
)

add_library(lib_utils STATIC
    ${LIB_UTILS_HEADERS}
    ${LIB_UTILS_SOURCES}
)

target_link_libraries(lib_utils
    PUBLIC
        fmt::fmt
)

target_include_directories(lib_utils
    PUBLIC
        ${COMMON_SOURCE_DIR}
)