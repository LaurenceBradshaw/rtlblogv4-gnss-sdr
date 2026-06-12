# Sources are grouped into folders (common/ io/ signals/ acquisition/ tracking/ engine/ pvt/),
# globbed recursively. Headers use flat includes ("channel.h", not "engine/channel.h"), so every
# source folder is added to the include path below - filenames are unique, so there is no ambiguity.
file(GLOB_RECURSE sdr_core_sources CONFIGURE_DEPENDS "*.cpp")

# The gui/ sources need Qt and are compiled only when SDR_GUI is on (see below); keep them out of
# the core/test sources so the headless build and the test binary never depend on Qt.
file(GLOB_RECURSE sdr_gui_sources CONFIGURE_DEPENDS "gui/*.cpp")
if(sdr_gui_sources)
    list(REMOVE_ITEM sdr_core_sources ${sdr_gui_sources})
endif()

set(sdr_include_dirs
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_CURRENT_SOURCE_DIR}/common
    ${CMAKE_CURRENT_SOURCE_DIR}/io
    ${CMAKE_CURRENT_SOURCE_DIR}/signals
    ${CMAKE_CURRENT_SOURCE_DIR}/signals/gps
    ${CMAKE_CURRENT_SOURCE_DIR}/signals/galileo
    ${CMAKE_CURRENT_SOURCE_DIR}/signals/beidou
    ${CMAKE_CURRENT_SOURCE_DIR}/acquisition
    ${CMAKE_CURRENT_SOURCE_DIR}/tracking
    ${CMAKE_CURRENT_SOURCE_DIR}/engine
    ${CMAKE_CURRENT_SOURCE_DIR}/pvt
)

find_package(cxxopts CONFIG REQUIRED)
find_package(FFTW3f CONFIG REQUIRED)
find_package(fmt CONFIG REQUIRED)
find_package(Eigen3 CONFIG REQUIRED)

# Live RTL-SDR front-end (Rtlsdr_device). Requires librtlsdr (system pkg-config package).
find_package(PkgConfig REQUIRED)
pkg_check_modules(RTLSDR REQUIRED librtlsdr)

add_executable(sdr ${sdr_core_sources} )
target_include_directories(sdr PRIVATE ${sdr_include_dirs} ${RTLSDR_INCLUDE_DIRS})
target_link_libraries(sdr PRIVATE
        cxxopts::cxxopts
        fmt::fmt
        Eigen3::Eigen
        FFTW3::fftw3f
        ${RTLSDR_LIBRARIES}
)

# Optional Qt GUI (run with: sdr --gui). Uses the system Qt6 (find_package), not vcpkg, so the
# default build stays light. Only the main `sdr` executable gets Qt + the gui/ sources + AUTOMOC;
# the test binary never does.
option(SDR_GUI "Build the Qt GUI (run with: sdr --gui)" OFF)
if(SDR_GUI)
    find_package(Qt6 COMPONENTS Widgets Network REQUIRED)
    target_sources(sdr PRIVATE ${sdr_gui_sources})
    target_include_directories(sdr PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/gui)
    target_link_libraries(sdr PRIVATE Qt6::Widgets Qt6::Network)
    target_compile_definitions(sdr PRIVATE SDR_GUI)
    set_target_properties(sdr PROPERTIES AUTOMOC ON)
endif()

if(ENABLE_UNIT_TESTS)
    add_executable(test_sdr ${sdr_core_sources})
    target_compile_definitions(test_sdr PRIVATE ENABLE_UNIT_TESTS)
    target_include_directories(test_sdr PRIVATE ${sdr_include_dirs} ${RTLSDR_INCLUDE_DIRS})
    target_link_libraries(test_sdr PRIVATE cxxopts::cxxopts fmt::fmt Eigen3::Eigen FFTW3::fftw3f ${RTLSDR_LIBRARIES} Catch2::Catch2WithMain)

    add_test(NAME test_sdr COMMAND test_sdr)
endif()
