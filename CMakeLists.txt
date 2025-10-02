cmake_minimum_required(VERSION 3.14)
project(xtc VERSION 1.0.0 LANGUAGES CXX)

# Set C++ standard
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Options
option(XTC_BUILD_EXAMPLES "Build example programs" OFF)
option(XTC_ENABLE_OPENMP "Enable OpenMP support" ON)

# Find OpenMP if enabled
if(XTC_ENABLE_OPENMP)
    find_package(OpenMP)
    if(OpenMP_CXX_FOUND)
        set(XTC_USE_OPENMP TRUE)
    else()
        message(STATUS "OpenMP not found, building without parallel support")
        set(XTC_USE_OPENMP FALSE)
    endif()
else()
    set(XTC_USE_OPENMP FALSE)
endif()

# Create library
add_library(xtc
    libxtc/xtc_unpack.cpp
    libxtc/xtc_unpack.h
    libxtc/xtc_reader.cpp
    libxtc/xtc_reader.h
    libxtc/xtc_writer.cpp
    libxtc/xtc_writer.h
)

# Set library properties
set_target_properties(xtc PROPERTIES
    VERSION ${PROJECT_VERSION}
    SOVERSION ${PROJECT_VERSION_MAJOR}
    PUBLIC_HEADER "libxtc/xtc_reader.h;libxtc/xtc_unpack.h"
)

# Include directories
target_include_directories(xtc
    PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
        $<INSTALL_INTERFACE:include>
)

# Link OpenMP if available
if(XTC_USE_OPENMP)
    target_link_libraries(xtc PUBLIC OpenMP::OpenMP_CXX)
    # Don't define _OPENMP - OpenMP already defines it
endif()

# Compiler-specific flags
if(MSVC)
    target_compile_options(xtc PRIVATE /W4)
else()
    target_compile_options(xtc PRIVATE -Wall -Wextra -O3)
endif()

# Add alias for namespaced usage
add_library(xtc::xtc ALIAS xtc)

# Installation
include(GNUInstallDirs)

install(TARGETS xtc
    EXPORT xtcTargets
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
    PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/libxtc
)

# Install export targets
install(EXPORT xtcTargets
    FILE xtcTargets.cmake
    NAMESPACE xtc::
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/xtc
)

# Create and install config files
include(CMakePackageConfigHelpers)

configure_package_config_file(
    ${CMAKE_CURRENT_SOURCE_DIR}/cmake/xtcConfig.cmake.in
    ${CMAKE_CURRENT_BINARY_DIR}/xtcConfig.cmake
    INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/xtc
)

write_basic_package_version_file(
    ${CMAKE_CURRENT_BINARY_DIR}/xtcConfigVersion.cmake
    VERSION ${PROJECT_VERSION}
    COMPATIBILITY SameMajorVersion
)

install(FILES
    ${CMAKE_CURRENT_BINARY_DIR}/xtcConfig.cmake
    ${CMAKE_CURRENT_BINARY_DIR}/xtcConfigVersion.cmake
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/xtc
)

# Export targets for build tree usage
export(EXPORT xtcTargets
    FILE ${CMAKE_CURRENT_BINARY_DIR}/xtcTargets.cmake
    NAMESPACE xtc::
)

# Build examples if requested
if(XTC_BUILD_EXAMPLES)
    add_subdirectory(examples)
endif()

# Print configuration summary
message(STATUS "")
message(STATUS "XTC Configuration Summary:")
message(STATUS "  Version:              ${PROJECT_VERSION}")
message(STATUS "  C++ Standard:         ${CMAKE_CXX_STANDARD}")
message(STATUS "  OpenMP Support:       ${XTC_USE_OPENMP}")
message(STATUS "  Build Examples:       ${XTC_BUILD_EXAMPLES}")
message(STATUS "  Install Prefix:       ${CMAKE_INSTALL_PREFIX}")
message(STATUS "")
