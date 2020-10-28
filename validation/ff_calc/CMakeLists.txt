cmake_minimum_required(VERSION 3.15...3.18)

project(ff_calc VERSION 1.1
                DESCRIPTION "Form factor calcuator"
                LANGUAGES CXX)

# Require out-of-source builds
file(TO_CMAKE_PATH "${PROJECT_BINARY_DIR}/CMakeLists.txt" LOC_PATH)
if(EXISTS "${LOC_PATH}")
    message(FATAL_ERROR "You cannot build in a source directory (or any directory with a CMakeLists.txt file). Please make a build subdirectory. Feel free to remove CMakeCache.txt and CMakeFiles.")
endif()

# Set required C++ standard version
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED on)
set(CMAKE_CXX_EXTENSIONS OFF)

# Find required packages
find_package(ROOT)

# Targets
add_library(ff_dstaunu SHARED src/ff_dstaunu.cpp inc/ff_dstaunu.hpp)

target_include_directories(
    ff_dstaunu
  PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/inc>
    $<INSTALL_INTERFACE:inc>
)

# Define install rules
include(GNUInstallDirs)

install(TARGETS ff_dstaunu
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
