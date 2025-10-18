cmake_minimum_required(VERSION 3.20)

set( BIFROST_VERSION "2.15.0.0")
set( BIFROST_LOCATION "C:/Program Files/Autodesk/Bifrost/Maya2026/${BIFROST_VERSION}/bifrost" )
set( PRJ_NAME "Snap" )

# Project Name
project(
    Snap 
    VERSION 1.0.0
    LANGUAGES CXX
)

# Bifrost setup and utilities...
include(${BIFROST_LOCATION}/sdk/cmake/setup.cmake)

# Verbose CMake
set( CMAKE_VERBOSE_MAKEFILE TRUE)

# Find Bifrost
find_package(Bifrost REQUIRED SDK)

# Build from sources
add_subdirectory(src)
