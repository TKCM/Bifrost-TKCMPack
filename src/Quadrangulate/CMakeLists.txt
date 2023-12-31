cmake_minimum_required(VERSION 3.17)
set( BIFROST_LOCATION "C:/Program Files/Autodesk/Bifrost/Maya2024/2.6.0.0-dev-202211041459-29b8754/bifrost")
set( PRJ_NAME "Quadrangulate" )

# Project Name
project(
    ${PRJ_NAME} 
    VERSION 1.0.0
    LANGUAGES CXX
)

include_directories(${BIFROST_LOCATION}/sdk/include)
include_directories(${BIFROST_LOCATION}/thirdparty/tbb_include)
include_directories("C:/D/GitHub/BifrostExts/eigen")

# Bifrost setup and utilities...
include(${BIFROST_LOCATION}/sdk/cmake/setup.cmake)

# Verbose CMake
set( CMAKE_VERBOSE_MAKEFILE TRUE)

# Find Bifrost
find_package(Bifrost REQUIRED SDK)

# Build from sources
add_subdirectory(src)
