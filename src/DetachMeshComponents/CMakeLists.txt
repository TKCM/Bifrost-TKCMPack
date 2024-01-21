cmake_minimum_required(VERSION 3.17)

set( BIFROST_VERSION "2.8.0.0")
set( BIFROST_LOCATION "C:/Program Files/Autodesk/Bifrost/Maya2024/${BIFROST_VERSION}/bifrost" )
set( PRJ_NAME "DetachMeshComponents" )

# Project Name
project(
    ${PRJ_NAME} 
    VERSION 1.0.0
    LANGUAGES CXX
)
include_directories( ${BIFROST_LOCATION}/sdk/include )

# Bifrost setup and utilities...
include( ${BIFROST_LOCATION}/sdk/cmake/setup.cmake )

# Verbose CMake
set( CMAKE_VERBOSE_MAKEFILE TRUE )

# Find Bifrost
find_package( Bifrost REQUIRED SDK )

# Build from sources
add_subdirectory( src )

# openMP
find_package(OpenMP REQUIRED)
if(OpenMP_FOUND)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()