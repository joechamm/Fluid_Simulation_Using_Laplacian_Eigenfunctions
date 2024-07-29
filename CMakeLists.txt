cmake_minimum_required(VERSION 2.8...3.30)

project(Fluid_Simulation_Using_Laplacian_Eigenfunctions)

project(Fluid_Simulation_Using_Laplacian_Eigenfunctions
    VERSION 0.2.1
    DESCRIPTION "This repository contains my implementation of the paper 'Fluid Simulation Using Laplacian Eigenfunctions' by Tyler De Witt, Adam Bargteil, and Ladislav Kavan."
    LANGUAGES CXX C)

# CMake Debugging
include(CMakePrintHelpers)

message(STATUS "BBuilding ${PROJECT_NAME}")

# use folders
set_property(GLOBAL PROPERTY USE_FOLDERS ON)

# Configure assets header file
configure_file(src/helpers/RootDir.h.in src/helpers/RootDir.h)
include_directories(${CMAKE_BINARY_DIR}/src)

if(WIN32)
    add_definitions(-D_CRT_SECURE_NO_WARNINGS)
endif()

find_package(glm CONFIG REQUIRED)
find_package(glfw3 CONFIG REQUIRED)
find_package(Armadillo CONFIG REQUIRED)
find_package(glad CONFIG REQUIRED)
find_package(gli CONFIG REQUIRED)
find_package(Stb REQUIRED)

# include_directories(.) ...
include_directories(.)
include_directories(./libs)

# add shared utils code 
add_subdirectory(libs/commonFramework)
add_subdirectory(libs/glFramework)
add_subdirectory(libs/fluidSim)

add_subdirectory(src/apps/demo)