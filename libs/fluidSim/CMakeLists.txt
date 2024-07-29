cmake_minimum_required(VERSION 2.8...3.30)

project(fluidSim 
    VERSION 0.2.0
    LANGUAGES CXX 
    DESCRIPTION "Library of common helper functions and classes for OpenGL projects. Derived mainly from 3D Graphics Rendering Cookbook"
    )

# Debugging cmake
include(CMakePrintHelpers)

message(STATUS "Building" ${PROJECT_NAME})

#INCLUDE ANY OF OUR CMAKE HELPERS
include(${CMAKE_SOURCE_DIR}/cmake/setup-helpers.cmake)

find_package(glm CONFIG REQUIRED)
find_package(Armadillo CONFIG REQUIRED)

file(GLOB_RECURSE SRC_FILES LIST_DIRECTORIES false RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} *.c??)
file(GLOB_RECURSE HEADER_FILES LIST_DIRECTORIES false RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} *.h??)

cmake_print_variables(CMAKE_CURRENT_SOURCE_DIR)

message(STATUS "Checking file globs...")

cmake_print_variables(SRC_FILES)
cmake_print_variables(HEADER_FILES)

SETUP_GROUPS("${HEADER_FILES}")
SETUP_GROUPS("${SRC_FILES}")

add_library(fluidSim)

# export all our symbols for now
set(CMAKE_WINDOWS_EXPORT_ALL_SYMBOLS ON)
# TODO: refine our exports as in Chapter 3 of CMake Best Practices
# set_property(TARGET fluidSim PROPERTY CXX_VISIBILITY_PRESET "hidden")
# set_property(TARGET fluidSim PROPERTY VISIBILITY_INLINES_HIDDEN TRUE)
# include(GenerateExportHeader)
# generate_export_header(fluidSim EXPORT_FILE_NAME export/fluidSim/fluidSim.hpp)
# target_include_directories(jcfluidSim PUBLIC "${CMAKE_CURRENT_BINARY_DIR} /export")

target_sources(fluidSim PRIVATE ${SRC_FILES} ${HEADER_FILES})

target_compile_features(fluidSim PUBLIC cxx_std_20)

target_include_directories(fluidSim 
  PRIVATE 
  ${CMAKE_CURRENT_SOURCE_DIR}/src 
  ${CMAKE_CURRENT_SOURCE_DIR}/src/fluidSim
  ${CMAKE_CURRENT_SOURCE_DIR}/include/fluidSim
  PUBLIC
  ${CMAKE_CURRENT_SOURCE_DIR}/include
#	${CMAKE_CURRENT_SOURCE_DIR}/include/fluidSim
#	${CMAKE_CURRENT_SOURCE_DIR}/include/fluidSim/scene
	)

target_link_libraries(fluidSim 
    PRIVATE
    glm::glm
    armadillo
)

# Create virtual folders for Visual Studio
if(MSVC_IDE)
	# Macro to preserve source files hierarchy in the IDE
	macro(GroupSources curdir)
		file(GLOB children RELATIVE ${PROJECT_SOURCE_DIR}/${curdir} ${PROJECT_SOURCE_DIR}/${curdir}/*)
		
		foreach(child ${children})
			if(IS_DIRECTORY ${PROJECT_SOURCE_DIR}/${curdir}/${child})
				GroupSources(${curdir}/${child})
			else() 
				string(REPLACE "/" "\\" groupname ${curdir})
				string(REPLACE "src" "Sources" groupname ${groupname})
				source_group(${groupname} FILES ${PROJECT_SOURCE_DIR}/${curdir}/${child})
			endif()
		endforeach()
	endmacro()
	
	# Run Macro
	GroupSources(src)
	GroupSources(include/fluidSim)
endif()