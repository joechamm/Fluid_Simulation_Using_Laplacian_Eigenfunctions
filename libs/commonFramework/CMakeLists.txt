cmake_minimum_required(VERSION 2.8...3.30)

project(commonFramework 
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

file(GLOB_RECURSE SRC_FILES LIST_DIRECTORIES false RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} *.c??)
file(GLOB_RECURSE HEADER_FILES LIST_DIRECTORIES false RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} *.h??)

cmake_print_variables(CMAKE_CURRENT_SOURCE_DIR)

message(STATUS "Checking file globs...")

cmake_print_variables(SRC_FILES)
cmake_print_variables(HEADER_FILES)

SETUP_GROUPS("${HEADER_FILES}")
SETUP_GROUPS("${SRC_FILES}")

add_library(commonFramework)

# export all our symbols for now
set(CMAKE_WINDOWS_EXPORT_ALL_SYMBOLS ON)
# TODO: refine our exports as in Chapter 3 of CMake Best Practices
# set_property(TARGET commonFramework PROPERTY CXX_VISIBILITY_PRESET "hidden")
# set_property(TARGET commonFramework PROPERTY VISIBILITY_INLINES_HIDDEN TRUE)
# include(GenerateExportHeader)
# generate_export_header(commonFramework EXPORT_FILE_NAME export/commonFramework/commonFramework.hpp)
# target_include_directories(jcCommonFramework PUBLIC "${CMAKE_CURRENT_BINARY_DIR} /export")

target_sources(commonFramework PRIVATE ${SRC_FILES} ${HEADER_FILES})

target_compile_features(commonFramework PUBLIC cxx_std_20)

target_include_directories(commonFramework 
  PRIVATE 
  ${CMAKE_CURRENT_SOURCE_DIR}/src 
  ${CMAKE_CURRENT_SOURCE_DIR}/src/commonFramework
  ${CMAKE_CURRENT_SOURCE_DIR}/include/jcCommonFramework
  PUBLIC
  ${CMAKE_CURRENT_SOURCE_DIR}/include
#	${CMAKE_CURRENT_SOURCE_DIR}/include/commonFramework
#	${CMAKE_CURRENT_SOURCE_DIR}/include/commonFramework/scene
	)

target_link_libraries(commonFramework 
    PRIVATE
    glm::glm
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
	GroupSources(include/commonFramework)
endif()