cmake_minimum_required(VERSION 3.14)

project(DEMO CXX C)

find_package(glad CONFIG REQUIRED)
find_package(armadillo CONFIG REQUIRED)
find_package(glfw3 CONFIG REQUIRED)
find_package(glm CONFIG REQUIRED)

include_directories(${CMAKE_SOURCE_DIR}/libs/commonFramework/include)
include_directories(${CMAKE_SOURCE_DIR}/libs/glFramework/include)

add_executable(demo)

set_property(TARGET demo PROPERTY FOLDER "apps")

target_compile_features(demo PRIVATE cxx_std_20)

target_sources(demo PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/src/main.cpp)

target_link_libraries(demo PRIVATE
    glad
    glfw3
    armadillo
)
