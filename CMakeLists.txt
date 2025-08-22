cmake_minimum_required(VERSION 3.16.3)
project(qobi)

set(CMAKE_BUILD_TYPE Debug)
set(CMAKE_BUILD_TYPE RelWithDebInfo)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

file(GLOB SRC CONFIGURE_DEPENDS "src/*")

if(WIN32)
include_directories(../INCLUDE ./include)
link_directories(../LIB)
else()
include_directories(./include)
find_package(glfw CONFIG REQUIRED)
endif()

add_executable(QOBI ${SRC} glad.c)

if(WIN32)
target_link_libraries(QOBI glfw3)
else()
target_link_libraries(QOBI glfw)
endif()

add_executable(TestNomenclature src/tests/nomenclature_test.cpp src/molecule.cpp src/atom.cpp)
add_test(NAME TestNomenclature COMMAND TestNomenclature)
