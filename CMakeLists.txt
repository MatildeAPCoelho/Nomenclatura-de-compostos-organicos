cmake_minimum_required(VERSION 3.16.3)
project(qobi)

set(CMAKE_BUILD_TYPE Debug)
#set(CMAKE_BUILD_TYPE RelWithDebInfo)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

file(GLOB SRC CONFIGURE_DEPENDS "src/*")

find_package(glfw3 CONFIG REQUIRED)

include_directories(./include)

add_executable(QOBI ${SRC} glad.c)
target_link_libraries(QOBI glfw)

add_executable(TestNomenclature src/tests/nomenclature_test.cpp src/molecule.cpp src/atom.cpp)
add_test(NAME TestNomenclature COMMAND TestNomenclature)