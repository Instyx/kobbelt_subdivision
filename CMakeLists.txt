cmake_minimum_required(VERSION 3.6..4.1)
project(kobbelt_subdivision)


# Path definitions:
# option(LIBIGL_WITH_TRIANGLE "Use Triangle for constrained DT" ON)
list(PREPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)
set(EXTERNAL_DIR "external")
# Configure general CMake variables
set(CMAKE_DEBUG_POSTFIX d)
set(CMAKE_BUILD_TYPE Debug)
#set(CMAKE_POSITION_INDEPENDENT_CODE ON)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON CACHE INTERNAL "")
set(CMAKE_POSITION_INDEPENDENT_CODE ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)


if (NOT EXISTS ${CMAKE_BINARY_DIR}/CMakeCache.txt)
  if (NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Release" CACHE STRING "" FORCE)
  endif()
endif()




add_subdirectory(${EXTERNAL_DIR}/geometry-central)

add_subdirectory(${EXTERNAL_DIR}/polyscope)

add_subdirectory(${EXTERNAL_DIR}/libigl)
# find_package(Boost REQUIRED)
# igl_include(predicates)
# igl_include(restricted triangle)

file(GLOB SRC_FILES src/*.cpp include/*.hpp)
add_executable(${PROJECT_NAME} ${SRC_FILES})

target_link_libraries(${PROJECT_NAME} PUBLIC
  igl::core
  # igl::predicates
  # igl_restricted::triangle
  polyscope
  geometry-central
  # Boost::boost
)
