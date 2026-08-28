cmake_minimum_required(VERSION 3.16)
project(plasma_prep LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -march=native -DNDEBUG")
set(CMAKE_CXX_FLAGS_DEBUG "-O0 -g -Wall -Wextra")

find_package(Eigen3 3.3 REQUIRED NO_MODULE)
find_package(PkgConfig REQUIRED)
pkg_check_modules(FFTW3 REQUIRED IMPORTED_TARGET fftw3)

add_library(core INTERFACE)
target_include_directories(core INTERFACE ${CMAKE_SOURCE_DIR}/include)
target_link_libraries(core INTERFACE Eigen3::Eigen)

file(GLOB_RECURSE FD_SOURCES CONFIGURE_DEPENDS "${CMAKE_SOURCE_DIR}/lib/fd/*.cpp")
add_library(fd STATIC ${FD_SOURCES})
target_include_directories(fd PUBLIC ${CMAKE_SOURCE_DIR}/include)
target_link_libraries(fd PUBLIC core)

file(GLOB EXPERIMENT_SOURCES CONFIGURE_DEPENDS "${CMAKE_SOURCE_DIR}/src/*.cpp")

foreach(src_file ${EXPERIMENT_SOURCES})
  get_filename_component(exe_name ${src_file} NAME_WE)
  add_executable(${exe_name} ${src_file})
  target_include_directories(${exe_name} PRIVATE ${CMAKE_SOURCE_DIR}/include)
  if(exe_name MATCHES "^01_")
    target_link_libraries(${exe_name} PRIVATE fd)
  elseif(exe_name MATCHES "^00_")
    target_link_libraries(${exe_name} PRIVATE Eigen3::Eigen PkgConfig::FFTW3)
  endif()
endforeach()

enable_testing()