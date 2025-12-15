cmake_minimum_required(VERSION 3.15)

project(toArp LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

add_subdirectory(external/AudioLib/Source)

find_package(PkgConfig REQUIRED)
pkg_check_modules(SNDFILE REQUIRED sndfile)

add_executable(toArp
  src/main.cpp
  src/find_peaks.cpp
)

target_include_directories(toArp PRIVATE
  ${SNDFILE_INCLUDE_DIRS}
  include
)

target_link_libraries(toArp
  Signal
  ${SNDFILE_LIBRARIES}
)

