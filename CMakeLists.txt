# uDALES (https://github.com/uDALES/u-dales).

# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

# Copyright (C) 2019 the uDALES Team.

cmake_minimum_required(VERSION 3.9)
project(uDALES Fortran)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake")

if (NOT CMAKE_BUILD_TYPE)
    message(STATUS "No build type selected, default to Release")
    set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Build type (default Release)" FORCE)
endif()

set(SKIP_UPDATE_EXTERNAL_PROJECTS OFF CACHE STRING "Update external projects")
if (SKIP_UPDATE_EXTERNAL_PROJECTS)
    message(STATUS "Skipping updating external projects")
endif()

find_package(MPI REQUIRED)

if(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
    set(CMAKE_Fortran_FLAGS "-g -fbacktrace -fdefault-real-8 -ffree-line-length-none -ffpe-trap=invalid,zero,overflow -fallow-argument-mismatch")
    set(CMAKE_Fortran_FLAGS_DEBUG "-O0 -finit-real=nan -fcheck=all -Wall -Wextra -Wuninitialized -Warray-bounds -Wconversion")
    set(CMAKE_Fortran_FLAGS_RELEASE "-O3")
elseif(CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
    # FIXME: `-heap-arrays 10` required as using Intel Fortran
    # means they are stored on the stack (by default), whereas GNU Fortran stores
    # them on the heap (https://github.com/uDALES/u-dales/issues/13).
    set(CMAKE_Fortran_FLAGS "-g -traceback -r8 -fpe0 -heap-arrays 10")
    set(CMAKE_Fortran_FLAGS_DEBUG "-O0 -init=snan -CB -check all")
    set(CMAKE_Fortran_FLAGS_RELEASE "-O3")
elseif(CMAKE_Fortran_COMPILER_ID MATCHES "Cray")
    # https://pubs.cray.com/content/S-3901/8.5/cray-fortran-reference-manual-85
    set(CMAKE_Fortran_FLAGS "-G0 -s real64 -N 1023 -K trap=divz,inv,ovf")
    set(CMAKE_Fortran_FLAGS_DEBUG "-O0 -ei -R bcdps -m 0")
    set(CMAKE_Fortran_FLAGS_RELEASE "-O3")
else()
    message(FATAL_ERROR "Only GNU, Intel, and Cray Fortran compilers are supported")
endif()

#set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/findFFTW")
configure_file(downloadFindFFTW.cmake.in findFFTW-download/CMakeLists.txt)
execute_process(COMMAND ${CMAKE_COMMAND} -G "${CMAKE_GENERATOR}" .
        RESULT_VARIABLE result
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/findFFTW-download )
if(result)
    message(FATAL_ERROR "CMake step for findFFTW failed: ${result}")
    else()
    message("CMake step for findFFTW completed (${result}).")
endif()
execute_process(COMMAND ${CMAKE_COMMAND} --build .
        RESULT_VARIABLE result
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/findFFTW-download )
if(result)
    message(FATAL_ERROR "Build step for findFFTW failed: ${result}")
endif()
set(findFFTW_DIR ${CMAKE_CURRENT_BINARY_DIR}/findFFTW-src)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${findFFTW_DIR}")
find_package(FFTW REQUIRED)
include_directories(${FFTW_INCLUDE_DIRS})

find_package(NetCDF REQUIRED COMPONENTS F90)
include_directories(${NETCDF_F90_INCLUDE_DIRS})

file(GLOB_RECURSE src_f90 "src/*.f90")

#add_subdirectory(external)

#ExternalProject_Get_Property(fishpack-cmake-project binary_dir)
#set(FISHPACK_LIBRARY ${binary_dir}/src/${CMAKE_STATIC_LIBRARY_PREFIX}fishpack_static${CMAKE_STATIC_LIBRARY_SUFFIX})

#ExternalProject_Get_Property(vfftpack-cmake-project binary_dir)
#set(VFFTPACK_LIBRARY ${binary_dir}/${CMAKE_STATIC_LIBRARY_PREFIX}vfftpack_static${CMAKE_STATIC_LIBRARY_SUFFIX})

option(DOUBLE_PRECISION "Build 2decomp_fft with double precision" ON)
if (DOUBLE_PRECISION)
  add_definitions("-DDOUBLE_PREC")
endif()

option(SINGLE_PRECISION_OUTPUT "Build 2decomp_fft with output in single precision" OFF)
if (SINGLE_PRECISION_OUTPUT)
  add_definitions("-DSAVE_SINGLE")
endif()

set(FFT_Choice "fftw3" CACHE STRING "FFT for 2decomp_fft (generic is the default)")
set_property(CACHE FFT_Choice PROPERTY STRINGS generic fftw3)

add_subdirectory(2decomp-fft)

add_executable(u-dales ${src_f90})
include_directories(2decomp-fft)

target_link_libraries(u-dales ${CMLIB_CMAKE_LIBRARY})
target_link_libraries(u-dales MPI::MPI_Fortran)
target_link_libraries(u-dales ${NETCDF_F90_LIBRARIES})
#target_link_libraries(u-dales ${FISHPACK_LIBRARY})
#target_link_libraries(u-dales ${VFFTPACK_LIBRARY})
target_link_libraries(u-dales ${FFTW_LIBRARIES})
target_link_libraries(u-dales 2decomp-fft)

target_include_directories(u-dales PUBLIC $<TARGET_PROPERTY:2decomp-fft,BINARY_DIR>)
target_include_directories(u-dales PUBLIC ${FFTW_INCLUDE_DIRS})

#add_dependencies(u-dales fishpack-cmake-project)
#add_dependencies(u-dales vfftpack-cmake-project)
add_dependencies(u-dales 2decomp-fft)
