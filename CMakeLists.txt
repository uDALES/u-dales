# uDALES (https://github.com/uDALES/u-dales).
# Copyright (C) 2019 D. Meyer.

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


cmake_minimum_required(VERSION 3.9)
project(uDALES Fortran)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake")

if (NOT CMAKE_BUILD_TYPE)
    message(STATUS "No build type selected, default to Release")
    set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Build type (default Release)" FORCE)
endif()

if(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fdefault-real-8 -ffree-line-length-none -cpp -ffpe-trap=inv,zero,ov")
    set(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} -fbounds-check -fno-unsafe-math-optimizations -frounding-math -fsignaling-nans -ggdb")
elseif(CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -r8 -fpp -WB -fpe0 -extend_source -mcmodel=medium")
    set(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} -g -traceback -warn unused -check bounds")
else()
    message(FATAL_ERROR "Only gfortran and ifort compilers supported")
endif()

find_package(MPI REQUIRED)

find_package(NetCDF REQUIRED COMPONENTS F90)
include_directories(${NETCDF_F90_INCLUDE_DIRS})

# FIXME: list files instead of using glob -- do during uDALES clean-up.
file(GLOB_RECURSE src_f77 "src/*.F")
file(GLOB_RECURSE src_f90 "src/*.f90")

add_executable(u-dales ${src_f77} ${src_f90})
target_link_libraries(u-dales MPI::MPI_Fortran)
target_link_libraries(u-dales ${NETCDF_F90_LIBRARIES})