if($ENV{SYST} STREQUAL "localpc_gfortran")
  set(CMAKE_Fortran_COMPILER "mpif90")
  add_definitions(-W -Wall -fdefault-real-8 -ffree-line-length-none -I/include -I/usr/include -I/usr/local/include)
elseif($ENV{SYST}  STREQUAL "HUYGENS")
  set(CMAKE_Fortran_COMPILER "mpfort")
  set(CMAKE_Fortran_FLAGS "-O5 -qnoipa -qsmallstack=dynlenonheap -qfree=F90 -qrealsize=8 -qwarn64 -qflttrap=en:ov:zero:inv:imp -qflag=w:e")
elseif($ENV{SYST}  STREQUAL "localpc_ifort")
  set(CMAKE_Fortran_COMPILER "mpif90")
  set(CMAKE_Fortran_FLAGS "-r8 -ftz -fpe0 -extend_source -O3")
endif($ENV{SYST} STREQUAL "localpc_gfortran")
PROJECT(DALES Fortran)
cmake_minimum_required(VERSION 2.6)

set(CMAKE_BUILD_TYPE "Release")

set(VERSION_MAJOR "3")
set(VERSION_MINOR "1")
set(VERSION_PATCH "81")
find_path(NETCDF_INCLUDE_DIR netcdf.mod NETCDF.mod
  DOC "NetCDF include directory (must contain netcdf.mod)"
  PATHS
  /usr/include
  $ENV{HOME}/include
  $ENV{SARA_NETCDF_INCLUDE}
PATH_SUFFIXES netcdf-3)
find_library(NETCDF_C_LIB netcdf
  DOC "NetCDF C library"
  PATHS
  /usr/lib
  /usr/lib64
  $ENV{HOME}/lib
  $ENV{HOME}/lib64
  $ENV{SARA_NETCDF_LIB}
)
find_library(NETCDF_FORTRAN_LIB netcdff
  DOC "NetCDF Fortran library"
  PATHS
  /usr/lib
  /usr/lib64
  $ENV{HOME}/lib
  $ENV{HOME}/lib64
  $ENV{SARA_NETCDF_LIB}
)

set(NETCDF_LIBS ${NETCDF_C_LIB})
if(NETCDF_FORTRAN_LIB)
  set(NETCDF_LIBS ${NETCDF_LIBS} ${NETCDF_FORTRAN_LIB})
endif(NETCDF_FORTRAN_LIB)
###documentation
set(doc_dir doc)
find_package(Doxygen)
if(DOXYGEN)
  configure_file(${CMAKE_SOURCE_DIR}/Doxyfile.cmake Doxyfile)
  add_custom_target(docs
    ${DOXYGEN_EXECUTABLE} ${PROJECT_BINARY_DIR}/Doxyfile)
endif(DOXYGEN)

ADD_SUBDIRECTORY( src)
