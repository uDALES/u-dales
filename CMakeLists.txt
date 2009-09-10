if($ENV{SYST} STREQUAL "localpc_gfortran")
  set(CMAKE_Fortran_COMPILER "mpif90")
  set(CMAKE_Fortran_FLAGS "-W -Wall -fdefault-real-8 -ffree-line-length-none -I/include -I/usr/include -I/usr/local/include" CACHE STRING "")
elseif($ENV{SYST}  STREQUAL "HUYGENS")
  set(CMAKE_Fortran_COMPILER "mpfort")
  set(CMAKE_Fortran_FLAGS "-O5 -qnoipa -qsmallstack=dynlenonheap -qfree=F90 -qrealsize=8 -qwarn64 -qflttrap=en:ov:zero:inv:imp -qflag=w:e" CACHE STRING "")
elseif($ENV{SYST}  STREQUAL "localpc_ifort")
  set(CMAKE_Fortran_COMPILER "mpif90")
  set(CMAKE_Fortran_FLAGS "-r8 -ftz -fpe0 -extend_source -O3" CACHE STRING "")
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

if(NETCDF_INCLUDE_DIR)
else(NETCDF_INCLUDE_DIR)
  MESSAGE(STATUS "WARNING: No NETCDF bindings are found.")
endif(NETCDF_INCLUDE_DIR)

if(NETCDF_C_LIB)
  set(NETCDF_LIBS ${NETCDF_C_LIB})
else(NETCDF_C_LIB)
  MESSAGE(STATUS "WARNING: No NETCDF bindings are found.")
endif(NETCDF_C_LIB)

if(NETCDF_FORTRAN_LIB)
  set(NETCDF_LIBS ${NETCDF_LIBS} ${NETCDF_FORTRAN_LIB})
else(NETCDF_FORTRAN_LIB)
  MESSAGE(STATUS "WARNING: No Fortran NETCDF bindings are found.")
endif(NETCDF_FORTRAN_LIB)
###documentation
INCLUDE(FindDoxygen)
IF (DOXYGEN)
    ADD_SUBDIRECTORY(utils/doc)
ELSE (DOXYGEN)
    MESSAGE(STATUS "WARNING: Doxygen not found - Reference manual will 
not be created")
ENDIF (DOXYGEN)


ADD_SUBDIRECTORY( src)
