if (NOT CMAKE_BUILD_TYPE)
  set (CMAKE_BUILD_TYPE RELEASE CACHE STRING
      "Choose the type of build, options are: None Debug Release."
      FORCE)
endif (NOT CMAKE_BUILD_TYPE)

if($ENV{SYST} STREQUAL "")
  set($ENV{SYST} STREQUAL "localpc_gfortran")
endif($ENV{SYST} STREQUAL "")

if($ENV{SYST} STREQUAL "localpc_gfortran")
  set(CMAKE_Fortran_COMPILER "mpif90")
  set(CMAKE_Fortran_FLAGS "-finit-real=nan -W -Wall -fdefault-real-8 -ffree-line-length-none -I/include -I/usr/include -I/usr/local/include" CACHE STRING "")
  set (CMAKE_Fortran_FLAGS_RELEASE "-funroll-all-loops -fno-f2c -O3")
  set (CMAKE_Fortran_FLAGS_DEBUG   "-fbounds-check -fbacktrace -fno-f2c -O0 -g -ffpe-trap=invalid,zero,overflow")
elseif($ENV{SYST}  STREQUAL "HUYGENS")
  set(CMAKE_Fortran_COMPILER "mpfort")
  set(CMAKE_Fortran_FLAGS "-qextname -qmsallstack=dynlenonheap -qfree=F90 -q64 -qarch=auto -qmaxmem=-1 -qrealsize=8  -qwarn64 -qflttrap=en:ov:zero:inv:imp -qflag=w:e -I$ENV{HOME}/include " CACHE STRING "")
  set(CMAKE_Fortran_FLAGS_RELEASE "-O5 -qnoipa -qtune=auto " CACHE STRING "")
  set(CMAKE_Fortran_FLAGS_DEBUG "-O2 -g -qfullpath -C -qsigtrap -qinitauto=ff" CACHE STRING "")
elseif($ENV{SYST}  STREQUAL "localpc_ifort")
  set(CMAKE_Fortran_COMPILER "mpif90")
  set(CMAKE_Fortran_FLAGS "-r8 -ftz -fpe0 -extend_source" CACHE STRING "")
  set (CMAKE_Fortran_FLAGS_RELEASE "-O3 -pg")
  set (CMAKE_Fortran_FLAGS_DEBUG   "-O0 -g -ffpe-trap=invalid,zero,overflow")
endif($ENV{SYST} STREQUAL "localpc_gfortran")
PROJECT(DALES Fortran)
cmake_minimum_required(VERSION 2.6)

set(VERSION_MAJOR "3")
set(VERSION_MINOR "1")
set(VERSION_PATCH "82")
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

if (CASE)
    FILE(GLOB usrfile "${CMAKE_CURRENT_SOURCE_DIR}/cases/${CASE}/moduser.f90")
    FILE(REMOVE ${CMAKE_CURRENT_SOURCE_DIR}/src/moduser.f90)
    if (usrfile STREQUAL "")
      set(usrfile "${CMAKE_CURRENT_SOURCE_DIR}/cases/standard/moduser.f90")
    endif (usrfile STREQUAL "")
      CONFIGURE_FILE(${usrfile} ${CMAKE_CURRENT_SOURCE_DIR}/src/moduser.f90 COPYONLY)
      MESSAGE(STATUS "Case " ${CASE} " uses " ${usrfile})
else(CASE)
   FILE(GLOB usrfile "${CMAKE_CURRENT_SOURCE_DIR}/src/moduser.f90")
   if (usrfile STREQUAL "")
     CONFIGURE_FILE( ${CMAKE_CURRENT_SOURCE_DIR}/cases/standard/moduser.f90 ${CMAKE_CURRENT_SOURCE_DIR}/src/moduser.f90 COPYONLY)
    endif (usrfile STREQUAL "")
endif(CASE)
#unset(CASE CACHE)


ADD_SUBDIRECTORY( src)
