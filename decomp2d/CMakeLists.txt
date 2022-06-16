#cmake_minimum_required(VERSION 3.9)
#project(decomp2d Fortran)

file(GLOB files_decomp decomp_2d.f90 glassman.f90 io.f90)

include_directories(${CMAKE_SOURCE_DIR}/decomp2d)

if(${FFT_Choice} MATCHES "generic")
  file(GLOB files_fft fft_generic.f90)
endif(${FFT_Choice} MATCHES "generic")

set(SRCFILES ${files_decomp} ${files_fft})

add_library(decomp2d STATIC ${SRCFILES})
if (MPI_FOUND)
  target_link_libraries(decomp2d PRIVATE MPI::MPI_Fortran)
endif (MPI_FOUND)
