cmake_minimum_required(VERSION 3.18)

project(uDALESPreprocessing LANGUAGES C Fortran)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib")
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib")

option(BUILD_PREPROCESSING_VIEW3D "Build the View3D executable" ON)
option(BUILD_PREPROCESSING_IBM "Build the IBM preprocessing executable" ON)
option(BUILD_PREPROCESSING_DIRECTSHORTWAVE_F2PY "Build the directshortwave f2py module" ON)
option(BUILD_PREPROCESSING_IBM_F2PY "Build the IBM preprocessing f2py module" ON)

set(PREPROCESSING_PYTHON_EXECUTABLE "" CACHE FILEPATH "Python interpreter used to build f2py preprocessing modules")

if(PREPROCESSING_PYTHON_EXECUTABLE)
  set(_python_executable "${PREPROCESSING_PYTHON_EXECUTABLE}")
else()
  find_package(Python3 REQUIRED COMPONENTS Interpreter)
  set(_python_executable "${Python3_EXECUTABLE}")
endif()

if(BUILD_PREPROCESSING_VIEW3D)
  add_subdirectory("${CMAKE_CURRENT_SOURCE_DIR}/../View3D" "${CMAKE_CURRENT_BINARY_DIR}/View3D")
endif()

if(BUILD_PREPROCESSING_IBM)
  add_executable(
    IBM_preproc
    "${CMAKE_CURRENT_SOURCE_DIR}/../python/fortran/ibm_preproc/in_mypoly_functions.f90"
    "${CMAKE_CURRENT_SOURCE_DIR}/../python/fortran/ibm_preproc/boundaryMasking.f90"
    "${CMAKE_CURRENT_SOURCE_DIR}/../python/fortran/ibm_preproc/matchFacetsCells.f90"
    "${CMAKE_CURRENT_SOURCE_DIR}/../python/fortran/ibm_preproc/IBM_preproc_io.f90"
    "${CMAKE_CURRENT_SOURCE_DIR}/../python/fortran/ibm_preproc/IBM_preproc_driver.f90"
    "${CMAKE_CURRENT_SOURCE_DIR}/../python/fortran/ibm_preproc/IBM_preproc_main.f90"
  )
  if(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
    target_compile_options(IBM_preproc PRIVATE -O3 -fopenmp)
    target_link_options(IBM_preproc PRIVATE -fopenmp)
  endif()
endif()

if(BUILD_PREPROCESSING_DIRECTSHORTWAVE_F2PY)
  execute_process(
    COMMAND "${_python_executable}" -c "import pathlib, sysconfig; header = pathlib.Path(sysconfig.get_paths().get('include', '')) / 'Python.h'; print(header); raise SystemExit(0 if header.is_file() else 1)"
    RESULT_VARIABLE _python_header_result
    OUTPUT_VARIABLE _python_header_path
    OUTPUT_STRIP_TRAILING_WHITESPACE
  )
  if(NOT _python_header_result EQUAL 0)
    message(FATAL_ERROR
      "Python development headers not found for ${_python_executable}. "
      "Expected Python.h at ${_python_header_path}.")
  endif()

  execute_process(
    COMMAND "${_python_executable}" -c "import sysconfig; print(sysconfig.get_config_var('EXT_SUFFIX') or '.so')"
    OUTPUT_VARIABLE _python_ext_suffix
    OUTPUT_STRIP_TRAILING_WHITESPACE
  )

  set(_directshortwave_target_dir "${CMAKE_CURRENT_SOURCE_DIR}/../python/udprep")
  set(_directshortwave_output "${_directshortwave_target_dir}/directshortwave_f2py${_python_ext_suffix}")
  set(_directshortwave_signature "${CMAKE_CURRENT_SOURCE_DIR}/../python/fortran/directshortwave_f2py.pyf")
  set(_directshortwave_source "${CMAKE_CURRENT_SOURCE_DIR}/../python/fortran/directShortwave.f90")
  set(_directshortwave_wrapper "${CMAKE_CURRENT_SOURCE_DIR}/../python/fortran/directshortwave_f2py_wrapper.f90")

  add_custom_command(
    OUTPUT "${_directshortwave_output}"
    COMMAND "${CMAKE_COMMAND}" -E make_directory "${_directshortwave_target_dir}"
    COMMAND "${CMAKE_COMMAND}" -E remove -f "${_directshortwave_target_dir}/directshortwave_f2py${_python_ext_suffix}"
    COMMAND "${_python_executable}" -m numpy.f2py -c
            "${_directshortwave_signature}"
            "${_directshortwave_source}"
            "${_directshortwave_wrapper}"
            --noopt
            --f90flags=-O3\ -ffree-line-length-none
    WORKING_DIRECTORY "${_directshortwave_target_dir}"
    DEPENDS "${_directshortwave_signature}" "${_directshortwave_source}" "${_directshortwave_wrapper}"
    COMMENT "Building directshortwave_f2py in ${_directshortwave_target_dir}"
    VERBATIM
  )

  add_custom_target(directshortwave_f2py ALL DEPENDS "${_directshortwave_output}")
endif()

if(BUILD_PREPROCESSING_IBM_F2PY)
  execute_process(
    COMMAND "${_python_executable}" -c "import sysconfig; print(sysconfig.get_config_var('EXT_SUFFIX') or '.so')"
    OUTPUT_VARIABLE _python_ext_suffix_ibm
    OUTPUT_STRIP_TRAILING_WHITESPACE
  )

  set(_ibm_f2py_target_dir "${CMAKE_CURRENT_SOURCE_DIR}/../python/udprep")
  set(_ibm_f2py_output "${_ibm_f2py_target_dir}/ibm_preproc_f2py${_python_ext_suffix_ibm}")
  set(_ibm_f2py_wrapper "${CMAKE_CURRENT_SOURCE_DIR}/../python/fortran/ibm_preproc_f2py.f90")
  set(_ibm_f2py_signature "${CMAKE_CURRENT_SOURCE_DIR}/../python/fortran/ibm_preproc_f2py.pyf")
  set(_ibm_f2py_sources
    "${CMAKE_CURRENT_SOURCE_DIR}/../python/fortran/ibm_preproc/in_mypoly_functions.f90"
    "${CMAKE_CURRENT_SOURCE_DIR}/../python/fortran/ibm_preproc/boundaryMasking.f90"
    "${CMAKE_CURRENT_SOURCE_DIR}/../python/fortran/ibm_preproc/matchFacetsCells.f90"
    "${CMAKE_CURRENT_SOURCE_DIR}/../python/fortran/ibm_preproc/IBM_preproc_io.f90"
    "${CMAKE_CURRENT_SOURCE_DIR}/../python/fortran/ibm_preproc/IBM_preproc_driver.f90"
  )

  add_custom_command(
    OUTPUT "${_ibm_f2py_output}"
    COMMAND "${CMAKE_COMMAND}" -E make_directory "${_ibm_f2py_target_dir}"
    COMMAND "${CMAKE_COMMAND}" -E remove -f "${_ibm_f2py_target_dir}/ibm_preproc_f2py${_python_ext_suffix_ibm}"
    COMMAND "${_python_executable}" -m numpy.f2py -c
            "${_ibm_f2py_signature}"
            "${CMAKE_CURRENT_SOURCE_DIR}/../python/fortran/ibm_preproc/in_mypoly_functions.f90"
            "${CMAKE_CURRENT_SOURCE_DIR}/../python/fortran/ibm_preproc/boundaryMasking.f90"
            "${CMAKE_CURRENT_SOURCE_DIR}/../python/fortran/ibm_preproc/matchFacetsCells.f90"
            "${CMAKE_CURRENT_SOURCE_DIR}/../python/fortran/ibm_preproc/IBM_preproc_io.f90"
            "${CMAKE_CURRENT_SOURCE_DIR}/../python/fortran/ibm_preproc/IBM_preproc_driver.f90"
            "${_ibm_f2py_wrapper}"
            --f90flags=-O3\ -fopenmp\ -ffree-line-length-none
            -lgomp
    WORKING_DIRECTORY "${_ibm_f2py_target_dir}"
    DEPENDS "${_ibm_f2py_signature}" "${_ibm_f2py_wrapper}" ${_ibm_f2py_sources}
    COMMENT "Building ibm_preproc_f2py in ${_ibm_f2py_target_dir}"
    VERBATIM
  )

  add_custom_target(ibm_preproc_f2py ALL DEPENDS "${_ibm_f2py_output}")
endif()

add_custom_target(preprocessing_tools)
if(TARGET view3d)
  add_dependencies(preprocessing_tools view3d)
endif()
if(TARGET IBM_preproc)
  add_dependencies(preprocessing_tools IBM_preproc)
endif()
if(TARGET directshortwave_f2py)
  add_dependencies(preprocessing_tools directshortwave_f2py)
endif()
if(TARGET ibm_preproc_f2py)
  add_dependencies(preprocessing_tools ibm_preproc_f2py)
endif()
