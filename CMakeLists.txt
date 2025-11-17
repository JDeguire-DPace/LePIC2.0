cmake_minimum_required(VERSION 3.23)
project(LePIC LANGUAGES Fortran)

# --- Compiler (override with -DCMAKE_Fortran_COMPILER=ifx, gfortran, etc.) ---
set(CMAKE_Fortran_COMPILER mpiifx CACHE STRING "Fortran compiler")

# --- Build type & flags ---
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Debug CACHE STRING "Build type" FORCE)
endif()

set(CMAKE_Fortran_FLAGS_RELEASE "-O3 -ipo -no-prec-div -m64 -xHost -qopenmp -auto -fp-model fast=2 -unroll -qmkl"
    CACHE STRING "" FORCE)
set(CMAKE_Fortran_FLAGS_DEBUG   "-O0 -g -traceback -check all -check bounds -fpe0 -qopenmp -qmkl"
    CACHE STRING "" FORCE)

# --- Fortran standard & output dirs ---
set(CMAKE_Fortran_STANDARD 2008)
set(CMAKE_Fortran_STANDARD_REQUIRED ON)
set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/mod)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

# --- Sources ---
file(GLOB_RECURSE MODULE_SOURCES    CONFIGURE_DEPENDS "${CMAKE_SOURCE_DIR}/Src/modules/*.f90")
file(GLOB_RECURSE SUBMODULE_SOURCES CONFIGURE_DEPENDS "${CMAKE_SOURCE_DIR}/Src/submodules/*.f90")
file(GLOB_RECURSE OTHER_SOURCES     CONFIGURE_DEPENDS "${CMAKE_SOURCE_DIR}/Src/*.f90")
list(REMOVE_ITEM OTHER_SOURCES ${MODULE_SOURCES} ${SUBMODULE_SOURCES})
add_executable(run_pic3D ${MODULE_SOURCES} ${SUBMODULE_SOURCES} ${OTHER_SOURCES})
target_include_directories(run_pic3D PRIVATE ${CMAKE_Fortran_MODULE_DIRECTORY})

# --- HDF5 (Fortran) ---
# Tip: set HDF5_ROOT if needed: -DHDF5_ROOT=/opt/hdf5-ifx
set(HDF5_ROOT "/opt/hdf5-ifx" CACHE PATH "HDF5 prefix (ifx build)")
find_package(HDF5 1.14 COMPONENTS Fortran REQUIRED)
target_include_directories(run_pic3D PRIVATE ${HDF5_Fortran_INCLUDE_DIRS})
target_link_libraries(run_pic3D  PRIVATE ${HDF5_Fortran_LIBRARIES} ${HDF5_LIBRARIES})

# --- oneAPI + HDF5 RPATH (run without LD_LIBRARY_PATH) ---
set(_ONEAPI_RPATH "")
if(DEFINED ENV{MKLROOT})
  list(APPEND _ONEAPI_RPATH "$ENV{MKLROOT}/lib/intel64")
else()
  list(APPEND _ONEAPI_RPATH "/opt/intel/oneapi/mkl/latest/lib/intel64")
endif()

if(DEFINED ENV{ONEAPI_ROOT})
  list(APPEND _ONEAPI_RPATH
    "$ENV{ONEAPI_ROOT}/compiler/latest/lib"
    "$ENV{ONEAPI_ROOT}/mpi/latest/lib/release"
    "$ENV{ONEAPI_ROOT}/mpi/latest/libfabric/lib")
else()
  list(APPEND _ONEAPI_RPATH
    "/opt/intel/oneapi/compiler/latest/lib"
    "/opt/intel/oneapi/mpi/latest/lib/release"
    "/opt/intel/oneapi/mpi/latest/libfabric/lib")
endif()

list(APPEND _ONEAPI_RPATH
  ${HDF5_Fortran_LIBRARY_DIRS}
  ${HDF5_LIBRARY_DIRS}
  "${HDF5_DIR}/lib")

set_target_properties(run_pic3D PROPERTIES
  BUILD_RPATH   "${_ONEAPI_RPATH}"
  INSTALL_RPATH "${_ONEAPI_RPATH}")
target_link_options(run_pic3D PRIVATE -Wl,-rpath,$<JOIN:${_ONEAPI_RPATH},:>)
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH ON)

# --- Convenience: mpirun -np 4 ---
add_custom_target(run
  COMMAND mpirun -np 4 $<TARGET_FILE:run_pic3D>
  DEPENDS run_pic3D
  WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
  COMMENT "Running with 4 MPI ranks")

# --- Summary ---
message(STATUS "Fortran compiler      : ${CMAKE_Fortran_COMPILER}")
message(STATUS "Build type            : ${CMAKE_BUILD_TYPE}")
message(STATUS "RPATH paths           : ${_ONEAPI_RPATH}")
message(STATUS "HDF5 Fortran includes : ${HDF5_Fortran_INCLUDE_DIRS}")
message(STATUS "HDF5 Fortran libs     : ${HDF5_Fortran_LIBRARIES}")
message(STATUS "HDF5 C libs           : ${HDF5_LIBRARIES}")
