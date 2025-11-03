cmake_minimum_required(VERSION 3.23)
project(LePIC LANGUAGES Fortran)

# ---------------- Compiler ----------------
# You can override on the command line:
#   cmake -S . -B build -DCMAKE_Fortran_COMPILER=mpiifx
set(CMAKE_Fortran_COMPILER mpiifx CACHE STRING "Fortran compiler (e.g., mpiifx)")

# ---------------- Build type & flags ----------------
set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)

# Release flags (fast; using MKL and OpenMP)
set(CMAKE_Fortran_FLAGS_RELEASE
    "-O3 -ipo -no-prec-div -m64 -xHost -qopenmp -auto -fp-model fast=2 -unroll -qmkl"
    CACHE STRING "Release flags" FORCE)

# Optional: a useful Debug profile
set(CMAKE_Fortran_FLAGS_DEBUG
    "-O0 -g -traceback -check all -check bounds -fpe0 -qopenmp -qmkl"
    CACHE STRING "Debug flags" FORCE)

# ---------------- Fortran standard & modules ----------------
set(CMAKE_Fortran_STANDARD 2008)
set(CMAKE_Fortran_STANDARD_REQUIRED ON)
set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/mod)

# ---------------- Output directories ----------------
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

# ---------------- Source collection ----------------
file(GLOB_RECURSE MODULE_SOURCES     CONFIGURE_DEPENDS "${CMAKE_SOURCE_DIR}/Src/modules/*.f90")
file(GLOB_RECURSE SUBMODULE_SOURCES  CONFIGURE_DEPENDS "${CMAKE_SOURCE_DIR}/Src/submodules/*.f90")
file(GLOB_RECURSE OTHER_SOURCES      CONFIGURE_DEPENDS "${CMAKE_SOURCE_DIR}/Src/*.f90")
list(REMOVE_ITEM OTHER_SOURCES ${MODULE_SOURCES})
list(REMOVE_ITEM OTHER_SOURCES ${SUBMODULE_SOURCES})
set(SRC_FILES ${MODULE_SOURCES} ${SUBMODULE_SOURCES} ${OTHER_SOURCES})

# ---------------- Target ----------------
add_executable(run_pic3D ${SRC_FILES})
target_include_directories(run_pic3D PRIVATE ${CMAKE_Fortran_MODULE_DIRECTORY})

# ---------------- HDF5 (Fortran) ----------------
# Prefix of HDF5 built with Intel ifx (change if you installed elsewhere).
set(HDF5_ROOT "/opt/hdf5-ifx" CACHE PATH "Prefix of HDF5 built with Intel ifx")
find_package(HDF5 1.14 COMPONENTS Fortran REQUIRED)

# .mod search path (hdf5.mod & friends)
target_include_directories(run_pic3D PRIVATE ${HDF5_Fortran_INCLUDE_DIRS})

# Link HDF5 Fortran (and underlying C libs)
target_link_libraries(run_pic3D PRIVATE ${HDF5_Fortran_LIBRARIES} ${HDF5_LIBRARIES})

# ---------------- oneAPI + HDF5 RPATH setup ----------------
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
    "$ENV{ONEAPI_ROOT}/mpi/latest/libfabric/lib"
  )
else()
  list(APPEND _ONEAPI_RPATH
    "/opt/intel/oneapi/compiler/latest/lib"
    "/opt/intel/oneapi/mpi/latest/lib/release"
    "/opt/intel/oneapi/mpi/latest/libfabric/lib"
  )
endif()

# Add HDF5 lib dir to RPATH so binary runs without LD_LIBRARY_PATH
list(APPEND _ONEAPI_RPATH "${HDF5_Fortran_LIBRARY_DIRS}")
list(APPEND _ONEAPI_RPATH "${HDF5_LIBRARY_DIRS}")
list(APPEND _ONEAPI_RPATH "${HDF5_DIR}/lib")  # harmless if empty

# Apply RPATH
set_target_properties(run_pic3D PROPERTIES
  BUILD_RPATH   "${_ONEAPI_RPATH}"
  INSTALL_RPATH "${_ONEAPI_RPATH}"
)
target_link_options(run_pic3D PRIVATE -Wl,-rpath,$<JOIN:${_ONEAPI_RPATH},:>)
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH ON)

# ---------------- Convenience run target ----------------
add_custom_target(run
  COMMAND mpirun -np 4 $<TARGET_FILE:run_pic3D>
  DEPENDS run_pic3D
  WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
  COMMENT "Running PIC with 4 MPI ranks (cwd = source root)"
)

# ---------------- Summary message ----------------
message(STATUS "Fortran compiler: ${CMAKE_Fortran_COMPILER}")
message(STATUS "Build type      : ${CMAKE_BUILD_TYPE}")
message(STATUS "RPATH paths     : ${_ONEAPI_RPATH}")
message(STATUS "HDF5 Fortran includes: ${HDF5_Fortran_INCLUDE_DIRS}")
message(STATUS "HDF5 Fortran libs    : ${HDF5_Fortran_LIBRARIES}")
message(STATUS "HDF5 C libs          : ${HDF5_LIBRARIES}")
