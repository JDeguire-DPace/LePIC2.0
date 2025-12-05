cmake_minimum_required(VERSION 3.23)
project(LePIC LANGUAGES Fortran)

# --- Compiler (override with -DCMAKE_Fortran_COMPILER=mpifort/gfortran/ifx) ---
# Default to mpifort so MPI modules/flags are handled automatically.
set(CMAKE_Fortran_COMPILER mpifort CACHE STRING "Fortran compiler (mpifort/gfortran/etc.)")

# --- Build type & flags ---
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Debug CACHE STRING "Build type" FORCE)
endif()

set(CMAKE_Fortran_FLAGS_RELEASE "-O3 -march=native -funroll-loops -fopenmp"
    CACHE STRING "" FORCE)
set(CMAKE_Fortran_FLAGS_DEBUG   "-O0 -g -fcheck=all -fbacktrace -fopenmp"
    CACHE STRING "" FORCE)

# --- Fortran standard & output dirs ---
set(CMAKE_Fortran_STANDARD 2008)
set(CMAKE_Fortran_STANDARD_REQUIRED ON)
set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/mod)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

# --- Sources ---
file(GLOB_RECURSE MODULE_SOURCES    CONFIGURE_DEPENDS "${CMAKE_SOURCE_DIR}/Src/modules/*.f90")
file(GLOB_RECURSE OTHER_SOURCES     CONFIGURE_DEPENDS "${CMAKE_SOURCE_DIR}/Src/*.f90")
list(REMOVE_ITEM OTHER_SOURCES ${MODULE_SOURCES})

add_executable(run_pic3D ${MODULE_SOURCES} ${OTHER_SOURCES})
target_include_directories(run_pic3D PRIVATE ${CMAKE_Fortran_MODULE_DIRECTORY})

target_compile_options(run_pic3D PRIVATE
  $<$<COMPILE_LANGUAGE:Fortran>:
    -ffree-line-length-none
    -Wno-line-truncation
    -Wno-error=line-truncation>
)

# --- MPI (Fortran) ---
# Works whether you use mpifort or plain gfortran.
find_package(MPI REQUIRED Fortran)
target_link_libraries(run_pic3D PRIVATE MPI::MPI_Fortran)

# --- HDF5 (Fortran) ---
# Tip: if needed, pass a custom prefix: -DHDF5_ROOT=/usr/lib/x86_64-linux-gnu/hdf5/gfortran
set(HDF5_ROOT "" CACHE PATH "Prefix of HDF5 with Fortran libs (optional)")
find_package(HDF5 COMPONENTS Fortran REQUIRED)
target_include_directories(run_pic3D PRIVATE ${HDF5_Fortran_INCLUDE_DIRS})
target_link_libraries(run_pic3D  PRIVATE ${HDF5_Fortran_LIBRARIES} ${HDF5_LIBRARIES})

# --- RPATH (so the binary runs without manual LD_LIBRARY_PATH) ---
set(_RPATH_CANDIDATES
  ${HDF5_Fortran_LIBRARY_DIRS}
  ${HDF5_LIBRARY_DIRS}
  "$ENV{MPI_HOME}/lib"
  "/usr/lib/x86_64-linux-gnu/openmpi/lib"        # common OpenMPI path on Debian/Ubuntu
  "/usr/lib/x86_64-linux-gnu"                    # generic fallback
)

# Drop empty entries and apply
set(_RPATH "")
foreach(p IN LISTS _RPATH_CANDIDATES)
  if(p)
    list(APPEND _RPATH "${p}")
  endif()
endforeach()

set_target_properties(run_pic3D PROPERTIES
  BUILD_RPATH   "${_RPATH}"
  INSTALL_RPATH "${_RPATH}")
target_link_options(run_pic3D PRIVATE -Wl,-rpath,$<JOIN:${_RPATH},:>)
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
message(STATUS "MPI                   : ${MPI_Fortran_FOUND}")
message(STATUS "HDF5 Fortran includes : ${HDF5_Fortran_INCLUDE_DIRS}")
message(STATUS "HDF5 Fortran libs     : ${HDF5_Fortran_LIBRARIES}")
message(STATUS "Extra RPATH paths     : ${_RPATH}")
