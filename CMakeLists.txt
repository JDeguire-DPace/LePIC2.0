cmake_minimum_required(VERSION 3.23)
project(LePIC LANGUAGES Fortran)

# --- Compiler (or pass -DCMAKE_Fortran_COMPILER=mpiifx) ---
set(CMAKE_Fortran_COMPILER mpiifx)

# --- Build type & flags ---
set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
set(CMAKE_Fortran_FLAGS_RELEASE "-O3 -ipo -no-prec-div -m64 -xHost -qopenmp -auto -fp-model fast=2 -unroll -qmkl")

# --- Fortran standard & .mod output ---
set(CMAKE_Fortran_STANDARD 2008)
set(CMAKE_Fortran_STANDARD_REQUIRED ON)
set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/mod)

# --- Keep executables inside the build tree ---
# If you use a multi-config generator (e.g., Visual Studio), add per-config dirs:
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

# --- Sources ---
file(GLOB_RECURSE MODULE_SOURCES   CONFIGURE_DEPENDS "${CMAKE_SOURCE_DIR}/Src/modules/*.f90")
file(GLOB_RECURSE SUBMODULE_SOURCES CONFIGURE_DEPENDS "${CMAKE_SOURCE_DIR}/Src/submodules/*.f90")
file(GLOB_RECURSE OTHER_SOURCES    CONFIGURE_DEPENDS "${CMAKE_SOURCE_DIR}/Src/*.f90")

# De-dup
list(REMOVE_ITEM OTHER_SOURCES ${MODULE_SOURCES})
list(REMOVE_ITEM OTHER_SOURCES ${SUBMODULE_SOURCES})

set(SRC_FILES ${MODULE_SOURCES} ${SUBMODULE_SOURCES} ${OTHER_SOURCES})

# --- Target ---
add_executable(run_pic3D ${SRC_FILES})
target_include_directories(run_pic3D PRIVATE ${CMAKE_Fortran_MODULE_DIRECTORY})

# --- (Optional) MPI run convenience target, no copying to root ---
# Working dir is the project root so relative input files still work.
add_custom_target(run
  COMMAND mpirun -np 4 $<TARGET_FILE:run_pic3D>
  DEPENDS run_pic3D
  WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
  COMMENT "Running PIC with 4 MPI ranks from build/bin (cwd = source root)"
)
