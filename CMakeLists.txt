cmake_minimum_required(VERSION 3.23)
project(LePIC LANGUAGES Fortran)

# --- Compiler (you can also pass -DCMAKE_Fortran_COMPILER=mpiifx in your config) ---
set(CMAKE_Fortran_COMPILER mpiifx)

# --- Build type & flags (ifx/mpiifx) ---
set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
set(CMAKE_Fortran_FLAGS_RELEASE "-O3 -ipo -no-prec-div -m64 -xHost -qopenmp -auto -fp-model fast=2 -unroll -qmkl")

# --- Fortran standard & .mod output ---
set(CMAKE_Fortran_STANDARD 2008)
set(CMAKE_Fortran_STANDARD_REQUIRED ON)
set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/mod)

# --- Safer: build exe in build/; we'll copy it to repo root after link ---
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})

# --- Sources (your Src/ tree with modules & submodules & main.f90) ---
file(GLOB_RECURSE SRC_FILES
     "${CMAKE_SOURCE_DIR}/Src/*.f90")

# --- Target ---
add_executable(run_pic3D ${SRC_FILES})


# --- After building, copy the exe to the project root so you can run ./run_pic3D ---
add_custom_command(TARGET run_pic3D POST_BUILD
  COMMAND ${CMAKE_COMMAND} -E copy_if_different
          "$<TARGET_FILE:run_pic3D>"
          "${CMAKE_SOURCE_DIR}/run_pic3D"
  COMMENT "Copying run_pic3D to project root"
)

# --- Convenience 'run' target (MPI 4 ranks as example) ---
add_custom_target(run
  COMMAND mpirun -np 4 "${CMAKE_SOURCE_DIR}/run_pic3D"
  DEPENDS run_pic3D
  WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
  COMMENT "Running PIC with 4 MPI ranks"
)
