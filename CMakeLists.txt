cmake_minimum_required(VERSION 3.23)
project(LePIC LANGUAGES Fortran)

# Choose compiler (if not already set in environment)
set(CMAKE_Fortran_COMPILER mpiifx)   # or mpiifort / ifx / ifort

# Optimization flags (analogous to your C++ line)
set(CMAKE_BUILD_TYPE Release)
set(CMAKE_Fortran_FLAGS_RELEASE "-O3 -xHost -qopenmp -ipo -fast -qmkl -Wall")

# Fortran standard
set(CMAKE_Fortran_STANDARD 2008)
set(CMAKE_Fortran_STANDARD_REQUIRED ON)

# Add include directory for .mod files
set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/mod)

# Gather all sources automatically
file(GLOB_RECURSE SRC_FILES src/*.F90)

# Create the executable
add_executable(pic3D ${SRC_FILES})

# Link to libraries (optional)
find_package(MPI REQUIRED)
find_package(OpenMP REQUIRED)

target_link_libraries(pic3D PRIVATE MPI::MPI_Fortran OpenMP::OpenMP_Fortran)

# If you use MKL
# target_link_options(pic1d PRIVATE -qmkl)
