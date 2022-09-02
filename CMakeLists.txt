cmake_minimum_required(VERSION 3.5)
project(MC3D)
set(CMAKE_C_FLAGS_RELEASE "-DNDEBUG -O3")
set(CMAKE_C_FLAGS_DEBUG "-Wall -g")

# This will allow the user to choose to compile the MC3D with UMFPACK or PARDISO
# Options are "PARDISO" or "UMFPACK" - all caps
set(SPARSE_SOLVER "PARDISO")

# This file finds the SuiteSparse UMFPack installation
# and adds all of the subdirectories in the MC3D project.
#
# Created by:
#
# Ilya Prigozhin
# ilyap@bu.edu
# Computational Electronics Group
# Boston University
##################################################################################
# Guidelines to be followed by MC3D version numbers:
# MAJOR: incremented when changes are made that affect whether existing
#        input files are compatible with the new release
# MINOR: incremented when functionality is added in a backwards compatible manner
# PATCH: incremented when bug fixes are implemented in a backwards compatible manner
set(MAJOR_RELEASE 2)
set(MINOR_RELEASE 13)
set(PATCH_RELEASE 0)
##################################################################################
# Set CMake module and output binary paths.
# CMake modules will go in "Source/cmake/Modules/"
# Binary output will go in "build/bin"
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules")
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set(CMAKE_VERBOSE_MAKEFILE ON)

# Add general compiler options for 64-bit compile
# and C11 standard with GNU extensions
# NOTE: This command only works if it comes before the target is added
add_compile_options(-m64 -std=gnu11)
##################################################################################
if(${SPARSE_SOLVER} MATCHES "UMFPACK")
	message("*****Compiling with UMFPACK*****")
	add_compile_definitions(UMF_SPARSE_SOLVER)
	find_package(UMFPACK)
	# List of modules that we want to include
	set(dirs SparseSolverUmf	SparseMatrixUmf)
elseif(${SPARSE_SOLVER} MATCHES "PARDISO")
    	message("*****Compiling with Pardiso*****")
	add_compile_definitions(PARDISO_SPARSE_SOLVER)
	find_package(PARDISO)
	# List of modules that we want to include
	set(dirs SparseSolverPardiso SparseMatrixPardiso)
else(${SPARSE_SOLVER} MATCHES "UMFPACK")
    	message("Sparse Solver input not recognized, compiling with UMFPACK")
	add_compile_definitions(PARDISO_SPARSE_SOLVER)
	find_package(PARDISO)
	# List of modules that we want to include
	set(dirs SparseSolverPardiso SparseMatrixPardiso)
endif(${SPARSE_SOLVER} MATCHES "UMFPACK")

set(dirs ${dirs} Control Version RndGen List Reader SpecFun
	Octree NLSolver
	Geometry Algebra Particle Material
	AnalyticBand TPMesh FullBand Scattering
	Poisson MEMC DEMC Main
)

# Find OpenMP installation and set compiler flags
include(FindOpenMP)
if(OPENMP_FOUND)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
else(OPENMP_FOUND)
    message("ERROR: OpenMP could not be found.")
endif(OPENMP_FOUND)

##################################################################################
# Create function to add all of the dirs to the project as modules
function(add_dirs)
    foreach(arg IN LISTS ARGV)
        message(STATUS "Adding directory: ${arg}")
        add_subdirectory(${arg})
    endforeach()
endfunction()

# Now load the library modules
add_dirs(${dirs})

##################################################################################
if(CMAKE_BUILD_TYPE MATCHES "Release")
    add_executable(MC3D_release Main/src/main.c)
    target_link_libraries(MC3D_release control octree scattering
                                particle demc memc dfise 
                                gmsh specfun version)
elseif(CMAKE_BUILD_TYPE MATCHES "Debug")
    add_executable(MC3D_debug Main/src/main.c)
    target_link_libraries(MC3D_debug control octree scattering 
                                particle demc memc dfise
                                gmsh specfun version)
endif()
