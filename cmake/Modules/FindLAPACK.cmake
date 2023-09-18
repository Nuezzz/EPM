# Finds the intel math kernel library (MKL) required to run pardiso
# also finds the intel packages required to perform sparse matrix operations
# for csc and coo matrices

#Check is MKLROOT environment variable is set, if not set a local variable with MKLROOT
#If compiling outside of the BU CSM, will need to set MKLROOT Env variable
if(DEFINED ENV{MKLROOT})
	set(MKLROOT $ENV{MKLROOT})
else(DEFINED ENV{MKLROOT})
	set(MKLROOT "/mnt/csmhome/software/intel/mkl")
endif(DEFINED ENV{MKLROOT})

find_path(L_INCLUDE
  NAMES
  mkl.h
  PATHS
  ${MKLROOT}
  PATH_SUFFIXES
  include
)

#find pardiso libraries
set(MKLROOT_INTEL64 "${MKLROOT}/lib/intel64")

find_library(MKL_CORE_LIB mkl_core PATHS ${MKLROOT_INTEL64})
if(NOT MKL_CORE_LIB)
	message(FATAL_ERROR "library ?/libmkl_core.so not found")
endif (NOT MKL_CORE_LIB)

find_library(MKL_ILP64_LIB mkl_intel_ilp64 PATHS ${MKLROOT_INTEL64})
if(NOT MKL_ILP64_LIB)
	message(FATAL_ERROR "library ?/libmkl_intel_ilp64.so not found")
endif (NOT MKL_ILP64_LIB)

find_library(MKL_GNU_LIB mkl_gnu_thread PATHS ${MKLROOT_INTEL64})
if(NOT MKL_GNU_LIB)
	message(FATAL_ERROR "library ?/libmkl_gnu_thread.so not found")
endif (NOT MKL_GNU_LIB)

#add compile line option -DMKL_ILP64
add_compile_definitions(MKL_ILP64)

#message(STATUS "Val: ${MKLROOT_INTEL64}")
#message(STATUS "Val: ${MKL_CORE_LIB}")
#message(STATUS "Val: ${MKL_ILP64_LIB}")
#message(STATUS "Val: ${MKL_GNU_LIB}")
set(CMAKE_EXE_LINKER_FLAGS "-Wl,--no-as-needed")

mark_as_advanced(LAPACK_INCLUDE MKL_CORE_LIB MKL_GNU_LIB MKL_ILP64_LIB)