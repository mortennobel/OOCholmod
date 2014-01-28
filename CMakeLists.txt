cmake_minimum_required (VERSION 2.6)
project (oocholmod)

ADD_DEFINITIONS(
    -std=c++11 # Or -std=c++0x
    # Other flags
)

# should we use our own math functions?
option (USE_ACML 
        "Use AMD Core Math Library (ACML)" ON) 

if (USE_ACML)
	set(ACML_LIBRARY "" CACHE STRING "ACML library")
	set(ACML_INCLUDE "" CACHE STRING "ACML header file location")
else (USE_ACML)	
	IF(APPLE)        
		FIND_LIBRARY(ACCELERATE_LIBRARY Accelerate)
		MARK_AS_ADVANCED (ACCELERATE_LIBRARY)
		SET(EXTRA_LIBS ${ACCELERATE_LIBRARY})
	endif (APPLE)
	set(BLAS_LIBRARY "" CACHE STRING "Blas library")	
	set(BLAS_INCLUDE "" CACHE STRING "Blas header files location")	
	set(LAPACK_LIBRARY "" CACHE STRING "Blas library")	
	set(LAPACK_INCLUDE "" CACHE STRING "Blas header files location")
	
 endif (USE_ACML)		

if (USE_ACML)
	if (ACML_INCLUDE)
		include_directories ("${ACML_INCLUDE}")
	endif ()
else (USE_ACML)
	if (BLAS_INCLUDE)
		include_directories ("${BLAS_INCLUDE}")
	endif ()
	if (LAPACK_INCLUDE)
		include_directories ("${LAPACK_INCLUDE}")
	endif ()
endif (USE_ACML)
 
set(CHOLMOD_INCLUDE "" CACHE STRING "Location of Cholmod include")
set(CHOLMOD_LIBRARY "" CACHE STRING "Location of Cholmod library")
 
 if (CHOLMOD_INCLUDE)
	include_directories ("${CHOLMOD_INCLUDE}")
 endif ()
 
# include header of scr
include_directories ("${PROJECT_SOURCE_DIR}/src")
# add subdirectory (cmakelists.txt)
add_subdirectory (src) 

# include headers of test
include_directories ("${PROJECT_SOURCE_DIR}/test")
add_subdirectory (test)



