# CMakeLists.txt for event package. It creates a library with dictionary and a main program
cmake_minimum_required(VERSION 3.0 FATAL_ERROR)
project(RootFftwWrapper)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/")

# You need to tell CMake where to find the ROOT installation. This can be done in a number of ways:
#   - ROOT built with classic configure/make use the provided $ROOTSYS/etc/cmake/FindROOT.cmake
#   - ROOT built with CMake. Add in CMAKE_PREFIX_PATH the installation prefix for ROOT
list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS})


if(DEFINED ENV{ANITA_UTIL_INSTALL_DIR})
    message("ANITA_UTIL_INSTALL_DIR is set to $ENV{ANITA_UTIL_INSTALL_DIR}")
    set(UTIL_LIB_DIR $ENV{ANITA_UTIL_INSTALL_DIR}/lib)
    set(UTIL_INC_DIR $ENV{ANITA_UTIL_INSTALL_DIR}/include)
    set(UTIL_BIN_DIR $ENV{ANITA_UTIL_INSTALL_DIR}/bin)
    set(LD_UTIL $ENV{ANITA_UTIL_INSTALL_DIR}/lib)
    set(INC_UTIL $ENV{ANITA_UTIL_INSTALL_DIR}/include)
elseif(DEFINED ENV{ARA_UTIL_INSTALL_DIR})
    message("ARA_UTIL_INSTALL_DIR is set to $ENV{ANITA_UTIL_INSTALL_DIR}")
    set(UTIL_LIB_DIR $ENV{ARA_UTIL_INSTALL_DIR}/lib)
    set(UTIL_INC_DIR $ENV{ARA_UTIL_INSTALL_DIR}/include)
    set(UTIL_BIN_DIR $ENV{ARA_UTIL_INSTALL_DIR}/bin)	
    set(LD_UTIL $ENV{ARA_UTIL_INSTALL_DIR}/lib)
    set(INC_UTIL $ENV{ARA_UTIL_INSTALL_DIR}/include)
else()
    message("Will install to /usr/local  ANITA_UTIL_INSTALL_DIR not set")
    set(UTIL_LIB_DIR /usr/local/lib)
    set(UTIL_INC_DIR /usr/local/include)
    set(UTIL_BIN_DIR /usr/local/bin)	
    set(LD_UTIL /usr/local/lib)
    set(INC_UTIL /usr/local/include)
endif()


#---Locate the ROOT package and defines a number of variables (e.g. ROOT_INCLUDE_DIRS)
find_package(ROOT REQUIRED COMPONENTS MathMore MathCore RIO Hist Tree Net)
find_package(FFTW REQUIRED)

#---Define useful ROOT functions and macros (e.g. ROOT_GENERATE_DICTIONARY)
include(${ROOT_USE_FILE})


include_directories(${CMAKE_SOURCE_DIR} ${ROOT_INCLUDE_DIRS} include ${FFTW_INCLUDES})
add_definitions(${ROOT_CXX_FLAGS})

file(GLOB OUR_HEADER
    "include/*.h"    
)

file(GLOB OUR_SOURCE
    "src/*.cxx"    
)

file(GLOB_RECURSE SOURCE_FILES src/*.cxx)

set(DICTNAME "G__FFtDict")

ROOT_GENERATE_DICTIONARY(${DICTNAME} ${OUR_HEADER} LINKDEF LinkDef.h)

#---Create a shared library with geneated dictionary
add_library(RootFftwWrapper SHARED ${SOURCE_FILES} ${DICTNAME}.cxx)
target_link_libraries(RootFftwWrapper ${ROOT_LIBRARIES} ${FFTW_LIBRARIES} MathMore Minuit2)

add_custom_target(${DICTNAME}.pcm DEPENDS ${DICTNAME})


message("UTIL_INC_DIR is set to ${UTIL_INC_DIR}")	   
message("OUR_SOURCE is set to ${OUR_SOURCE}")
message("DICTNAME is set to ${DICTNAME}" )
message("CMAKE_LIBRARY_OUTPUT_DIRECTORY is set to ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}")
message("PROJECT_BINARY_DIR is set to ${PROJECT_BINARY_DIR}")

install (FILES ${OUR_HEADER} DESTINATION ${UTIL_INC_DIR})
install (FILES ${PROJECT_BINARY_DIR}/${DICTNAME}_rdict.pcm DESTINATION ${UTIL_LIB_DIR})
install (TARGETS RootFftwWrapper
	ARCHIVE DESTINATION ${UTIL_LIB_DIR}
	LIBRARY DESTINATION ${UTIL_LIB_DIR}
	RUNTIME DESTINATION ${UTIL_BIN_DIR})
