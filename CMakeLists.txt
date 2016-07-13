# CMakeLists.txt for event package. It creates a library with dictionary and a main program
cmake_minimum_required(VERSION 2.8.10)

project(RootFftwWrapper)
set(libname "RootFftwWrapper")

set (LIBRARY_VERSION_MAJOR 3)
set (LIBRARY_VERSION_MINOR 0)
set (LIBRARY_VERSION_PATCH 1)
set (LIBRARY_VERSION_STRING ${LIBRARY_VERSION_MAJOR}.${LIBRARY_VERSION_MINOR}.${LIBRARY_VERSION_PATCH})
 
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/modules/")

# You need to tell CMake where to find the ROOT installation. This can be done in a number of ways:
#   - ROOT built with classic configure/make use the provided $ROOTSYS/etc/cmake/FindROOT.cmake
#   - ROOT built with CMake. Add in CMAKE_PREFIX_PATH the installation prefix for ROOT
list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS})


if(DEFINED ENV{ANITA_UTIL_INSTALL_DIR})
    message("ANITA_UTIL_INSTALL_DIR is set to $ENV{ANITA_UTIL_INSTALL_DIR}")
    set(UTIL_LIB_DIR $ENV{ANITA_UTIL_INSTALL_DIR}/lib)
    set(UTIL_INC_DIR $ENV{ANITA_UTIL_INSTALL_DIR}/include)
    set(UTIL_BIN_DIR $ENV{ANITA_UTIL_INSTALL_DIR}/bin)
    set(UTIL_SHARE_DIR $ENV{ANITA_UTIL_INSTALL_DIR}/share)
    set(LD_UTIL $ENV{ANITA_UTIL_INSTALL_DIR}/lib)
    set(INC_UTIL $ENV{ANITA_UTIL_INSTALL_DIR}/include)
elseif(DEFINED ENV{ARA_UTIL_INSTALL_DIR})
    message("ARA_UTIL_INSTALL_DIR is set to $ENV{ANITA_UTIL_INSTALL_DIR}")
    set(UTIL_LIB_DIR $ENV{ARA_UTIL_INSTALL_DIR}/lib)
    set(UTIL_INC_DIR $ENV{ARA_UTIL_INSTALL_DIR}/include)
    set(UTIL_BIN_DIR $ENV{ARA_UTIL_INSTALL_DIR}/bin)
    set(UTIL_SHARE_DIR $ENV{ARA_UTIL_INSTALL_DIR}/share)
    set(LD_UTIL $ENV{ARA_UTIL_INSTALL_DIR}/lib)
    set(INC_UTIL $ENV{ARA_UTIL_INSTALL_DIR}/include)
else()
    message("Will install to /usr/local  ANITA_UTIL_INSTALL_DIR not set")
    set(UTIL_LIB_DIR /usr/local/lib)
    set(UTIL_INC_DIR /usr/local/include)
    set(UTIL_BIN_DIR /usr/local/bin)
    set(UTIL_SHARE_DIR /usr/local/share)
    set(LD_UTIL /usr/local/lib)
    set(INC_UTIL /usr/local/include)
endif()


#---Locate the ROOT package and defines a number of variables (e.g. ROOT_INCLUDE_DIRS)
find_package(ROOT REQUIRED COMPONENTS MathMore MathCore RIO Hist Tree Net Minuit2)
find_package(FFTW REQUIRED)

#---Define useful ROOT functions and macros (e.g. ROOT_GENERATE_DICTIONARY)

message("ROOT_VERSION is set to ${ROOT_VERSION}")
if( ${ROOT_VERSION} VERSION_GREATER "5.99/99")
    message("Using ROOT_VERSION 6")
    include(${ROOT_USE_FILE})
else()
#    include(RootNewMacros)
     add_definitions(${ROOT_DEFINITIONS})	
endif()


message("ROOT_INCLUDE_DIRS is set to ${ROOT_INCLUDE_DIRS}")
message("FFTW_INCLUDES is set to ${FFTW_INCLUDES}")

include_directories(${PROJECT_SOURCE_DIR} ${ROOT_INCLUDE_DIRS} ${PROJECT_SOURCE_DIR}/include ${FFTW_INCLUDES})
add_definitions(${ROOT_CXX_FLAGS})

file(GLOB HEADER_FILES
    "include/*.h"    
)

file(GLOB_RECURSE SOURCE_FILES src/*.cxx)
#file(GLOB_RECURSE THIS_LD LinkDef.h)


set(DICTNAME G__${libname})

ROOT_GENERATE_DICTIONARY(${DICTNAME} ${HEADER_FILES} LINKDEF LinkDef.h)

#---Create a shared library with geneated dictionary
add_library(${libname} SHARED ${SOURCE_FILES} ${DICTNAME}.cxx)
set_target_properties(${libname} PROPERTIES VERSION ${LIBRARY_VERSION_STRING}
  SOVERSION ${LIBRARY_VERSION_MAJOR})


#target_link_libraries(${libname} ${ROOT_LIBRARIES} ${FFTW_LIBRARIES} MathMore Minuit2)
target_link_libraries(${libname} ${ROOT_LIBRARIES} ${FFTW_LIBRARIES} )


if( ${ROOT_VERSION} VERSION_GREATER "5.99.99")
    add_custom_target(${DICTNAME}.pcm DEPENDS ${DICTNAME})
endif()


message("UTIL_INC_DIR is set to ${UTIL_INC_DIR}")
message("DICTNAME is set to ${DICTNAME}" )
message("PROJECT_BINARY_DIR is set to ${PROJECT_BINARY_DIR}")

install (FILES ${HEADER_FILES} DESTINATION ${UTIL_INC_DIR})
install (FILES cmake/modules/FindFFTW.cmake DESTINATION ${UTIL_SHARE_DIR}/cmake/modules)
install (TARGETS ${libname}
	ARCHIVE DESTINATION ${UTIL_LIB_DIR}
	LIBRARY DESTINATION ${UTIL_LIB_DIR}
	RUNTIME DESTINATION ${UTIL_BIN_DIR})
#Only needed for ROOT6
if( ${ROOT_VERSION} VERSION_GREATER "5.99.99")
  install (FILES ${PROJECT_BINARY_DIR}/${DICTNAME}_rdict.pcm DESTINATION ${UTIL_LIB_DIR})
#  install (FILES ${PROJECT_BINARY_DIR}/lib${libname}.rootmap DESTINATION ${UTIL_LIB_DIR})
endif() 
