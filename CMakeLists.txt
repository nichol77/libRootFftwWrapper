# CMakeLists.txt for libRootFftwWrapper.
# 
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
find_package(ROOT REQUIRED COMPONENTS MathMore MathCore RIO Hist Tree Net
  Minuit  Minuit2 Spectrum )
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


file (GLOB CMAKE_FILES "cmake/modules/*.cmake")

install (FILES ${HEADER_FILES} DESTINATION ${UTIL_INC_DIR})
install (FILES ${CMAKE_FILES} DESTINATION ${UTIL_SHARE_DIR}/cmake/modules)
install (TARGETS ${libname}
	ARCHIVE DESTINATION ${UTIL_LIB_DIR}
	LIBRARY DESTINATION ${UTIL_LIB_DIR}
	RUNTIME DESTINATION ${UTIL_BIN_DIR})
#Only needed for ROOT6
if( ${ROOT_VERSION} VERSION_GREATER "5.99.99")
  install (FILES ${PROJECT_BINARY_DIR}/${DICTNAME}_rdict.pcm DESTINATION ${UTIL_LIB_DIR})
#  install (FILES ${PROJECT_BINARY_DIR}/lib${libname}.rootmap DESTINATION ${UTIL_LIB_DIR})
endif() 



###############################################################################################
### Binaries
### I'm a CMAKE n00b so ther's probably a much better way to do this? 
###############################################################################################

set (BINSRCDIR test) 
macro( do_binary binary_name ) 
  add_executable(${binary_name} ${BINSRCDIR}/${binary_name}.cxx) 
  target_link_libraries(${binary_name}  ${ROOT_LIBRARIES} ${FFTW_LIBRARIES} ${libname})
endmacro() 

do_binary(testFFTtools) 
do_binary(testSubtract) 



#################################################################################3
### The default build configuration is INSANE. No optimization? What is this 1971? 

message (" Adding new build type") 

set(CMAKE_CXX_FLAGS_DEFAULT 
  "-Os -g -Wall" 
  CACHE STRING "c++ Flags used during default libRootFftwWrapper builds" 
  FORCE ) 

set(CMAKE_C_FLAGS_DEFAULT 
  "-Os -g -Wall"
  CACHE STRING "c Flags used during default libRootFftwWrapper builds" 
  FORCE ) 

set(CMAKE_EXE_LINKER_FLAGS_DEFAULT 
  "-g"
  CACHE STRING "ld Flags used during default libRootFftwWrapper builds" 
  FORCE ) 

set(CMAKE_SHARED_LINKER_FLAGS_DEFAULT 
  "-g"
  CACHE STRING "ld Flags used during default libRootFftwWrapper builds" 
  FORCE ) 


mark_as_advanced ( CMAKE_CXX_FLAGS_DEFAULT  CMAKE_C_FLAGS_DEFAULT CMAKE_EXE_LINKER_FLAGS_DEFAULT CMAKE_SHARED_LINKER_FLAGS_DEFAULT) 

if (NOT CMAKE_BUILD_TYPE) 
  set (CMAKE_BUILD_TYPE Default
    CACHE STRING "Choose tye type of build: None Debug Release RelWithDebInfo MinSizeRel Default"
    FORCE ) 
endif()



#################################################################################
###  Manual vectorization option 
###  Maybe in the future this can be made a CMake module or something 
################################################################################# 

option(VECTORIZE "Enable Manual SIMD Vectorization" ON) 
set (VECTORDIR vectorclass) 
if (VECTORIZE) #defined in anitaBuildTool, default on
  add_definitions( -DENABLE_VECTORIZE ) 

  if(CMAKE_COMPILER_IS_GNUCXX) 

    ### someone should do this for clang if they want it to be as fast as possible 
    if (CMAKE_CXX_COMPILER_VERSION VERSION_GREATER 4.1)
      if(NATIVE_ARCH) #define in anitaBuildTool, default on
	add_definitions( -march=native)
      endif()
    endif()
    add_definitions(  -fabi-version=0 ) 
  endif() 

  if (EXISTS ${VECTORDIR}/vectorclass.zip )
    message (" Found vectorclass already downloaded... ") 

  else () 
    file(MAKE_DIRECTORY ${VECTORDIR} ) 
    message (" Downloading vector class...")
    file(DOWNLOAD http://www.agner.org/optimize/vectorclass.zip ${VECTORDIR}/vectorclass.zip) 

    #TODO there is probably a more cross platform way of doing this... 
    execute_process(COMMAND unzip ${VECTORDIR}/vectorclass.zip -d ${VECTORDIR} )
    execute_process(COMMAND unzip ${VECTORDIR}/special.zip -d ${VECTORDIR} ) 
    
    #### TERRIBLE TERRIBLE TERRIBLE HACK FOR MISBEHAVING CLANG VERSION REPORTING 
    FILE(READ ${VECTORDIR}/vectorf128.h INFILE ) 
    STRING(REGEX REPLACE "6020000" "1000000000" MINFILE ${INFILE} )
    FILE (WRITE ${VECTORDIR}/vectorf128.h ${MINFILE} )

  endif() 


  ## is this right? who knows 
  include_directories( ${CMAKE_BINARY_DIR}/${VECTORDIR} ) 


endif (VECTORIZE) 



##################################################################
###  Linear algebra options 
##################################################################

option(USE_EIGEN "Use Eigen3 for certain linear algebra options" OFF) 
set (EIGEN3_INCLUDE_DIR "/usr/include/eigen3" CACHE String "Eigen3 include path") 
option(USE_ARMADILLO "Use Armadillo for certain linear algebra options" OFF) 

if (USE_EIGEN AND USE_ARMADILLO) 

  message(FATAL_ERROR "USE_EIGEN and USE_ARMADILLO are mutually exclusive") 

elseif(USE_EIGEN) 
  add_definitions ( -DUSE_EIGEN) 
  include_directories (  ${EIGEN3_INCLUDE_DIR} ) 
elseif(USE_ARMADILLO) 

  add_definitions( -DUSE_ARMADILLO) 
  target_link_libraries(${libname} armadillo) 

endif() 


###############################################################################
#### Multithread options 
################################################################################
option(FFTTOOLS_ENABLE_OPENMP "Enable OpenMP support (experimenta, mutually exclusive with FFTTOOLS_ENABLE_THREAD_SAFE) for doFft and doInvFft" OFF) 
option(FFTTOOLS_ENABLE_THREAD_SAFE "Make fft methods threadsafe (mutually exclusive with FFTTOOLS_ENABLE_OPEMP)" OFF) 

if (FFTTOOLS_ENABLE_OPENMP AND FFTTOOLS_ENABLE_THREAD_SAFE)

  message(FATAL_ERROR "FFTTOOLS_ENABLE_THREAD_SAFE and FFTTOOLS_ENABLE_OPENMP are mutually exclusive") 

elseif(FFTTOOLS_ENABLE_THREAD_SAFE) 
  add_definitions( -DFFTTOOLS_THREAD_SAFE )
elseif(FFTTOOLS_ENABLE_OPENMP) 

  FIND_PACKAGE( OpenMP REQUIRED)
  if(OPENMP_FOUND)
    message("OPENMP FOUND")
    add_definitions (-DFFTTOOLS_USE_OMP) 
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")

    do_binary(testOpenMP) 
  else() 
    message(FATAL_ERROR "you tried to use openmp but Cmake couldn't find compiler support") 
  endif() 
endif() 


#### api compatibility 

option ( FORCE_OLD_GPP_ABI " Force old g++ ABI; this might be necessary if using new g++ with ROOT compiled with older g++ or other similar situations" OFF ) 
if (FORCE_OLD_GPP_ABI) 
  add_definitions( -D_GLIBCXX_USE_CXX11_ABI=0 ) 
endif() 
    


### Let's see if this works 
macro(stupid_option option_name option_description) 
  option (${option_name} ${option_description} OFF) 
  if (${option_name}) 
    add_definitions ( -D${option_name} )
  endif()
endmacro() 




#### Random Sine Subtraction options #### 
stupid_option(SINE_SUBTRACT_USE_FLOATS  "Use floats for vectorized sine subtraction") 
mark_as_advanced(SINE_SUBTRACT_USE_FLOATS)
stupid_option(SINE_SUBTRACT_PROFILE  "Enable Sine Subtraction profiling (you probably don't want to do this)") 
mark_as_advanced(SINE_SUBTRACT_PROFILE)

#option ( SINE_SUBTRACT_USE_FLOATS "Use floats for vectorized Sine Subtraction" OFF)
#mark_as_advanced( SINE_SUBTRACT_USE_FLOATS) 
#if ( SINE_SUBTRACT_USE_FLOATS) 
#  add_definitions ( -DSINE_SUBTRACT_USE_FLOATS)
#endif()
#
#option ( SINE_SUBTRACT_PROFILE "Enable Sine Subtraction profiling (you probably don't want to do this) " OFF)
#mark_as_advanced( SINE_SUBTRACT_PROFILE) 
#if ( SINE_SUBTRACT_PROFILE) 
#  add_definitions ( -DSINE_SUBTRACT_PROFILE)
#endif()
#
#

#### Misc options  #### 

stupid_option(FFTW_USE_PATIENT  "Use FFTW Patient plans... not recommended unless you are good about saving wisdom")
mark_as_advanced(FFTW_USE_PATIENT) 

#option(USE_PATIENT_PLANS "Use FFTW_PATIENT" OFF) 
#mark_as_advanced( USE_PATIENT_PLANS) 
#
#if (USE_PATIENT_PLANS) 
#  add_definitions ( -DFFTW_USE_PATIENT ) 
#endif() 







