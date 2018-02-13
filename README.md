# libRootFftwWrapper 

libRootFftwWrapper -- A ROOT wrapper for FFTW3
Ryan Nichol (rjn@hep.ucl.ac.uk)  -- December 2007
###############################################################################

Introduction
===============
This is a simple wrapper library that provides an interface to FFTW for use with ROOT. 

Prerequisities
===============
ROOT (http://root.cern.ch)
FFTW3 (http://www.fftw.org/)


Making  and Installing
======================
i) Ensure both ROOT and FFTW3 are on your system and that ROOTSYS is correctly set (or root-config is in your path) and that the directory containing libfftw3 is in either the system or your own LD_LIBRARY_PATH.
ii) Make sure that cmake is installed on your system
iii) Set the environmental variable ANITA_UTIL_INSTALL_DIR to point to your desired installation dir (libraries will be installed in a lib subdirectory from this location)
iii) To build the code: 


  Right now, we recommend you build using cmake: 

   (optional) make configure
   make 
   make install 

  To use the old build system:  
   (optional) modify Makefile.config to suit your needs 
    make legacy 
    make legacy-install 

  Trying to mix the two build systems might get things confused 




Mac Users
=========
For Mac OS X I would recommned usinge homebrew to install FFTW, ROOT and cmake. See http://brew.sh After installing brew you can get all the prerequisites by doing
brew update
brew install cmake
brew install fftw
brew tap homebrew/science
brew install root6



Usage from within ROOT
=======================
a) Ensure that the directory containing the libRootFftwWrapper library is in your LD_LIBRARY_PATH (e.g. ANITA_UTIL_INSTALL_DIR/lib )
b) To load the library you can (note the lack of .so or .dylib):
gSystem->Load("libRootFftwWrapper");
c) The FFTWComplex and FFTtools classes should now be accessible from the ROOT prompt or ROOT macro.

FFTtools Description
======================
Most of the member functions have fairly descriptive names, the documentation is at present somewhat lacking.


Compatibility Mode (cmake only) 
===========================

At some point, FFTtools changed from a static-method-only class to a namespace
so that other things could be added to the FFTtools namespace.  This was
thought not to adversely affect anyone, since the FFTtools::doFFT, etc. usage
would remain the same. However, apparently some users had code which did things
like: FFTtools tools; tools.doFFT( );. Even though the instance is extraneous,
it may be a pain to go through a lot of legacy code to fix it. For this reason, 
a libRootFftwWrapperCompat shared library is also created, which should be compatible with 
old uses of the library. To compile old code against it, you will need to define FFTTOOLS_COMPAT_MODE somehow, either like:

#define FFTTOOLS_COMPAT_MODE
#include "FFTtools.h" 

or by adding -DFFTTOOLS_COMPAT_MODE to your compilation flags (either with CXX_FLAGS or in the build system) 


