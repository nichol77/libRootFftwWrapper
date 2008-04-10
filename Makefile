#############################################################################
## Makefile -- New Version of my Makefile that works on both linux
##              and mac os x
## Ryan Nichol <rjn@hep.ucl.ac.uk>
##############################################################################
include Makefile.arch

#Site Specific  Flags (adjust to local site)
SYSINCLUDES	= -I/sw/include
SYSLIBS         = -L/sw/lib

ifdef ANITA_UTIL_INSTALL_DIR
ANITA_UTIL_LIB_DIR=${ANITA_UTIL_INSTALL_DIR}/lib
ANITA_UTIL_INC_DIR=${ANITA_UTIL_INSTALL_DIR}/includes
else
ANITA_UTIL_LIB_DIR=
ANITA_UTIL_INC_DIR=
endif

#Generic and Site Specific Flags
CXXFLAGS     += $(ROOTCFLAGS) $(SYSINCLUDES) 
LDFLAGS      += $(ROOTLDFLAGS) 
LIBS          = $(ROOTLIBS) -lMathMore $(SYSLIBS)
GLIBS         = $(ROOTGLIBS) $(SYSLIBS)


#Now the bits we're actually compiling


#ROOT stuff

ROOT_LIBRARY = libRootFftwWrapper.${DLLSUF}
LIB_OBJS =  FFTWComplex.o FFTtools.o fftDict.o
CLASS_HEADERS =   FFTWComplex.h FFTtools.h

all : $(ROOT_LIBRARY) 

fftDict.C: $(CLASS_HEADERS)
	@echo "Generating dictionary ..."
	@ rm -f *Dict* 
	rootcint $@ -c $(CXXFLAGS) $(CLASS_HEADERS) LinkDef.h


#The library
$(ROOT_LIBRARY) : $(LIB_OBJS) 
	@echo "Linking $@ ..."
ifeq ($(PLATFORM),macosx)
# We need to make both the .dylib and the .so
	$(LD) $(SOFLAGS) $^ $(OutPutOpt) $@
ifeq ($(MACOSX_MINOR),4)
	ln -sf $@ $(subst .$(DLLSUF),.so,$@)
else
	$(LD) -bundle -undefined $(UNDEFOPT) $(LDFLAGS) $^ \
	 $(OutPutOpt) $(subst .$(DLLSUF),.so,$@)
endif
else
	$(LD) $(SOFLAGS) $(LDFLAGS) $(LIB_OBJS) -o $@
endif

%.$(OBJSUF) : %.$(SRCSUF)
	@echo "<**Compiling**> "$<
	$(CXX) $(CXXFLAGS) -c $< -o  $@

%.$(OBJSUF) : %.C
	@echo "<**Compiling**> "$<
	$(CXX) $(CXXFLAGS) $ -c $< -o  $@


clean:
	@rm -f *Dict*
	@rm -f *.${OBJSUF}
	@rm -f $(LIBRARY)
	@rm -f $(ROOT_LIBRARY)
	@rm -f $(subst .$(DLLSUF),.so,$(ROOT_LIBRARY))	
	@rm -f $(TEST)


install: $(ROOT_LIBRARY)
ifeq ($(PLATFORM),macosx)
	@cp $(ROOT_LIBRARY) $(subst .$(DLLSUF),.so,$(ROOT_LIBRARY)) $(ANITA_UTIL_LIB_DIR)
else
	@cp $(ROOT_LIBRARY) $(ANITA_UTIL_LIB_DIR)
endif
	@cp  $(CLASS_HEADERS) $(ANITA_UTIL_INC_DIR)
