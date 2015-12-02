#############################################################################
## Makefile -- New Version of my Makefile that works on both linux
##              and mac os x
## Ryan Nichol <rjn@hep.ucl.ac.uk>
##############################################################################
include Makefile.arch


#If you have the 64-bit version of fftw installed, try this to help CINT out.
#CINTFLAGS=-DFFTW_64_BIT


#Site Specific  Flags (adjust to local site)
OUTLIBDIR          = 
OUTINCDIR          = 
SYSINCLUDES	= -I/usr/local/include
SYSLIBS         = -L/usr/local/lib
DLLSUF = ${DllSuf}
OBJSUF = ${ObjSuf}
SRCSUF = ${SrcSuf}

ifdef ANITA_UTIL_INSTALL_DIR
UTIL_LIB_DIR=${ANITA_UTIL_INSTALL_DIR}/lib
UTIL_INC_DIR=${ANITA_UTIL_INSTALL_DIR}/include
LD_UTIL=-L$(ANITA_UTIL_LIB_DIR)
INC_UTIL=-I$(ANITA_UTIL_INC_DIR)
else 
ifdef ARA_UTIL_INSTALL_DIR
UTIL_LIB_DIR=${ARA_UTIL_INSTALL_DIR}/lib
UTIL_INC_DIR=${ARA_UTIL_INSTALL_DIR}/include
LD_UTIL=-L$(ARA_UTIL_LIB_DIR)
INC_UTIL=-I$(ARA_UTIL_INC_DIR)
else
UTIL_LIB_DIR=/usr/local/lib
UTIL_INC_DIR=/usr/local/include
endif
endif

#Generic and Site Specific Flags
CXXFLAGS     += $(ROOTCFLAGS) $(SYSINCLUDES) 
LDFLAGS      += $(ROOTLDFLAGS) 
LIBS          = $(ROOTLIBS) -lMathMore $(SYSLIBS) -lfftw3
GLIBS         = $(ROOTGLIBS) $(SYSLIBS)


#Now the bits we're actually compiling


#ROOT stuff

ROOT_LIBRARY = libRootFftwWrapper.${DllSuf}
LIB_OBJS =  FFTWComplex.o FFTtools.o RFSignal.o RFFilter.o FFTtoolsWisdomManager.o fftDict.o
CLASS_HEADERS =   FFTWComplex.h FFTtools.h RFSignal.h FFTtoolsWisdomManager.h RFFilter.h
BINARIES = testFFTtools

all : $(ROOT_LIBRARY) $(BINARIES)

$(BINARIES): %: %.$(SRCSUF) $(ROOT_LIBRARY) 
	@echo "<**Compiling**> "
	@echo $<
	$(LD) $(CXXFLAGS) $(LDFLAGS) $(LIBS) $< $(ROOT_LIBRARY) -o $@

fftDict.C: $(CLASS_HEADERS)
	@echo "Generating dictionary ..."
	@ rm -f *Dict* 
	rootcint $@ -c -p -I$(shell $(RC) --incdir) $(SYSINCLUDES) $(CINTFLAGS) $(CLASS_HEADERS) LinkDef.h


#The library
$(ROOT_LIBRARY) : $(LIB_OBJS) 
	@echo "Linking $@ ..."
ifeq ($(PLATFORM),macosx)
# We need to make both the .dylib and the .so
		$(LD) $(SOFLAGS)$@ $(LDFLAGS) $^ $(LIBS) $(OutPutOpt) $@
ifneq ($(subst $(MACOSX_MINOR),,1234),1234)
ifeq ($(MACOSX_MINOR),4)
		ln -sf $@ $(subst .$(DllSuf),.so,$@)
else
		$(LD) -bundle -undefined $(UNDEFOPT) $(LDFLAGS) $(LIBS) $^ \
		   $(OutPutOpt) $(subst .$(DllSuf),.so,$@)
endif
endif
else
	$(LD) $(SOFLAGS) $(LDFLAGS) $(LIBS) $(LIB_OBJS) -o $@
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
	@rm -f $(subst .$(DllSuf),.so,$(ROOT_LIBRARY))	
	@rm -f $(TEST)


install: $(ROOT_LIBRARY)
	install -d $(UTIL_INC_DIR)
	install -d $(UTIL_LIB_DIR)
ifeq ($(PLATFORM),macosx)
	@install -c -m 755 $(ROOT_LIBRARY) $(subst .$(DllSuf),.so,$(ROOT_LIBRARY)) $(UTIL_LIB_DIR)
else
	install -c -m 755 $(ROOT_LIBRARY) $(UTIL_LIB_DIR)
endif
	install -c -m 644 $(CLASS_HEADERS) $(UTIL_INC_DIR)
