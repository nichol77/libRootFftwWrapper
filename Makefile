#############################################################################
## Makefile -- New Version of my Makefile that works on both linux
##              and mac os x
## Ryan Nichol <rjn@hep.ucl.ac.uk>
##############################################################################
include Makefile.arch

 ### all configuration moved to this file 
include Makefile.config


### Package subdirectories
LIBDIR=lib
BUILDDIR=build
INCLUDEDIR=include
BINDIR=bin
VECTORDIR=vectorclass






#Generic and Site Specific Flags
CXXFLAGS     += $(ROOTCFLAGS) $(SYSINCLUDES)
LDFLAGS      += $(ROOTLDFLAGS) -L$(LIBDIR) -L$(UTIL_LIB_DIR)
LIBS          = $(ROOTLIBS) -lMathMore -lMinuit2 $(SYSLIBS) -lfftw3 
GLIBS         = $(ROOTGLIBS) $(SYSLIBS)






#Now the bits we're actually compiling
.PHONY: all clean install doc


#ROOT stuff
ROOT_LIBRARY = $(LIBDIR)/libRootFftwWrapper.${DllSuf}

LIB_OBJS =  $(addprefix $(BUILDDIR)/, FFTWComplex.o FFTtools.o\
																			RFSignal.o RFFilter.o \
						                          FFTWindow.o SineSubtract.o \
																			DigitalFilter.o RFInterpolate.o\
																		 	AnalyticSignal.o Averager.o \
																			CWT.o fftDict.o) 

CLASS_HEADERS =   $(addprefix $(INCLUDEDIR)/, FFTWComplex.h FFTtools.h \
																							RFSignal.h RFFilter.h\
																							FFTWindow.h SineSubtract.h \
																							RFInterpolate.h DigitalFilter.h\
																							Averager.h AnalyticSignal.h CWT.h) 

BINARIES = $(addprefix $(BINDIR)/, testFFTtools)


all : $(ROOT_LIBRARY) $(BINARIES)


#force recompile if any part of Makefile updated 
Makefile: Makefile.config Makefile.arch 
	touch Makefile


$(BINDIR)/%: test/%.$(SRCSUF) $(ROOT_LIBRARY) Makefile | $(BINDIR) 
	@echo "<**Compiling**> "
	@echo $<
	$(LD) -I$(INCLUDEDIR)  $(CXXFLAGS) $(LDFLAGS) $(LIBS) -lRootFftwWrapper $< -o $@

$(LIB_OBJS): | $(BUILDDIR) 

$(BINDIR): 
	mkdir -p $(BINDIR)

$(BUILDDIR): 
	mkdir -p $(BUILDDIR)


$(LIBDIR): 
	mkdir -p $(LIBDIR)

#### Download and unzip Agner Fog's VCL vectorization class  
$(VECTORDIR): 
	mkdir -p $(VECTORDIR) 
	curl http://www.agner.org/optimize/vectorclass.zip > $(VECTORDIR)/vectorclass.zip 
	unzip $(VECTORDIR)/vectorclass.zip -d $(VECTORDIR) 


$(BUILDDIR)/fftDict.C: $(CLASS_HEADERS) LinkDef.h | $(BUILDDIR) 
	@echo "Generating dictionary ..."
	rootcint -f $@ -c -p -I$(shell $(RC) --incdir) $(SYSINCLUDES) $(CINTFLAGS) $(CLASS_HEADERS) LinkDef.h


#The library
$(ROOT_LIBRARY) : $(LIB_OBJS) | $(LIBDIR) 
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

$(BUILDDIR)/%.$(OBJSUF) : src/%.$(SRCSUF) $(CLASS_HEADERS) Makefile | $(BUILDDIR) $(VECTORIZE) 
	@echo "<**Compiling**> "$<
	$(CXX) -I$(INCLUDEDIR) $(CXXFLAGS)  -c $< -o  $@

$(BUILDDIR)/%.$(OBJSUF) : $(BUILDDIR)/%.C
	@echo "<**Compiling**> "$<
	$(CXX) -I$(INCLUDEDIR) -I./ $(CXXFLAGS) -c $< -o  $@


clean:
	rm -rf $(BUILDDIR) 
	rm -rf $(BINDIR) 
	rm -rf $(LIBDIR) 
	rm -rf doc/html
	rm -rf doc/latex

doc: 
	doxygen doc/doxynew.config

install: $(ROOT_LIBRARY)
	install -d $(UTIL_INC_DIR)
	install -d $(UTIL_LIB_DIR)
ifeq ($(PLATFORM),macosx)
	@install -c -m 755 $(ROOT_LIBRARY) $(subst .$(DllSuf),.so,$(ROOT_LIBRARY)) $(UTIL_LIB_DIR)
else
	install -c -m 755 $(ROOT_LIBRARY) $(UTIL_LIB_DIR)
endif
	install -c -m 644 $(CLASS_HEADERS) $(UTIL_INC_DIR)

# Try installing pcm file with library.
# This is a new thing for ROOT6 and at the present time doesn't seem to matter if this fails.
# So we will suppress the warning and continue if it fails.

	if [ -e $(BUILDDIR)/fftDict_rdict.pcm ];  then install -c -m 755 $(BUILDDIR)/fftDict_rdict.pcm $(UTIL_LIB_DIR); fi; 
