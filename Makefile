### Makefile that delegates building either to cmake or legacy Makefile 
### Cmake is default, if you don't want to use CMake, you can do make legacy   
### ( or move this file and rename Makefile.legacy to Makefile,
###  or modify the all/clean/install targets below  ) 
###


.PHONY: all configure clean cleaner install legacy legacy-clean legacy-install cmake-build cmake-clean cmake-install

all: cmake-build 
clean: cmake-clean
install: cmake-install 


### TODO add doxygen into CMakelists 
doc: legacy-doc

cmake-build: build/Makefile 
	@+make -C  ./build

legacy-doc: 
	@make -f Makefile.legacy doc 

legacy: 
	@+make -f Makefile.legacy 

legacy-clean: 
	@make -f Makefile.legacy clean 

legacy-install: 
	@make -f Makefile.legacy install 

configure: build/Makefile 
	@ccmake . build 

cmake-install: 
	@make -C ./build install 

build/Makefile: 
	@echo "Setting up cmake build. Use "make legacy" to use legacy Makefile" 
	@mkdir -p build 
	@cd build && cmake ../ 

distclean: 
	@echo "Removing cmake directory" 
	@rm -rf build 

cmake-clean: build/Makefile 
	@make -C ./build clean 
