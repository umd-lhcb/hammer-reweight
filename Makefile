.PHONY: dev-shell clean clean-nix clean-general patch build

BINPATH	:=	bin
VPATH	:=	src

# System env
PWD=$(shell pwd)

# Compiler settings
COMPILER	:=	$(shell root-config --cxx)
#CXXFLAGS	:=	$(shell root-config --cflags)
CXXFLAGS	:=	-pthread -std=c++14 -m64
LINKFLAGS	:=	$(shell root-config --libs)
ADDLINKFLAGS	:=	-lHammerTools -lHammerBase -lHammerCore -lFormFactors -lAmplitudes -lRates

dev-shell:
	@nix-shell --pure -E "with import <nixpkgs> { overlays = [(import ./nix/overlay)]; }; callPackage ./nix/overlay/hammer-phys {}"

clean:
	@rm -rf ./bin/*

clean-nix:
	@sudo rm -rf ./Hammer-*
	@sudo rm -rf ./out

clean-general:
	@rm -rf ./out ./Hammer/build

patch:
	@cd Hammer && \
	git reset --hard && \
	git apply ../nix/overlay/hammer-phys/add_missing_header.patch

build: patch
	@mkdir -p out
	@cd Hammer && \
	cmake -S . -B build \
		-DCMAKE_INSTALL_PREFIX=../out \
		-DWITH_PYTHON=ON -DWITH_ROOT=ON \
		-DWITH_EXAMPLES=ON \
		-DINSTALL_EXTERNAL_DEPENDENCIES=ON \
		-DMAX_CXX_IS_14=OFF && \
	cd build && \
	make && \
	make install
	@echo ""
	@echo "Please add the following lines to your shell config:"
	@echo "export=CPLUS_INCLUDE_PATH=$(PWD)/out/include:"'$$CPLUS_INCLUDE_PATH'
	@echo "export=LD_LIBRARY_PATH=$(PWD)/out/lib:"'$$LD_LIBRARY_PATH'


####################
# Generic patterns #
####################

%: %.cpp
	$(COMPILER) $(CXXFLAGS) -o $(BINPATH)/$@ $< $(LINKFLAGS) $(ADDLINKFLAGS)
