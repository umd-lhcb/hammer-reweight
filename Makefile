.PHONY: dev-shell clean clean-nix clean-general patch build

BINPATH	:=	bin
VPATH	:=	utils:src:$(BINPATH)

export PATH := utils:$(BINPATH):$(PATH)

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
	@rm -rf ./gen/*

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


###############
# Sample plot #
###############

sample-plots: gen/el.png gen/q2.png gen/mm2.png

gen/el.png gen/q2.png gen/mm2.png &: \
	samples/rdst-run1.root \
	gen/rdst-run1-ff_w.root \
	plot_ratio.py
	$(word 3, $^) -d $< -w $(word 2, $^) -t mc_dst_tau -T mc_dst_tau_ff_w

gen/rdst-run1-ff_w.root: \
	samples/rdst-run1.root \
	rdx-run1-sample
	$(word 2, $^) $< $@


####################
# Generic patterns #
####################

%: %.cpp
	$(COMPILER) $(CXXFLAGS) -o $(BINPATH)/$@ $< $(LINKFLAGS) $(ADDLINKFLAGS)
