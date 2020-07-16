.PHONY: dev-shell clean-nix clean-general patch build

PWD=$(shell pwd)

dev-shell:
	@nix-shell --pure -E "with import <nixpkgs> { overlays = [(import ./nix/overlay)]; }; callPackage ./nix/overlay/hammer-phys {}"

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
