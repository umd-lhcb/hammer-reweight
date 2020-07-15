dev-shell:
	@nix-shell --pure -E "with import <nixpkgs> { overlays = [(import ./nix/overlay)]; }; callPackage ./nix/overlay/hammer-phys {}"

clean:
	@sudo rm -rf ./Hammer-*
	@sudo rm -rf ./out
