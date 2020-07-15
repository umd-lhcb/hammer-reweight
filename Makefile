dev-shell:
	@nix-shell --pure -E "with import <nixpkgs> { overlays = [(import ./nix/overlay)]; }; callPackage ./nix/overlay/hammer-phys {}"

cleanup:
	@sudo rm -rf ./Hammer-*
