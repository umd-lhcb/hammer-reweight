{
  description = "Code for FF reweighting in  HAMMER.";

  inputs = {
    root-curated.url = "github:umd-lhcb/root-curated";
    nixpkgs.follows = "root-curated/nixpkgs";
    flake-utils.follows = "root-curated/flake-utils";
  };

  outputs = { self, nixpkgs, flake-utils, root-curated }:
    {
      overlay = import ./nix/overlay.nix;
    } //
    flake-utils.lib.eachDefaultSystem (system:
      let
        pkgs = import nixpkgs {
          inherit system;
          overlays = [ root-curated.overlay self.overlay ];
        };
      in
      rec {
        packages = flake-utils.lib.flattenTree {
          dev-shell = devShell.inputDerivation;
        };
        devShell = pkgs.mkShell {
          name = "hammer-reweight";
          buildInputs = with pkgs; [
            root
            hammer-phys
            ff_calc
            python3
            clang-tools
          ];

          shellHook = ''
            export PATH=$(pwd)/bin:$(pwd)/utils:$PATH
          '';
        };
      });
}
