{
  description = "Code for FF reweighting in  HAMMER.";

  inputs = {
    root-curated.url = "github:umd-lhcb/root-curated";
    nixpkgs.follows = "root-curated/nixpkgs";
    flake-utils.follows = "root-curated/flake-utils";

    pyTuplingUtils.url = "github:umd-lhcb/pyTuplingUtils";
  };

  outputs = { self, nixpkgs, flake-utils, root-curated, pyTuplingUtils }:
    {
      overlay = import ./nix/overlay.nix;
    } //
    flake-utils.lib.eachDefaultSystem (system:
      let
        pkgs = import nixpkgs {
          inherit system;
          config = { allowUnfree = true; };
          overlays = [ root-curated.overlay self.overlay pyTuplingUtils.overlay ];
        };
        python = pkgs.python3;
        pythonPackages = python.pkgs;
      in
      rec {
        packages = flake-utils.lib.flattenTree {
          dev-shell = devShell.inputDerivation;
        };
        devShell = pkgs.mkShell {
          name = "hammer-reweight";
          buildInputs = (with pkgs; with pythonPackages; [
            # Dev tools
            clang-tools

            root
            hammer-phys
            ff_calc

            # Python stack
            pythonPackages.pyTuplingUtils
          ]);

          FONTCONFIG_FILE = pkgs.makeFontsConf {
            fontDirectories = with pkgs; [
              corefonts
            ];
          };

          shellHook = ''
            export PATH=$(pwd)/bin:$(pwd)/utils:$PATH
          '';
        };
      });
}
