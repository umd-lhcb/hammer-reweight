{
  description = "Code for FF reweighting in HAMMER.";

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
          overlays = [ root-curated.overlay pyTuplingUtils.overlay self.overlay ];
        };
        python = pkgs.python3;
        pythonPackages = python.pkgs;
      in
      rec {
        packages = flake-utils.lib.flattenTree {
          dev-shell = devShell.inputDerivation;
          hammer-reweight = pkgs.hammer-reweight;
        };
        devShell = pkgs.mkShell rec {
          name = "hammer-reweight-dev";
          buildInputs = (with pkgs; with pythonPackages; [
            # Dev tools
            clang-tools

            root
            # hammer-phys
            hammer-phys-dev
            ff_calc
            cxxopts
            boost

            # Python stack
            pyyaml
            pylint
            pythonPackages.pyTuplingUtils
          ]);

          FONTCONFIG_FILE = pkgs.makeFontsConf {
            fontDirectories = with pkgs; [
              gyre-fonts
            ];
          };

          shellHook = ''
            export PATH=$(pwd)/bin:$(pwd)/utils:$PATH
            export MPLBACKEND=agg  # the backend w/o a UI
            export MPLCONFIGDIR=$(pwd)/.matplotlib
          '';
        };
      });
}
