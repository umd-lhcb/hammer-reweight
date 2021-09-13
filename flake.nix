{
  description = "Code for FF reweighting in  HAMMER.";

  inputs = {
    root-curated.url = "github:umd-lhcb/root-curated";
    nixpkgs.follows = "root-curated/nixpkgs";
    flake-utils.follows = "root-curated/flake-utils";

    #pyTuplingUtils.url = "github:umd-lhcb/pyTuplingUtils";
    # FIXME: Can't use the flake above because of a "follows" problem
    #        For more details, see this:
    #          https://github.com/NixOS/nix/issues/3602
    #        This problem is fixed in a very recent nix release (cira Sep 2021),
    #        but it's a hassle to update so let's use a workaround instead
  };

  outputs = { self, nixpkgs, flake-utils, root-curated }:
    {
      overlay = import ./nix/overlay.nix;
    } //
    flake-utils.lib.eachDefaultSystem (system:
      let
        pkgs = import nixpkgs {
          inherit system;
          config = { allowUnfree = true; };
          #overlays = [ root-curated.overlay self.overlay pyTuplingUtils.overlay ];
          overlays = [ root-curated.overlay self.overlay ];
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
            hammer-phys
            ff_calc

            # Python stack
            #pythonPackages.pyTuplingUtils
            virtualenvwrapper
          ]);

          FONTCONFIG_FILE = pkgs.makeFontsConf {
            fontDirectories = with pkgs; [
              corefonts
            ];
          };

          shellHook = ''
            export PATH=$(pwd)/bin:$(pwd)/utils:$PATH

            # Allow the use of wheels.
            SOURCE_DATE_EPOCH=$(date +%s)

            if test -d $HOME/build/python-venv; then
              VENV=$HOME/build/python-venv/${name}
            else
              VENV=./.virtualenv
            fi

            if test ! -d $VENV; then
              virtualenv $VENV
            fi
            source $VENV/bin/activate

            # allow for the environment to pick up packages installed with virtualenv
            export PYTHONPATH=$VENV/${python.sitePackages}/:$PYTHONPATH

            # fix libstdc++.so not found error
            export LD_LIBRARY_PATH=${pkgs.stdenv.cc.cc.lib}/lib:$LD_LIBRARY_PATH
          '';
        };
      });
}
