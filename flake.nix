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
        python = pkgs.python3;
        pythonPackages = python.pkgs;
      in
      rec {
        packages = flake-utils.lib.flattenTree {
          dev-shell = devShell.inputDerivation;
        };
        devShell = pkgs.mkShell rec {
          name = "hammer-reweight";
          buildInputs = (with pkgs; with pythonPackages; [
            root
            hammer-phys
            ff_calc
            clang-tools

            # Python virtualenv
            jedi
            flake8
            pylint
            virtualenvwrapper
          ]);

          shellHook = ''
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

            export PATH=$(pwd)/bin:$(pwd)/utils:$PATH
          '';
        };
      });
}
