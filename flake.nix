{
  description = "MC reweighting code with HAMMER.";

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
      {
        devShell = pkgs.mkShell {
          name = "hammer-redist";
          buildInputs = with pkgs; [
            root
            hammer-phys
            #ff_calc
          ];
        };
      });
}
