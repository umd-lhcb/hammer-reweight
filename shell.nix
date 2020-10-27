let
  pkgs = import <nixpkgs> { overlays = [(import ./nix/overlay)]; };
  stdenv = pkgs.stdenv;
in

pkgs.mkShell {
  name = "hammer-redist";
  propagatedBuildInputs = with pkgs; [
    hammer-phys
    boost
    libyamlcpp
    root
    ff_calc
  ];
}
