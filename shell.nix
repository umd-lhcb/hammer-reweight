let
  pkgs = import <nixpkgs> { overlays = [(import ./nix/overlay)]; };
  stdenv = pkgs.stdenv;
in

pkgs.mkShell {
  name = "hammer-redist";
  buildInputs = with pkgs; [
    hammer-phys
  ];

  propagatedBuildInputs = with pkgs; [
    boost
    libyamlcpp
    root
  ];
}
