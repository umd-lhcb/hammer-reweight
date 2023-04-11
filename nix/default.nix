{ stdenv
, root
# , hammer-phys
, hammer-phys-dev
, cxxopts
}:

stdenv.mkDerivation {
  pname = "hammer-reweight";
  version = "0.5.0";

  src = builtins.path { path = ./..; name = "hammer-reweight"; };

  buildInputs = [ root hammer-phys-dev cxxopts ];

  installPhase = ''
    mkdir -p $out/bin
    cp bin/* $out/bin
  '';
}
