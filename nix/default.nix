{ stdenv
, root
, hammer-phys
, cxxopts
}:

stdenv.mkDerivation {
  pname = "hammer-reweight";
  version = "0.4.2";

  src = builtins.path { path = ./..; name = "hammer-reweight"; };

  buildInputs = [ root hammer-phys cxxopts ];

  installPhase = ''
    mkdir -p $out/bin
    cp bin/* $out/bin
  '';
}
