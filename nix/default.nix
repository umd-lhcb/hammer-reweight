{ stdenv
, root
, hammer-phys
}:

stdenv.mkDerivation {
  pname = "hammer-reweight";
  version = "0.3.5";

  src = builtins.path { path = ./..; name = "hammer-reweight"; };

  buildInputs = [ root hammer-phys ];

  installPhase = ''
    mkdir -p $out/bin
    cp bin/* $out/bin
  '';
}
