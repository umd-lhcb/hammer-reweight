{ stdenv
, makeWrapper, cmake, root }:

stdenv.mkDerivation rec {
  pname = "ff_calc";
  version = "1.1";

  src = builtins.path { path = ./../../../validation/ff_calc; name = "ff_calc"; };

  postInstall = ''
    cp -r ${src}/inc $out/include
  '';

  nativeBuildInputs = [ makeWrapper cmake root ];
}
