{ stdenv, cmake, root }:

stdenv.mkDerivation rec {
  pname = "ff_calc";
  version = "1.2";

  src = builtins.path { path = ./../../validation/ff_calc; name = pname; };

  postInstall = ''
    cp -r ${src}/inc $out/include
  '';

  nativeBuildInputs = [ cmake root ];
}
