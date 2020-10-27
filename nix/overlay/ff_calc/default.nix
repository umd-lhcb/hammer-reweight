{ stdenv , makeWrapper }:

stdenv.mkDerivation {
  pname = "ff_calc";
  version = "2020-10-27";

  src = builtins.path { path = ./../../../validation/ff_calc; name = "ff_calc"; };

  postInstall = ''
    cp -r build/lib $out/lib
    cp -r inc $out/include
  '';

  nativeBuildInputs = [ makeWrapper ];
}
