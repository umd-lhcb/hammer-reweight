{ stdenv, fetchgit, boost, python, root }:

stdenv.mkDerivation rec {
  pname = "hammer-phys";
  version = "868fade5";  # Updated on Jun 14, 2020.

  src = fetchgit {
    url = "https://gitlab.com/mpapucci/Hammer.git";
    ref = "developmnt";
    rev = version;
    sha256 = "14j09s87926z4mffi9rdjqp4xnryml4q4vs7y8c02d7hahlwfn9f";
  };

  nativeBuildInputs = [ cmake makeWrapper ];
  buildInputs = [
    boost
    python
    root
  ];

  cmakeFlags = [
    "-Drpath=ON"
    "-DCMAKE_INSTALL_LIBDIR=lib"
    "-DCMAKE_INSTALL_INCLUDEDIR=include"
    "-DWITH_PYTHON=ON"
    "-DWITH_ROOT=ON"
  ];
}
