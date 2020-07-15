{ stdenv, fetchgit
, cmake, makeWrapper, pkgconfig, libyamlcpp
, boost, python3, root }:

stdenv.mkDerivation rec {
  pname = "hammer-phys";
  version = "868fade5";  # Updated on Jun 14, 2020.

  src = fetchgit {
    url = "https://gitlab.com/mpapucci/Hammer.git";
    branchName = "developmnt";
    rev = version;
    sha256 = "14j09s87926z4mffi9rdjqp4xnryml4q4vs7y8c02d7hahlwfn9f";
  };

  nativeBuildInputs = [ cmake makeWrapper pkgconfig ];
  buildInputs = [
    boost
    libyamlcpp
    python3
    root
  ];
  propagatedBuildInputs = with python3.pkgs; [
    cython
    cymove
    numpy
    matplotlib
  ];

  cmakeFlags = [
    "-Drpath=ON"
    "-DCMAKE_INSTALL_LIBDIR=lib"
    "-DCMAKE_INSTALL_INCLUDEDIR=include"
    "-DWITH_PYTHON=ON"
    "-DWITH_ROOT=ON"
    "-DMAX_CXX_IS_14=ON"  # on NixOS, ROOT is compiled with c++11 flag.
  ];
}
