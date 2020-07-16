{ stdenv, fetchgit
, cmake, makeWrapper, pkgconfig
, boost, python3, libyamlcpp, root }:

stdenv.mkDerivation rec {
  pname = "hammer-phys";
  version = "f7827bf";  # Updated on Jun 17, 2020, 03:17 CST.

  src = fetchgit {
    url = "https://gitlab.com/mpapucci/Hammer.git";
    branchName = "developmnt";
    rev = version;
    sha256 = "0dnv2vwh1fz61gl60x7bymkd4mxi16523jygwj7b1135lqsbms2z";
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
    setuptools
  ];

  patches = [ ./add_missing_header.patch ];

  cmakeFlags = [
    "-DCMAKE_INSTALL_PREFIX=$out"
    "-DBUILD_SHARED_LIBS=ON"
    "-DWITH_PYTHON=ON"
    "-DWITH_ROOT=ON"
    "-DINSTALL_EXTERNAL_DEPENDENCIES=OFF"
    "-DMAX_CXX_IS_14=ON"  # on NixOS, ROOT is compiled with c++11 flag.
  ];
}
