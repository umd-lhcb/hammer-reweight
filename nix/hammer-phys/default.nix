{ stdenv
, cmake
, makeWrapper
, pkgconfig
, boost
, libyamlcpp
, root
}:

stdenv.mkDerivation rec {
  pname = "hammer-phys";
  version = "f7827bf"; # Updated on Jun 17, 2020, 03:17 CST.

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
    root
  ];

  patches = [ ./add_missing_header.patch ];

  cmakeFlags = [
    "-DCMAKE_INSTALL_LIBDIR=lib"
    "-DCMAKE_INSTALL_INCLUDEDIR=include"
    "-DBUILD_SHARED_LIBS=ON"
    "-DWITH_PYTHON=OFF"
    "-DWITH_ROOT=ON"
    "-DINSTALL_EXTERNAL_DEPENDENCIES=OFF"
    "-DMAX_CXX_IS_14=ON" # on NixOS, ROOT is compiled with c++11 flag.
  ];

  # Move the .so files to the lib folder so the output looks like this:
  #   lib/*.so
  # instead of:
  #   lib/Hammer/*.so
  postFixup = ''
    mv $out/lib/Hammer/* $out/lib
    rm -rf $out/lib/Hammer
  '';
}
