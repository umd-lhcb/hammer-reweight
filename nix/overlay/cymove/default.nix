{ stdenv, buildPythonPackage, fetchPypi, cython }:

buildPythonPackage rec {
  pname = "cymove";
  version = "1.0.2";

  src = fetchPypi {
    inherit pname version;
    extension = "tar.gz";
    sha256 = "79c1350db2f1f92a459b87ee9072ec0790faab233bfaeb73bf78a5caadd5aaa8";
  };

  propagatedBuildInputs = [
    cython
  ];

  doCheck = false;
}
