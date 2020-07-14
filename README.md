# hammer-redist
HAMMER-related info.

## To compile Hammer
Hammer uses `CMake` to configure out-of-source build. System requirements and
configurations are described in https://hammer.physics.lbl.gov/readme.html.

It's been compiled successfully with `gcc 10.1.0`, with interfaces to `Python 3.8.3` and `ROOT 6.20/04`.


### Install Boost
If `boost` is not installed on the system, install it first.

```
tar xzvf tar.gz
cd boost_1_73_0/
./bootstrap.sh --prefix=path/to/install
./b2 install
```

And environment variable is defined:
`export BOOST_DIR=boost-install/lib`


### Build Hammer
Get source code from https://gitlab.com/mpapucci/Hammer.git.

You may want to switch to the `development` branch to pull in the most recent
developments.

At the time of documenting this, one minor edit of the source code was necessary:
add `#include <string>` to `Hammer-source/include/Hammer/Math/Units.hh`.

Then configure with CMake. The following command includes installing external
dependencies, Python binding and ROOT interface.

Later ROOT versions can be compiled with C++17 spec so we don't need to enforce
C++14.

```
mkdir Hammer-build
cd Hammer-build
cmake -DCMAKE_INSTALL_PREFIX=../Hammer-install -DWITH_PYTHON=ON -DWITH_ROOT=ON -DWITH_EXAMPLES=ON -DINSTALL_EXTERNAL_DEPENDENCIES=ON -DMAX_CXX_IS_14=OFF ../Hammer-source
make
make install
```

You should then be able to find example code and their executables under
`Hammer-install/share/Hammer/examples`.
