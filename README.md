# hammer-redist
HAMMER-related info.

## Compile Hammer in a recent Linux distribution
Hammer uses `CMake` to configure out-of-source build. System requirements and
configurations are described in [Hammer documentation](https://hammer.physics.lbl.gov/readme.html).

It's been compiled successfully with `gcc 10.1.0`, with interfaces to `Python
3.8.3` and `ROOT 6.20/04`, on Arch Linux.

### Install dependencies on Arch Linux
```
sudo pacman -S boost yaml-cpp root
```

### Build Hammer
Get source code from Hammer's official repository. We prefer the `development`
branch for now:
```
git clone -b development --single-branch https://gitlab.com/mpapucci/Hammer.git
```

At the time of documenting this, one minor edit of the source code was necessary:
add `#include <string>` to `Hammer-source/include/Hammer/Math/Units.hh`.

Then configure with CMake. The following command includes installing external
dependencies, Python binding and ROOT interface.

**Note**: The C++ flags need to be in sync with ROOT's flags!

```
mkdir Hammer-build
cd Hammer-build
cmake \
    -DCMAKE_INSTALL_PREFIX=../Hammer-install \
    -DWITH_PYTHON=ON -DWITH_ROOT=ON \
    -DWITH_EXAMPLES=ON \
    -DINSTALL_EXTERNAL_DEPENDENCIES=ON \  # Install required Python packages
    -DMAX_CXX_IS_14=OFF \  # on Arch, ROOT is compiled with c++17 flag, on NixOS, it is c++11.
    ../Hammer-source
make
make install
```

Note that here we'll install a local copy of `yaml-cpp`, regardless of if it's
available globally.

You should then be able to find example code and their executables under
`Hammer-install/share/Hammer/examples`.
