# hammer-redist
HAMMER-related info.

Hammer uses `CMake` to configure out-of-source build. System requirements and
configurations are described in [Hammer documentation](https://hammer.physics.lbl.gov/readme.html).

Hammer has been compiled successfully on:

- Arch Linux, with `gcc 10.1.0`, `Python 3.8.3`, and `ROOT 6.20/04`.
- macOS with `nix`, with `clang 7.1.0`, `Python 3.7.6`, and `ROOT 6.18/04`.


## Install dependencies
Hammer requires the following dependencies, at least in our case:

- `boost`
- `yaml-cpp`
- `root`

### Arch Linux
In Arch, these dependencies can be installed via:

```
sudo pacman -S boost yaml-cpp root
```


## Build Hammer
Get source code from Hammer's official repository. We prefer the `development`
branch for now:
```
git clone -b development --single-branch https://gitlab.com/mpapucci/Hammer.git
```

Now, you just need to type:
```
make build
```


## More details on building Hammer

### A patch to make Hammer compile with newer `gcc`
At the time of documenting this, one minor edit of the source code was
necessary:

- Add `#include <string>` to `Hammer-source/include/Hammer/Math/Units.hh`.

This editing is included as a patch in:
```
./nix/overlay/hammer-phys/default.nix
```

And can be applied with:
```
cd Hammer
git apply ../nix/overlay/hammer-phys/add_missing_header.patch
```

### Explanations on configuration flags
See the comments below:

```
mkdir -p out
cd Hammer && \
cmake -S . -B build \
	-DCMAKE_INSTALL_PREFIX=../out \
	-DWITH_PYTHON=ON -DWITH_ROOT=ON \
	-DWITH_EXAMPLES=ON \
	-DINSTALL_EXTERNAL_DEPENDENCIES=ON \  # Install required Python packages
	-DMAX_CXX_IS_14=OFF && \  # on Arch, ROOT is compiled with c++17 flag, on NixOS, it is c++11.
cd build && \
make && \
make install
```

**Note**: The C++ flags need to be in sync with ROOT's flags!

You should then be able to find example code and their executables under
`./out/share/Hammer/examples`.
