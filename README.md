# hammer-redist
MC reweighting code with HAMMER.

Hammer uses `CMake` to configure out-of-source build. System requirements and
configurations are described in [Hammer documentation](https://hammer.physics.lbl.gov/readme.html).

For old instructions on installing Hammer on non-NixOS Linux distributions,
please take a look at [version `0.1`](https://github.com/umd-lhcb/hammer-redist/tree/0.1).


## More details on building Hammer

### A patch to make Hammer compile with newer `gcc`
At the time of documenting this, one minor edit of the source code was
necessary:

Add:
```cpp
#include <cstring>
#include <climits>
```
to `include/Hammer/Math/Units.hh`.

This editing is included as a patch in:
```
./nix/overlay/hammer-phys/add_missing_header.patch
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
