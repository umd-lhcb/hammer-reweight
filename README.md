# hammer-redist
Redistribute `HAMMER` and related tools.

## Included packages

### `HAMMER`

Hammer uses `CMake` to configure out-of-source build. System requirements and
configurations are described in [Hammer documentation](https://hammer.physics.lbl.gov/readme.html).

For old instructions on installing `HAMMER` on non-NixOS Linux distributions,
please take a look at [version `0.1`](https://github.com/umd-lhcb/hammer-redist/tree/0.1).

Currently, [this patch](./nix/hammer-phys/add_missing_header.patch) is required
for `HAMMER` to compile. To apply that patch:

```
git apply <path_to_the_patch>
```
