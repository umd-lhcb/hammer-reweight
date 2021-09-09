# hammer-reweight
Code for FF reweighting in HAMMER.


## Installation

Install `nix` with flake support, then type in `nix develop`. In the resulting
shell prompt you can do:

```
make
```

to compile programs (which will be put in the `bin` folder) and

```
make sample-plots
```

to generate sample reweighting plots (which will be put in the `gen` folder).


## Usage

### `PrintMCDecay`

This program will print all partially truth-matched MC candidates. For example:
```
PrintMCDecay ./samples/rdx-run2-Bd2DstMuNu.root TupleB0/DecayTree
```

Will give you an output like this:
```
Total number of candidates: 23841
Truth-matched candidates: 23339
Truth-matched fraction: 0.978944
======
The following decay has 23314 candidates.
Is Tau decay: 0
B meson ID: B0 (511)
First daughter ID: D*+ (413)
  First G-daughter ID: D0 (421)
  Second G-daughter ID: pi+ (211)
======
The following decay has 25 candidates.
Is Tau decay: 0
B meson ID: B0 (511)
First daughter ID: D*+ (413)
  First G-daughter ID: D0 (421)
  Second G-daughter ID: pi+ (211)
  Third G-daughter ID: gamma (22)
```

Note that all particle IDs are taken to be absolute value so they may be unphysical.

### `ReweightRDX`

Form factor reweighter for RDX analyses. You can do something like this:
```
ReweightRDX samples/rdx-run2-Bd2DstMuNu.root output.root TupleB0/DecayTree run2
```

The `output.root` contains a `w_ff` branch and some other debugging branches.

### `ValidateRDX`

This is used to generate some toy data to validate HAMMER reweighting for RDX
with an independently implemented FF calculator.

To use it:
```
ValidateRDX test.root
```
