# hammer-reweight
Code for FF reweighting in HAMMER.


## Installation

Install `nix` with flake support, then type in `nix develop`. In the resulting
shell prompt you can do:

```
make
```

To generate weight ntuples based on sample ntuples from RDX run 2 analysis:

```
make rdx-run2-ntuples
```

Note that the `stdout` from the reweighter will also be saved as log files.

To compile programs (which will be put in the `bin` folder) and

```
make sample-plots
```

To generate sample reweighting plots (which will be put in the `gen` folder).

You can generate some reweighting validation plots with:
```
make validation-plots
```


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


## HAMMER tips

### General HAMMER reweighting workflow

1. Configure the allowed decay modes and input/output FF parameterizations
2. Define `Hammer::Particle` and adding these particles to a newly created
    `Hammer::Process`, denote `proc.addParticle(p)`
3. Specify the decay process tree with `proc.addVertex(mother, {list_of_daughters})`
4. Initialize a new HAMMER event object with `Hammer::initEvent`
5. Add a `proc` to the newly initialized event with `Hammer::addProcess(proc)`.
    This includes:
    1. Prune soft photons (the 4-momenta will be modified in this step)
        by adding them to nearest charged particle in the same decay tree level
        (polar-angle wise) then remove the photons from the tree

        This process may alter the kinematics in a bad way such that HAMMER
        won't process the event further (either a `NaN` in the cosine or a
        negative invariant mass)
    2. Compute the decay tree hash ID
    3. Cache particle dependencies so the decay tree only need to be traversed
        once (a technicality)
6. Compute the amplitude tensor with `Hammer::processEvent()`, then get the weight
    with `Hammer::getWeight("FF_scheme")`:
    - A bit more details on `processEvent`: it's defined in `src/Core/Hammer.cc`
    - Then `Event.calc()` is called, which is defined in `src/core/Event.cc`
    - Inside `Event.calc()`, it calls `p.second.calc()` for each element in `_process`.

        Note: `_process` is of type `ProcIdDict`, which is `std::map<ProcessUID, T>`.
    - The actual calculation happens `src/core/ProcManager.cc`.

### Soft photons (very)

HAMMER doesn't like very soft photons (for energy on the order of `1e-10`),
and for these photons, their effect are negligible as well. So don't add these photons
to the HAMMER's decay tree in the first place!

### Vertex level momentum conservation

HAMMER doesn't require momentum to be conserved in general; it only requires each particle
has a non-negative invariant mass.

However, for soft photon removal, the vertex-level momentum conservation is enforced. This
sometimes can lead to undesired effect of negative invariant mass

### Negative invariant mass

The invariant mass squared is computed in `Hammer::FourMomentum.mass2()`. This
is just `E^2 - p^2`. Due to float-point arithmetic, sometimes this is slightly
negative, which is quite common for particles with very small invariant masses.

HAMMER provides a more robust way, The `Hammer::FourMomentum.mass()` method:
```cpp
double FourMomentum::mass() const {
    if (fuzzyLess(mass2() , 0.0)) { // less than 'a little bit negative'
        throw Error("Rosebud! Negative mass^2: " + to_string(mass2()));
    }
    return sqrt(fabs(mass2()));
}
```

The `fuzzyLess` is defined as:
```cpp
static const double precision = 0.001;

inline bool fuzzyLess(const double val1, const double val2) {
    return (val1-val2 < -1.*std::max(precision, std::numeric_limits<double>::min()));
}
```

Only if the invariant mass is strictly negative, HAMMER would complains about it.

### HAMMER produces a segmentation fault

This has happened to us on **some of** the `D*_2` candidates. The reason was
that we also added the daughters of the `D*_2` to the decay tree, and HAMMER
really doesn't like them.

By removing them from the tree, HAMMER runs fine. So it's probably **a good idea
to not adding unneeded particles to the tree**.
