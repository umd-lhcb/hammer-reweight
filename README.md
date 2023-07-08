# hammer-reweight
Code for FF reweighting in HAMMER.


## Installation

Install `nix` with flake support, then type in `nix develop`. All dependencies
will be installed automatically for you.

To build all binaries:

```
make
```

Once this is done, to generate weight ntuples (using the nominal reweighting script 
`src/ReweightRDX.cpp`) based on sample ntuples (stored in `samples/`) from RDX run 2 analysis:

```
make rdx-run2-ntuples
```

Note that the `stdout` from the reweighter will also be saved as log files.

To plot resulting sample `q2` (stored in `gen/`):

```
make sample-plots
```

To generate some reweighting validation plots (stored in `gen`):
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

Note that all particle IDs are taken to be absolute value.

### `ReweightRDX`

To do the actual form factor reweighting for RDX analyses, you can do something like this:
```
ReweightRDX samples/rdx-run2-Bd2DstMuNu.root output.root TupleB0/DecayTree run2
```
This will generate weight ntuples, stored in the home directory (and named `output.root`). Notably,
you can also run other reweighting scripts (stored in `src/`) or the nominal reweighting with more
printouts for debugging with similar commands (eg. `ReweightRDXDstNoCorr`, or `'ReweightRDXDebug`
for debugging, instead of `ReweightRDX`). The additional reweighting scripts are mostly used for 
studies.

The `output.root` contains a `w_ff` branch and some other debugging branches.

NOTE: currently, the FF parameters/errors/variations set in the nominal reweighter `ReweightRDX`
have

- For $B\rightarrow D$: parameter values and errors taken from [1606.08030](https://arxiv.org/abs/1606.08030),
with one variation (along the direction of the nominal FF parameters, ie. the rescaling direction)
removed. Correlations from the paper are taken into account.
- For $B\rightarrow D^{\*}$: parameter values and errors taken from [2105.14019](https://arxiv.org/abs/2105.14019)
(v3, Table 12- "Lattice + both"), with one variation (along the rescaling direction) removed. Correlations
from the paper are taken into account.
- For $B\rightarrow D\_{(s)}^{\*\*}$: parameter values and errors taken from [1711.03110](https://arxiv.org/abs/1711.03110),
with one variation (in the direction of parameter that mainly effect normalization/rescaling: $\tau (1)$
and $\zeta (1)$) removed for both the wide/narrow states (which are parametrized separately). The $D_s^{\*\*}$
FFs are taken to be the same as for the $D^{\*\*}$. Previous fit results indicated that the paper errors were
too constraining, and potentially the central values ought to be shifted; currently, the errors used for
variations are _double_ the paper errors, and the central values are _not_ shifted. Correlations from
the paper are _not_ taken into account; the variations are exactly $\mathbf{f}\_{nom} +- \mathbf{sigma}\_i$
(where $\mathbf{sigma}\_i$ is a vector of 0's and the $i^{th}$ FF param error in the $i^{th}$ position).

### `ValidateRDX`

This is used to generate some toy data to validate HAMMER reweighting for RDX
with an independently implemented FF calculator.

To use it:
```
ValidateRDX test.root
```

### Compute FF variation parameters

To produce code that can be pasted into a reweighting script (eg. `ReweightRDX.cpp`) that specifies
the parameter values as well as the FF variations, use

```
make ff-params-RDX
```

See [`utils/gen_ham_params.py`](./utils/gen_ham_params.py) and [`spec/rdx-run2.yml`](./spec/rdx-run2.yml)
for how these are computed.

To produce the variations with rescaling removed from the variations (not nominally used for 
$D_{(s)}^{\*\*}$), use

```
make ff-params-RDX-no-rescale
```

The motivation and math underlying removing the rescaling from the variations can be found in
[`utils/gen_ham_params_no_rescale.py`](./utils/gen_ham_params_no_rescale.py).

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

### Weight computation

1. The amplitude tensor is computed
2. `Hammer::getWeight("FF_scheme")` is called.
3. In that method, find the base weights of the event: `double result = _event->getEventBaseWeight()`.

    This is by default set to 1 and we typically don't change it

4. For each _process_, multiply the base weight by `result *= _event->getWeight(scheme, elem)`.

    Here _process_ means decay processes added by user, like `B -> D* l nu`.
    Also, the amplitude tensor is contracted with the external eigenvectors.

    This is the step where the FF variation is effected.

### Form factor variation

1. Find/compute the covariance matrix $U$ of a FF parameterization from some reference paper
2. Find the eigenvectors and eigenvalues of $U$
3. Define a matrix $M$, with each _row_ being one of the eigenvectors
4. Define a matrix $C$, with $C_{ij} \equiv \delta_{ij} \frac{1}{\lambda_i}$
5. The transformation $A \equiv CM$ maps the parameterization basis to an
   _orthonormal_ error eigen basis.
6. We can either define $A^{-1}$ in HAMMER, or define variations in terms of $A^{-1} \delta$,
    where $\delta$ is some variation in the error eigen basis.

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
