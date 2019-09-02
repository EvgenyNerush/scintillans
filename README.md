## Note

The latest version 0.3.0 includes synchrotron emission only; work on pair photoproduction is in
progress.

## About

The code _Scintillans_ simulates evolution of electron, photon and positron energy spectra (fₑ, fᵧ
and fₚ of ϵ) in a constant magnetic field. _Scintillans_ provides solver of 1D Boltzmann's
equations which take into account synchrotron emission and pair photoproduction according to
well-known QED formulas.  For details, see [I I Artemenko et al 2019 Plasma Phys. Control. Fusion
61 074003](https://doi.org/10.1088/1361-6587/ab1712) (or see freely available
[preprint](https://www.researchgate.net/publication/332283915_Global_constant_field_approximation_for_radiation_reaction_in_collision_of_high-intensity_laser_pulse_with_electron_beam)).
Source-code documentation (see below) describes the code in details as well.

## How to use

The easiest way to build Scintillans is through
[Stack](https://docs.haskellstack.org/en/stable/README/) which, in turn, can be installed directly
with package managers in most Linux operating systems. Thus, the following commands are enough in
most cases:

    git clone https://github.com/EvgenyNerush/scintillans.git
    cd scintillans
    stack build

Then a separate executable can be made. See `scintillans/examples/PhotonEmissionTest.hs` (or other
tests in the `examples` directory) which can be built and run with

	stack ghc -- -O2 --make PhotonEmissionTest.hs
    ./PhotonEmissionTest

Alternatively, one can include _Scintillans_ as a dependency in `*.cabal` and `stack.yaml` files of
his/her own project.

## Documentation

Source-code documentation can be produced in HTML format with Haddock:

    stack haddock

which will take some time for the first run. This command also prints the path to the generated
docs, e.g.

    ...
    Updating Haddock index for local packages in
    .../scintillans/.stack-work/install/.../8.4.3/doc/index.html

which can be then open with an internet browser.

## Acknowledgments

Development of this code is supported by the Russian Science Foundation through Grant No.
18-72-00121
