# jE /dʒeɪ iː/ ![jay logo](doc/figures/jay-logo.png)

This code implements:

  * Metropolis algorithm
  * Baier-Katkov general formula for the photon emission by charged particles with spin 1/2. The
    refractive index, radiation recoil and spin flips can be taken into account
  * Formulas for fast computation of synchrotron-Cherenkov radiation
  * A number of tests

## How to

_jE_ provides header files (see `src` directory) which can be included easily into one's project.
Note however that _jE_ depend on [Erich Niebler's ranges](https://github.com/ericniebler/range-v3).
Therefore, to build the code one need to get _ranges_ first, e.g. to build and run the tests try:

    cd jE
    git clone https://github.com/ericniebler/range-v3.git
    cd tests
    make

## Documentation

Html documentation can be produced with Doxygen:

    cd doc
    doxygen doxygen.conf

and then can be found in doc/html/index.html. See online version for OLD code (v 1.0.0) [here](https://evgenynerush.github.io/jE-doc/).

## Licence

BSD-3
