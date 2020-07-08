# jE /dʒeɪ iː/

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

and then can be found in doc/html/index.html

## Licence

BSD-3

## Acknowledgments

Development of this code up to v1.0.0 was supported by the Russian Science Foundation through Grant
No.  18-72-00121.
