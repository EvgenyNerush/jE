#include <cmath>
#include <iostream>
#include <vector>
#include <tuple>
#include "../../src/histogram.hpp"
#include "../../src/metropolis.hpp"
#include "../../src/rng.hpp"

using namespace std;

double td(double x){
    if (x > 0 and x <= 1) {
        return 1 / sqrt(x);
    } else {
        return 0;
    }
}

void proposal_density(pair<uint64_t, double>& rx) {
    pm_rng(rx.first);
    double d = pm_cast_to_01(rx.first);
    rx.second = rx.second + 0.1 * (2 * d - 1);
}

vector<double> xs =
    metropolis<double, uint64_t, uint64_t> ( pm_rng
                                           , 12345ul
                                           , pm_cast_to_01
                                           , function<void(pair<uint64_t, double>&)>
                                             (proposal_density)
                                           , make_pair(54321ul, 0.5)
                                           , td
                                           , 16'000'000);

vector<double> bs = histogram_1D(1'600'000, xs);

extern "C" {
size_t get_n() {
    return xs.size();
}

size_t get_nb() {
    return bs.size();
}

// random points generated with metropolis algorithm
double* points() {
    double* ys = new double[xs.size()];
    for (size_t i = 0; i < xs.size(); ++i) {
        ys[i] = xs[i];
    }
    return ys;
}

// boundaries of bins for the given random points
double* bins() {
    double* ys = new double[bs.size()];
    for (size_t i = 0; i < bs.size(); ++i) {
        ys[i] = bs[i];
    }
    return ys;
}
}
