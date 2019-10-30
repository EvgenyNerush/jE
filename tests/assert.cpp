/* A collection of assertions */

#include <cassert>
#include <vector>
#include <iostream>
#include <random>
#include "../src/rng.hpp"
#include "../src/histogram.hpp"
#include "../src/radiation.hpp"

using namespace std;

int main(){
    
    { // Park-Miller RNG is the same as minstd_rand
        size_t n = 1000;

        vector<uint64_t> v1(n);
        uint64_t r = 1;
        for (auto &v: v1) {
            pm_rng(r);
            v = r;
        }

        vector<uint64_t> v2(n);
        minstd_rand g(1);
        for (auto &v: v2) {
            v = static_cast<uint64_t>(g());
        }

        assert(v1 == v2);
    }

    { // Examples of histogram_1D
        vector<double> xs1({0,2});
        vector<double> bs1({-1, 3});
        assert(histogram_1D(2, xs1) == bs1);

        vector<double> bs2({-1, 1, 3});
        assert(histogram_1D(1, xs1) == bs2);

        vector<double> xs2({6, 4, 2, 0});
        vector<double> bs3({-1, 3, 7});
        assert(histogram_1D(2, xs2) == bs3);
    }

    { // Test of trapezoidal rule of integration
        vector<double> t_nodes({0, 1, 2});
        auto f1 = function<double(double)>( [](double _) { return 1; } );
        assert(trap_rule(f1, t_nodes) == 2);

        auto f2 = function<double(double)>( [](double t) { return t; } );
        assert(trap_rule(f2, t_nodes) == 2);
    }

    cout << "assert.cpp: \x1b[32mpassed\x1b[0m\n";
    return 0;
}
