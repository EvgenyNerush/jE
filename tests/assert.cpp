/* A collection of assertions */

#include <cassert>
#include <vector>
#include <iostream>
#include <random>
#include "../src/rng.hpp"
#include "../src/histogram.hpp"
#include "../src/radiation.hpp"
#include "../src/proposal_density.hpp"

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

    { // Bisection to find the roots of f(x) = 0
        auto f = std::function<double(const double&)>( [](double x){ return x * x - 1; } );
        assert(bisection(f, -0.5, 0.5, 1) == std::nullopt);
        assert(bisection(f, -1, 0, 0).value() == -1);
        assert(bisection(f,  0, 1, 0).value() ==  1);
        assert(bisection(f,  0, 4, 0).value() ==  2);
        assert(bisection(f,  0, 4, 1).value() ==  2);
        assert(bisection(f,  0, 4, 2).value() ==  1);
        assert(bisection(f,  0, 4, 3).value() ==  1);
        assert(fabs(bisection(f, 0, 3, 10).value() - 1) < 3.0 / pow(2, 10));
    }

    { // Properties of vacuum refractive index in strong magnetic field; see figure 9 in
      // [McDonald K.T. et al., Proposal for experimental studies of nonlinear quantum
      // electrodynamics, Princeton U. preprint DOE ER, 1986]. Not that McDonald uses chi which is
      // a half of chi we use.
        double prec = 0.1; // precision
        auto n = [](double b, double omega) {
            return 4 * M_PI / (alpha * b * b) * (vacuum_refractive_index(b, omega) - 1);
        };
        assert(fabs(n(0.1, 0.1) * 45 / 14 - 1)                 < prec);
        assert(fabs(n(1, 70) / (-0.278 * pow(70, -4/3.0)) - 1) < prec);
        assert(fabs(n(1, 0.5) / 0.35 - 1)                      < prec);
        auto f = [=](double chi) { return n(1, chi); };
        assert(fabs(bisection(f, 0, 40, 10).value() / 17 - 1)  < 2 * prec);
    }

    { // Test zipWith from proposal_density.hpp
        auto t1 = std::make_tuple(1, 2, 3);
        auto t2 = std::make_tuple(1, 1, 1);
        auto z1 = zipWith<2, int, int, int>();
        z1( [](int x, int y){ return x + y; }, t1, t2 );
        assert(t1 == std::make_tuple(2, 3, 4));

        t1 = std::make_tuple(1, 2, 3);
        t2 = std::make_tuple(4, 5, 6);
        auto z2 = zipWith<1, int, int, int>();
        z2( [](int x, int y){ return x * y; }, t1, t2 );
        assert(t1 == std::make_tuple(4, 10, 3));
    }

    cout << "assertions: \x1b[32mpassed\x1b[0m\n";
    return 0;
}
