/* In this test *n* x-coordinates with 1/sqrt(x) distribution are generated, then the average
 * value of x, x^2 and x^3 are computed and compared with theoretical values.
 * This test takes about 2.2 seconds with Xeon X5550 processor.
 */

#include <cmath>
#include <iostream>
#include <chrono>
#include "../src/metropolis.hpp"
#include "../src/rng.hpp"

using namespace std;
using namespace std::chrono;

// Target distribution used in this test
double td(double x){
    if (x > 0 and x <= 1) {
        return 1 / sqrt(x);
    } else {
        return 0;
    }
}

// Function to update (rng_state, x), with new x evenly distributed in old x +- 0.1
void proposal_density(pair<uint64_t, double>& rx) {
    pm_rng(rx.first);
    double d = pm_cast_to_01(rx.first);
    rx.second = rx.second + 0.1 * (2 * d - 1);
}

// mean value of x^p, i.e. int_0^1 x^(p - 1/2) dx / \int_0^1 x^(-1/2) dx. p shoud be greater than
// -1/2.
double theor_mean_x_pow(int p){
    return 0.5 / (p + 0.5);
}

// mean value of x^p, computed numerically
double num_mean_x_pow(vector<double> xs, int p){
    double acc = 0;
    for (auto x: xs) {
        acc += pow(x, p);
    }
    return acc / xs.size();
}

int main() {
    size_t n = 16'000'000; // should be much less than pm_randmax
    double acceptable_accuracy = 4 / sqrt(static_cast<double>(n));

    steady_clock::time_point t1 = steady_clock::now();
    vector<double> xs
        = metropolis<double, uint64_t, uint64_t>
            ( pm_rng
            , 12345ul
            , pm_cast_to_01
            , function<void(pair<uint64_t, double>&)>(proposal_density)
            , make_pair(54321ul, 0.5)
            , td
            , n);

    vector<int> powers {1, 2, 3};
    vector<double> theor_vs;
    vector<double> num_vs;
    vector<double> absolute_errors;
    bool acc = true;

    for (auto p: powers) {
        double theor_v = theor_mean_x_pow(p);
        double num_v   = num_mean_x_pow(xs, p);
        theor_vs.push_back(theor_v);
        num_vs.push_back(num_v);

        double absolute_error = fabs(num_v - theor_v);
        absolute_errors.push_back(absolute_error);
        acc = acc and absolute_error < acceptable_accuracy;
    }
    steady_clock::time_point t2 = steady_clock::now();

    if (acc) {
        cout << "metropolis test: \x1b[32mpassed\x1b[0m\n";
    } else {
        cout << "metropolis test: \x1b[1;31mfailed\x1b[0m\n";
        cout << "acceptable absolute error value = " << acceptable_accuracy << '\n';
        cout << "power"          << '\t'
                  << "theor value"    << '\t'
                  << "num value"      << '\t'
                  << "absolute error" << '\n';
        for (size_t i = 0; i < powers.size(); ++i) {
            cout << powers[i]          << '\t'
                      << theor_vs[i]        << '\t'
                      << num_vs[i]          << '\t'
                      << absolute_errors[i] << '\n';
        }
    }

    duration<double> time_span = duration_cast<duration<double>>(t2 - t1);
    cout << "test takes " << time_span.count() << " seconds\n";
    return 0;
}
