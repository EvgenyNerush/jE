/* Let's check that the synchrotron spectrum computed numerically with #emission_probability_s from
 * ../src/radiation.hpp coincides with the classical synchrotron formula.
 * This test takes about 0.9 second with Xeon X5550 processor.
 */

#include <iostream>
#include <vector>
#include <chrono>
#include "../src/radiation.hpp"

using namespace std;
using namespace std::chrono;

int main() {
    // error relative to value of the maximum of the synchrotron spectrum
    double accepted_error = 2e-5;
    double gamma_e = 1'000;
    double b = 1e-4;
    double oc = omega_c(gamma_e);
    steady_clock::time_point t1 = steady_clock::now();
    // approximate value of the maximum of the synchrotron spectrum
    double max = jackson1483(b, gamma_e, 0, 0.42 * oc);

    vector<double> thetas;
    for (int i = 0; i < 3;  ++i) {
        thetas.push_back(static_cast<double>(i) / gamma_e);
    }
    vector<double> omegas;
    for (double x = 0.025; x < 2; x *= 4) {
        omegas.push_back(x * oc);
    }

    bool acc = true;
    vector<double> theor_vals;
    vector<double> num_vals;
    vector<double> errs;
    for (auto theta: thetas) {
        for (auto omega: omegas) {
            double theor_val = jackson1483    (b, gamma_e, theta, omega);
            double num_val   = jackson1483_num(b, gamma_e, theta, omega);
            double err = abs(num_val - theor_val) / max;
            acc = acc and (err < accepted_error);
            theor_vals.push_back(theor_val);
            num_vals.push_back(num_val);
            errs.push_back(err);
        }
    }
    steady_clock::time_point t2 = steady_clock::now();

    if (acc) {
        cout << "synchrotron_radiation test: \x1b[32mpassed\x1b[0m\n";
    } else {
        cout << "synchrotron_radiation test: \x1b[1;31mfailed\x1b[0m\n";
        cout << "accepted error = " << accepted_error << '\n';
        int i = 0;
        for (auto theta: thetas) {
            cout << "theta = " << theta * gamma_e << " / gamma_e" << '\n';
            cout << "om/om_c\tJackson \tJacksonNum\terror\n";
            for (auto omega: omegas) {
                cout << omega / oc    << '\t'
                     << theor_vals[i] << '\t'
                     << num_vals[i]   << '\t'
                     << errs[i]       << '\n';
                ++i;
            }
        }
    }

    duration<double> time_span = duration_cast<duration<double>>(t2 - t1);
    cout << "test takes " << time_span.count() << " seconds\n";
    return 0;
}
