/* Let's check that the synchrotron spectrum computed numerically with #emission_probability_s from
 * ../src/radiation.hpp coincides with the classical synchrotron formula.
 */

#include <iostream>
#include <vector>
#include <iomanip>
#include "../src/radiation.hpp"

using namespace std;

int main() {
    double accepted_relative_error = 1e-2;
    double gamma_e = 10'000;
    double b = 1e-4; // formally this yield chi = 1, but it is ok here because we compare two
                     // classical formulas with each other
    double oc = omega_c(gamma_e);

    vector<double> thetas;
    for (int i = 0; i < 3;  ++i) {
        thetas.push_back(0.5 * static_cast<double>(i) / gamma_e);
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
            double err = abs(num_val - theor_val) / theor_val;
            acc = acc and (err < accepted_relative_error);
            theor_vals.push_back(theor_val);
            num_vals.push_back(num_val);
            errs.push_back(err);
        }
    }

    if (acc) {
        cout << "synchrotron_radiation test: \x1b[32mpassed\x1b[0m\n";
    } else {
        cout << "synchrotron_radiation test: \x1b[1;31mfailed\x1b[0m\n";
        cout << "accepted error = " << accepted_relative_error << '\n';
        int i = 0;
        for (auto theta: thetas) {
            cout << "theta = " << theta * gamma_e << " / gamma_e" << '\n';
            cout << "om/om_c     Jackson     JacksonNum  rel.error\n";
            for (auto omega: omegas) {
                cout << setw(12) << left << omega / oc  
                     << setw(12) << left << theor_vals[i]
                     << setw(12) << left << num_vals[i]
                     << setw(12) << left << errs[i]
                     << '\n';
                ++i;
            }
        }
    }

    return 0;
}
