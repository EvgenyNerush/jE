/* Here dI/d\omega computed numerically with #bks_td from ../src/radiation.hpp is compared with the
 * analytical result found in [V. B.  Berestetskii, E.  M. Lifshitz and L. P.  Pitaevskii, Quantum
 * Electrodynamics, Pergamon, New York, 1982].
 *
 * The parameters are chosen such that the spin term is especially important.
 */

#include <iostream>
#include "../src/radiation.hpp"

using namespace std;

int main() {
    double acc_err = 0.01; // acceptable error
    double gamma_e = 1e5;
	double chi     = 10;
    double b       = chi / gamma_e;
    double xi      = 0.9;              // hbar omega / (mc^2 gamma_e)
    double omega   = xi * gamma_e / b;

    // emission power per theta and omega unit intervals, normalized to mc^2
    auto p = function<double(double)> (
        [=](double theta) {
            auto f = bks_synchrotron_td(1, 1, b, gamma_e);
            return omega * b * f(make_tuple(theta, omega)) / (2 * M_PI * r(gamma_e));
        }
    );

    size_t n_theta = 50; // number of nodes for integration over theta
    double theta_m = 4 / gamma_e;
    double dtheta = theta_m / static_cast<double>(n_theta);
    double dPdomega_num = 0; // emission power per unit frequency interval, normalized to mc^2
    for (size_t i = 0; i < n_theta; ++i) {
        dPdomega_num += 2 * dtheta * p(0.5 * dtheta + dtheta * static_cast<double>(i));
    }

    double dPdomega_theor = 1.128e-07; /* value computed with external program from equation found
                                          in [V. B.  Berestetskii, E.  M. Lifshitz and L. P.
                                          Pitaevskii, Quantum Electrodynamics, Pergamon, New York,
                                          1982]
                                        */

    if (fabs(dPdomega_num / dPdomega_theor - 1) < acc_err) {
        cout << "spin_term test: \x1b[32mpassed\x1b[0m\n";
    } else {
        cout << "spin_term test: \x1b[1;31mfailed\x1b[0m\n";
        cout << "acc. error  = " << acc_err                                 << '\n'
             << "rel. error  = " << fabs(dPdomega_num / dPdomega_theor - 1) << '\n'
             << "theor value = " << dPdomega_theor                          << '\n'
             << "num. value  = " << dPdomega_num                            << '\n';
    }

    return 0;
}
