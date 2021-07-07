// TODO: muon_emission (kappa on axis - wrong) -> examples
//
// (1) Emission spectrum of muon computed with classical_synchrotron_power_kernel (from
// radiation.hpp), in comparison with bks_synchrotron_td. Both the refractive index n = 1 and
// vacuum refractive index are used; gamma * B / B_S = 30,  kappa = 2, see also Fig. 4(b) in [I.I.
// Artemenko et al., New J. Phys. 22 (2020) 093072].
//
// (2) The probability computed with classical_synchrotron_power_kernel is
// integrated over the angle for a given frequency. Emission probability shouldn't be negative, and the contribution from the
// points where the probability is negative is compared with the overall probability.

#include <iostream>
#include <vector>
#include "../src/radiation.hpp"

using namespace std;

// linspace function NOT similar to one in python3
std::vector<double> linspace(double a, double b, size_t m) {
    double h = (b - a) / static_cast<double>(m);
    std::vector<double> xs(m);
	typename std::vector<double>::iterator x;
    double val;
    for (x = xs.begin(), val = a + 0.5 * h; x != xs.end(); ++x, val += h) {
        *x = val;
    }
    return xs;
}

//==== (1) ====//
double err = 0.001; // acceptable relative error

double mu = 207; // muon mass in electron masses
double b = 1e-3;
double chi  = 30;
double kappa = 2;

double gamma_mu = chi / b;
double omega_1 = kappa / (b * b);

double theta_max = 0.0001;
size_t n_theta = 35; // number of dots in angle range
vector<double> thetas( linspace(0, theta_max, n_theta) );

double dP_bks_1(double theta) {
    double ri = 1;
    auto f = bks_synchrotron_td(ri, mu, b, gamma_mu);
    return mu * omega_1 * b * f(make_tuple(theta, omega_1)) / (2 * M_PI * mu * r(gamma_mu));
}

double dP_bks_vac(double theta) {
    double ri = 0.5 * real( vacuum_refractive_index_perp    (b, omega_1)
                          + vacuum_refractive_index_parallel(b, omega_1) );
    auto f = bks_synchrotron_td(ri, mu, b, gamma_mu);
    return mu * omega_1 * b * f(make_tuple(theta, omega_1)) / (2 * M_PI * mu * r(gamma_mu));
}

double dP_ker_1(double theta) {
    double ri = 1;
    return classical_synchrotron_power_kernel(ri, mu, b, gamma_mu, theta, omega_1);
}

double dP_ker_vac(double theta) {
    double ri = 0.5 * real( vacuum_refractive_index_perp    (b, omega_1)
                          + vacuum_refractive_index_parallel(b, omega_1) );
    return classical_synchrotron_power_kernel(ri, mu, b, gamma_mu, theta, omega_1);
}

//==== (2) ====//
double err_2 = 0.01; // acceptable relative error

double chi_e = 1e7;
double kappa_2 = 200;
double theta_2 = 2e-4;

double gamma_e = chi_e / b;
double omega_2 = kappa_2 / (b * b);
vector<double> thetas_2( linspace(0, theta_2, n_theta) );

double dP_ker_vac_2(double theta) {
    std::complex<double> ri = 0.5 * ( vacuum_refractive_index_perp    (b, omega_1)
                                    + vacuum_refractive_index_parallel(b, omega_1) );
    return classical_synchrotron_power_kernel(ri, 1, b, gamma_e, theta, omega_2);
}

int main() {

    double acc_bks_1   = 0;
    double acc_bks_vac = 0;
    double acc_ker_1   = 0;
    double acc_ker_vac = 0;

    for (auto theta: thetas) {
        acc_bks_1   += dP_bks_1  (theta);
        acc_bks_vac += dP_bks_vac(theta);
        acc_ker_1   += dP_ker_1  (theta);
        acc_ker_vac += dP_ker_vac(theta);
    }

    // ==== //

    double acc_2 = 0;
    double acc_2_neg = 0;

    for (auto theta: thetas_2) {
        double dP = dP_ker_vac_2(theta);
        acc_2 += dP;
        if (dP < 0) {
            acc_2_neg -= dP;
        }
    }

    // ==== //

    if (  fabs(acc_bks_1   / acc_ker_1   - 1) < err
       && fabs(acc_bks_vac / acc_ker_vac - 1) < err
       && acc_2_neg / acc_2 < err_2 ) {
        cout << "classical_synchrotron_power_kernel test: \x1b[32mpassed\x1b[0m\n";
    } else {
        cout << "classical_synchrotron_power_kernel test: \x1b[1;31mfailed\x1b[0m"
             << '\n' << acc_bks_1
             << '\n' << acc_bks_vac
             << '\n' << acc_ker_1
             << '\n' << acc_ker_vac
             << '\n' << acc_2 << '\t' << acc_2_neg
             << '\n';
    }
    return 0;
}
