// The kernel (classical_synchrotron_power_kernel) is compared with values computed with different
// numerical method.

#include <iostream>
#include <iomanip>
#include "../src/radiation.hpp"

using namespace std;
using namespace std::literals::complex_literals;

double err = 0.001; // acceptable relative error
double err2a = 0.01; // acceptable relative error; in some cases the integrals converge too
                     // slowly and the acceptable error is bigger for the sake of performance
double chi = 1e7;
double b = 0.1;
double gamma_p = chi / b;

// calculations with only real, and both real and imaginary parts of the refractive index taken
// into account
enum class Ntype { Re, ReIm };


void print_params(Ntype ntype, double kappa, double theta) {
    double k = kappa / (b * b);
    complex<double> ri = vacuum_refractive_index_perp(b, k);
    double eta = real(ri);
    double nu;
    if (ntype == Ntype::Re) {
        nu = 0;
    } else {
        nu = imag(ri);
    }
    double signed_reverse_tau_parallel = k / (2 * M_PI) * ( theta * theta / 2 + eta /
            norm(ri) - 1 + 1 / (2 * gamma_p * gamma_p) );
    double tau_perp = gamma_p * pow(48 * M_PI / (k * r(gamma_p)), 1/3.0);
    double reverse_tau_d = k * nu;
    // coefficient in front of the integral
    double mult = alpha * b / M_PI * k * k * theta;
    cout << setprecision(11)
         << "gamma_p = " << gamma_p << '\n'
         << "theta = " << theta << '\n'
         << "mult = " << mult << '\n'
         << "reverse_tau_parallel = " << signed_reverse_tau_parallel << '\n'
         << "tau_perp = " << tau_perp << '\n'
         << "reverse_tau_d = " << reverse_tau_d << '\n';
}

double kernel(Ntype ntype, double kappa, double theta) {
    double k = kappa / (b * b);
    complex<double> ri = vacuum_refractive_index_perp(b, k);
    double eta = real(ri);
    double nu;
    if (ntype == Ntype::Re) {
        nu = 0;
    } else {
        nu = imag(ri);
    }
    return classical_synchrotron_power_kernel(eta + 1i * nu, 1, b, gamma_p, theta, k);
}

int main() {
    /**** 1 ****/
    double kappa = 2;

    // (a)
    Ntype ntype = Ntype::Re;
    double theta = 1.6e-3;

    //print_params(ntype, kappa, theta);
    double expected_res = 0.0044542;
    double res = kernel(ntype, kappa, theta);
    bool passed_1a;
    if (fabs(res / expected_res - 1) < err) {
        passed_1a = true;
    } else {
        passed_1a = false;
        cout << "classical_synchrotron_power_kernel test: \x1b[1;31mfailed\x1b[0m\n";
        cout << "acc. error   = " << err          << '\n'
             << "expected res = " << expected_res << '\n'
             << "computed res = " << res          << '\n';
    }

    // (b)
    ntype = Ntype::ReIm;
    theta = 4e-3;

    //print_params(ntype, kappa, theta);
    expected_res = 7.6298e-5;
    res = kernel(ntype, kappa, theta);
    bool passed_1b;
    if (fabs(res / expected_res - 1) < err) {
        passed_1b = true;
    } else {
        passed_1b = false;
        cout << "classical_synchrotron_power_kernel test: \x1b[1;31mfailed\x1b[0m\n";
        cout << "acc. error   = " << err          << '\n'
             << "expected res = " << expected_res << '\n'
             << "computed res = " << res          << '\n';
    }

    /**** 2 ****/
    kappa = 50;

    // (a)
    ntype = Ntype::Re;
    theta = 1e-4;

    //print_params(ntype, kappa, theta);
    expected_res = 8.7862e-6;
    res = kernel(ntype, kappa, theta);
    bool passed_2a;
    if (fabs(res / expected_res - 1) < err2a) {
        passed_2a = true;
    } else {
        passed_2a = false;
        cout << "classical_synchrotron_power_kernel test: \x1b[1;31mfailed\x1b[0m\n";
        cout << "acc. error   = " << err2a        << '\n'
             << "expected res = " << expected_res << '\n'
             << "computed res = " << res          << '\n';
    }

    // (b)
    ntype = Ntype::ReIm;
    theta = 5e-3;

    //print_params(ntype, kappa, theta);
    expected_res = 7.2779e-5;
    res = kernel(ntype, kappa, theta);
    bool passed_2b;
    if (fabs(res / expected_res - 1) < err) {
        passed_2b = true;
    } else {
        passed_2b = false;
        cout << "classical_synchrotron_power_kernel test: \x1b[1;31mfailed\x1b[0m\n";
        cout << "acc. error   = " << err          << '\n'
             << "expected res = " << expected_res << '\n'
             << "computed res = " << res          << '\n';
    }

    /**** 3 ****/
    kappa = 8;

    // (a)
    ntype = Ntype::Re;
    theta = 1.7e-4;

    //print_params(ntype, kappa, theta);
    expected_res = 2.2962e-4;
    res = kernel(ntype, kappa, theta);
    bool passed_3a;
    if (fabs(res / expected_res - 1) < err) {
        passed_3a = true;
    } else {
        passed_3a = false;
        cout << "classical_synchrotron_power_kernel test: \x1b[1;31mfailed\x1b[0m\n";
        cout << "acc. error   = " << err          << '\n'
             << "expected res = " << expected_res << '\n'
             << "computed res = " << res          << '\n';
    }

    // (b)
    ntype = Ntype::ReIm;
    theta = 0.02;

    //print_params(ntype, kappa, theta);
    expected_res = 2.1934e-5;
    res = kernel(ntype, kappa, theta);
    bool passed_3b;
    if (fabs(res / expected_res - 1) < err2a) {
        passed_3b = true;
    } else {
        passed_3b = false;
        cout << "classical_synchrotron_power_kernel test: \x1b[1;31mfailed\x1b[0m\n";
        cout << "acc. error   = " << err2a          << '\n'
             << "expected res = " << expected_res << '\n'
             << "computed res = " << res          << '\n';
    }

    /**** 4 ****/
    kappa = 100'000;

    // (a)
    ntype = Ntype::Re;
    theta = 1.1e-5;

    //print_params(ntype, kappa, theta);
    expected_res = 3.9308e-2;
    res = kernel(ntype, kappa, theta);
    bool passed_4a;
    if (fabs(res / expected_res - 1) < err) {
        passed_4a = true;
    } else {
        passed_4a = false;
        cout << "classical_synchrotron_power_kernel test: \x1b[1;31mfailed\x1b[0m\n";
        cout << "acc. error   = " << err          << '\n'
             << "expected res = " << expected_res << '\n'
             << "computed res = " << res          << '\n';
    }

    // (b)
    ntype = Ntype::ReIm;
    theta = 1.1e-5;

    //print_params(ntype, kappa, theta);
    expected_res = 3.8468e-2;
    res = kernel(ntype, kappa, theta);
    bool passed_4b;
    if (fabs(res / expected_res - 1) < err) {
        passed_4b = true;
    } else {
        passed_4b = false;
        cout << "classical_synchrotron_power_kernel test: \x1b[1;31mfailed\x1b[0m\n";
        cout << "acc. error   = " << err          << '\n'
             << "expected res = " << expected_res << '\n'
             << "computed res = " << res          << '\n';
    }

    /**** ****/

    if (  passed_1a && passed_1b && passed_2a && passed_2b && passed_3a && passed_3b
       && passed_4a && passed_4b ) {
        cout << "classical_synchrotron_power_kernel test: \x1b[32mpassed\x1b[0m\n";
    }

    return 0;
}
