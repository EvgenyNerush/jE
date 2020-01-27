/* Here we compare results of bks_synchrotron_emission_probability with the results computed
 * with WolframCloud, in the case there Cherenkov synchronism is important, i.e. for the phase
 * proportional to -t/tau_parallel + t^3 / tau_perp^3. See also
 * https://www.wolframcloud.com/obj/4df105db-f7ae-4f5c-b700-1661857f91d7
 */
#include <vector>
#include <iostream>
#include <iomanip>
#include "../src/radiation.hpp"

using namespace std;

int main(){
    double acceptable_relative_error = 0.03;
    double omega_s_fraction = 0.5; // \hbar \omega' / mc^2 \gamma_e
    double chi = omega_s_fraction / M_PI * sqrt(3 / 16.0); // ...thus tau_parallel = tau_perp,
                                                           // for theta = 0 and delta_ri = 0
    double gamma_e = 1e4;
    double b       = chi / gamma_e;
    double omega   = gamma_e / b * omega_s_fraction / (1 + omega_s_fraction);
    double omega_s = gamma_e * omega_s_fraction / b;

    // ri1, ri2, etc. lead to the chenge of the signum of the linear term in the phase
    double ri0 = 1 + 0.1  / pow(gamma_e, 2); // => tau_perp / tau_parallel = 0.8, varsigma > 0
    double ri1 = 1 + 0.5  / pow(gamma_e, 2); // => tau_perp / tau_parallel =   0, varsigma < 0
    double ri2 = 1 + 0.75 / pow(gamma_e, 2); // => tau_perp / tau_parallel = 0.5, varsigma < 0
    double ri3 = 1 + 2    / pow(gamma_e, 2); // => tau_perp / tau_parallel =   3, ...
    double ri4 = 1 + 6.5  / pow(gamma_e, 2); // => tau_perp / tau_parallel =  12
    double ri5 = 1 + 24.5 / pow(gamma_e, 2); // => tau_perp / tau_parallel =  48

    double w0 = bks_synchrotron_emission_probability(ri0, 1, b, gamma_e, 0, omega);
    double w1 = bks_synchrotron_emission_probability(ri1, 1, b, gamma_e, 0, omega);
    double w2 = bks_synchrotron_emission_probability(ri2, 1, b, gamma_e, 0, omega);
    double w3 = bks_synchrotron_emission_probability(ri3, 1, b, gamma_e, 0, omega);
    double w4 = bks_synchrotron_emission_probability(ri4, 1, b, gamma_e, 0, omega);
    double w5 = bks_synchrotron_emission_probability(ri5, 1, b, gamma_e, 0, omega);

    double tau_perp = gamma_e * pow(12 * M_PI / (omega_s * gamma_e), 1/3.0);
    // coefficients a1 and a2 arises in W_m 
    double a1 = (1 + pow((gamma_e - omega * b) / gamma_e, 2))
              * pow(tau_perp * tau_perp / gamma_e, 2);
    double a2 = pow(omega * b * b * tau_perp / (gamma_e * gamma_e), 2);
    // f is a function that computes value of W_m from values of the integrals I and J, see ref. to
    // wolframcloud above; actually for chi used here the spin term which is proportional to I^2,
    // is negligible in comparison with the main term proportional to J^2 See also
    // bks_emission_probability from ../src/radiation.hpp, note that coefficients c3 / tau_perp and
    // c1 * gamma_e / (tau_perp^2) there should coincide with I and J, respectively
    auto f = [=](double I, double J) {
                 return alpha * M_PI / (8 * omega)
                      * pow(omega_s / omega, 2) * (a1 * J * J + a2 * I * I);
             };

    double f0 = f(0.0975, 0.0549);
    double f1 = f(0.824 , 0.230 );
    double f2 = f(1.25  , 0.0850);
    double f3 = f(0.574 , 0.572 );
    double f4 = f(0.408 , 0.816 );
    double f5 = f(0.201 , 0.810 ) * 2;

    vector<double> ws({w0, w1, w2, w3, w4, w5});
    vector<double> fs({f0, f1, f2, f3, f4, f5});
    vector<double> errs;
    for (int i = 0; i < 5; ++i) {
        errs.push_back(fabs(ws[i] - fs[i]) / fs[i]);
    }
    bool is_ok = true;
    for (int i = 0; i < 5; ++i) {
        is_ok = is_ok and errs[i] < acceptable_relative_error;
    }

    if (!is_ok) {
        cout << "synchrotron_Cherenkov_radiation test: \x1b[1;31mfailed\x1b[0m\n";
        cout << setw(12) << left << "i           w num       w theor     rel err\n";
        for (int i = 0; i < 5; ++i) {
            cout << setw(12) << left << i
                 << setw(12) << left << ws[i]
                 << setw(12) << left << fs[i]
                 << setw(12) << left << errs[i] << '\n';
        }
    } else {
        cout << "synchrotron_Cherenkov_radiation test: \x1b[32mpassed\x1b[0m\n";
    }

    return 0;
}
