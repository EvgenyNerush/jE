/* This test is quite complicated. It rely on metropolis, histogram_1D, trap_rule and
 * bks_synchrotron_emission_probability functions.
 *
 * Here we check classical asymptotic of Baier-Katkov-Strakhovenko formula for synchrotron emission
 * probability for refractive index equal to 1. We use the Metropolis algorithm here to generate
 * the distribution of the emitted photons in (theta, omega) space. The Metropolis algorithm
 * produces the shape of the distribution only, and one need normalization to compute the overall
 * photon energy from this distribution. To find the normalization, we compute the probability of
 * the photon emission per frequency interval for a certain frequency, namely $ dW / d\omega $ (for
 * some \omega_0), by two different methods: from the photon distribution and integrating
 * numerically $ d^2 W / d\omega d\theta $ over the angle. The results should be the same for both
 * methods that allows to find the right normalization.
 *
 * This test consumes about 6 MB of memory and with Xeon X5550 processor takes about 130 s to
 * execute.
 */

#include <iostream>
#include <vector>
#include <chrono>
#include "../src/histogram.hpp"
#include "../src/metropolis.hpp"
#include "../src/proposal_density.hpp"
#include "../src/radiation.hpp"
#include "../src/rng.hpp"

using namespace std;
using namespace std::chrono;

int main() {
    /* If the magnetic field is weak (and quantum parameter $\chi$ is small), the emitted energy
     * computed with BKS formula is $ I \approx (1 - 6 \chi + 48 \chi^2) I_{cl} $, where $I_{cl}$
     * is the emitted energy computed with classical synchrotron formula. See [V. B. Berestetskii
     * and E.  M. Lifshitz and L. P.  Pitaevskii, Quantum Electrodynamics, Pergamon, New York,
     * 1982].
     */
    double    b        = 1e-5;
    double    gamma_e  = 2e3;
    double    chi      = b * gamma_e;      // hence (1 - 6 chi + 48 chi^2) = (1 - 0.12 + 0.0192)
    double    om_c     = omega_c(gamma_e); // chi << 1, thus om_c is the scale of the spectrum
    double    om_0     = 0.3 * om_c;       // frequency which is used in the normalization
    long long n        = 100'000;          // number of the photons to emit
    size_t    n_in_bin = 1'000;            // number of particles in a bin for histogram_1D
    double    acc_err  = 0.02;             // accepted relative error

    // target distributions
    auto bks_td = bks_synchrotron_td(1, 1, b, gamma_e);
    // proposal_density(pair(rng_state, tuple(theta, omega)))
    auto pd = make_proposal_density< uint64_t, double, double >
                  ( pm_rng
                  , make_tuple(1 / gamma_e, 0.5 * om_c) // width of the proposal density
                  , pm_cast_with_amplitude
                  );

    steady_clock::time_point t1 = steady_clock::now();

    // emitted photons in (theta, omega) space
    vector< tuple<double, double> > xs
        = metropolis< tuple<double, double>, uint64_t, uint64_t >
              ( pm_rng
              , 2468
              , pm_cast_to_01
              , pd
              , make_pair(8642, make_tuple(0, om_c)) // initial point
              , bks_td
              , n
              );

    // photon distribution in (omega) space
    vector<double> xs_omega;
    xs_omega.reserve(xs.size());
    for (size_t i = 0; i < xs.size(); ++i) {
        xs_omega.emplace_back(get<1>(xs[i]));
    }

    auto bin_boundaries = histogram_1D(n_in_bin, xs_omega);
    // slight correction of om_0 to be in the middle between two bin boundaries
    size_t i = 0;
    for ( ; i < bin_boundaries.size() - 1 and bin_boundaries[i] < om_0; ++i ) {
    };
    om_0 = 0.5 * (bin_boundaries[i] + bin_boundaries[i - 1]);
    // non-normalized photon distribution function $ dW/d\omega $ at om_0
    double f_0 = static_cast<double>(n_in_bin) / (bin_boundaries[i] - bin_boundaries[i - 1]);

    // Numerical computation of \int_{-\pi}^\pi (d^2 W / d\omega d\theta) \, d\theta, for omega =
    // om_0
    long long ntheta = 100;
    double theta_max = 5 / gamma_e;
    double dtheta    = 2 * theta_max / static_cast<double>(ntheta - 1);
    auto theta_nodes = ranges::v3::iota_view(0, ntheta)
                     | ranges::v3::views::transform(
                           [=](long long i){ return -theta_max + dtheta * static_cast<double>(i); }
                       );
    // d^2 W / d\theta d\omega at omega = om_0
    auto g = function<double(double)>(
         [=](double theta) { return bks_td(tuple(theta, om_0)); }
    );
    // photon distribution function $ dW/d\omega $ at om_0, for a single full-circle path of an
    // electron
    double f_fc = trap_rule(g, theta_nodes);

    // Computation of the overall photon energy
    double norm = f_fc / f_0;
    double bks_I = 0;
    for (auto omega: xs_omega) {
        bks_I += b * omega * norm; // in t_{rf} normalization, $\hbar \omega$ is $b \omega$.
    }

    steady_clock::time_point t2 = steady_clock::now();

    // I_cl / 2 \pi r = 2 e^4 B^2 \gamma_e^2 / 3 m^2
    double cl_I = 4 * M_PI / 3 * alpha * b * r(gamma_e) * pow(gamma_e, 2);
    double corrected_I = cl_I * (1 - 6 * chi + 48 * chi * chi);

    if (fabs(bks_I / corrected_I - 1) < acc_err) {
        cout << "bks_low_chi test: \x1b[32mpassed\x1b[0m\n";
    } else {
        cout << "bks_low_chi test: \x1b[1;31mfailed\x1b[0m\n";
        cout << "acc. error  = " << acc_err                       << '\n'
             << "rel. error  = " << fabs(bks_I / corrected_I - 1) << '\n'
             << "classical I = " << cl_I                          << '\n'
             << "corrected I = " << corrected_I                   << '\n'
             << "BKS I       = " << bks_I                         << '\n';
    }

    duration<double> time_span = duration_cast<duration<double>>(t2 - t1);
    cout << "test takes " << time_span.count() << " seconds\n";

    return 0;
}
