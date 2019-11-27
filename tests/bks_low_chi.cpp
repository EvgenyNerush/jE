/* This test is quite complicated. It rely on metropolis, histogram_1D, trap_rule and
 * bks_synchrotron_emission_probability functions.
 *
 * Here we check classical and quantum asymptotics of Baier-Katkov-Strakhovenko formula for
 * synchrotron emission probability for refractive index equal to 1. We use the Metropolis
 * algorithm here to generate the distribution of the emitted photons in (theta, omega) space. The
 * Metropolis algorithm produces the shape of the distribution only, and one need normalization to
 * compute the overall photon energy from this distribution. To find the normalization, we compute
 * the probability of photon emission per frequency interval for a certain frequency, namely $ dW /
 * d\omega(some \omega_0) $, by two different methods: from the photon distribution and integrating
 * numerically $ d^2 W / d\omega d\theta $ over the angle. The results should be the same for both
 * methods that allows to find the right normalization.
 *
 * Emission probability per frequency interval and per theta interval can be obtained analogous to
 * jackson1483_num from ../src/radiation.hpp (see comments there). Note that the generated photons
 * are evenly distributed around the electron trajectory (which is
 * a circle), i.e. they are evenly distributed in the angle $ \phi $ which corresponds to the direction
 * of the photon wavevector in the plane of the electron trajectory. Thus, the solid angle is $
 * d\Omega = \cos \theta \, d\phi d\theta $, and for the emission probability we get $ d^2 W /
 * d\theta d\omega = 2 W_m \omega^2 / \pi^2 $ with $ W_m $ the emission probability computed with
 * bks_synchrotron_emission_probability function from ../src/radiation.hpp.
 *
 * This test takes about one minute (!) with Xeon X5550 processor.
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
    double    acc_err  = 0.01;             // accepted relative error
    long long n        = 1'000'000;        // number of the photons to emit
    size_t    n_in_bin = 10'000;           // number of particles in a bin for histogram_1D
    double    om_0     = 0.5 * om_c;       // frequency which is used in the normalization
    double    om_guard = 0.0 * om_c;       // photons with omega < om_guard are not generated;
                                           // their contribution to the overall energy is small,
                                           // but they take a lot of time to be generated because
                                           // the emission probability is singular at omega = 0

    // target distributions
    auto bks_td = function<double(tuple<double, double>)>(
        [=](tuple<double, double> theta_omega) {
            double theta = get<0>(theta_omega);
            double omega = get<1>(theta_omega);
            if (omega > om_guard and omega * b < gamma_e) {
                return 2 * pow(omega / M_PI, 2) * cos(theta)
                         * bks_synchrotron_emission_probability(1, b, gamma_e, theta, omega);
            } else {
                return 0.0;
            }
        }
    );
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
              , make_pair(8642, make_tuple(0, 2 * om_guard)) // initial point
              , bks_td
              , n
              );

    vector<double> xs_omega; // photon distribution in (omega) space
    xs_omega.reserve(xs.size());
    for (size_t i = 0; i < xs.size(); ++i) {
        xs_omega.emplace_back(get<1>(xs[i]));
    }
    auto hist = histogram_1D(n_in_bin, xs_omega);

    // slight correction of om_0 to be in the middle between two hist bin boundaries
    size_t i = 0;
    for ( ; i < hist.size() - 1 and hist[i] < om_0; ++i ) { };
    om_0 = 0.5 * (hist[i] + hist[i - 1]);
    // non-normalized photon distribution function $ dW/d\omega $ at om_0
    double f_0 = static_cast<double>(n_in_bin) / (hist[i] - hist[i - 1]);

    // Numerical computation of \int_{-\pi}^\pi (d^2 W / d\omega d\theta) \, d\theta, for omega =
    // om_0
    long long ntheta = 100;
    double theta_max = 5 / gamma_e;
    double dtheta    = 2 * theta_max / static_cast<double>(ntheta - 1);
    auto theta_nodes = ranges::v3::iota_view(0, ntheta)
                     | ranges::v3::views::transform(
                           [=](long long i){ return -theta_max + dtheta * static_cast<double>(i); }
                       );
    auto g = function<double(double)>(
         [=](double theta) { return bks_td(tuple(theta, om_0)); }
    );
    double dWdomega_0 = trap_rule(g, theta_nodes);

    // Computation of the overall photon energy
    double norm = dWdomega_0 / f_0;
    double bks_I = 0;
    for (auto x: xs) {
        double omega = get<1>(x);
        bks_I += b * omega * norm;
    }
    // I_cl / 2 \pi r = 2 e^4 B^2 \gamma_e^2 / 3 m^2
    double cl_I = 4 * M_PI / 3 * alpha * b * r(gamma_e) * pow(gamma_e, 2);
    double corrected_I = cl_I * (1 - 6 * chi + 48 * chi * chi);

    steady_clock::time_point t2 = steady_clock::now();
    duration<double> time_span = duration_cast<duration<double>>(t2 - t1);
    cout << "test takes " << time_span.count() << " seconds\n";

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

    return 0;
}
