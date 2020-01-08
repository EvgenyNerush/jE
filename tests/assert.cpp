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
      // electrodynamics, Princeton U. preprint DOE ER, 1986]. Note that McDonald uses chi which is
      // a half of chi we use.
        double precision = 0.1;
        auto n = [](double b, double omega) {
            return 4 * M_PI / (alpha * b * b) * (vacuum_refractive_index(b, omega) - 1);
        };
        assert(fabs(n(0.1, 0.1) * 45 / 14 - 1)                 < precision);
        assert(fabs(n(1, 70) / (-0.278 * pow(70, -4/3.0)) - 1) < precision);
        assert(fabs(n(1, 0.5) / 0.35 - 1)                      < precision);
        auto f = [=](double chi) { return n(1, chi); };
        assert(fabs(bisection(f, 0, 40, 10).value() / 15 - 1)  < precision);
    }

    { // Test of zipWith from proposal_density.hpp
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

    { // Classical synchrotron emission probability and BKS emission probability should be the same
      // if the quantum parameter $ \chi $ is small, $ chi << 1 $.
        double b       = 1e-6;
        double gamma_e = 1e3;
        double chi     = b * gamma_e;
        double om_c    = omega_c(gamma_e);
        double w_cl = synchrotron_emission_probability(b, gamma_e,  true, 0, om_c)
                    + synchrotron_emission_probability(b, gamma_e, false, 0, om_c);
        double w_bks = bks_synchrotron_emission_probability(1, 1, b, gamma_e, 0, om_c);
        assert(fabs(w_bks / w_cl - 1) < 10 * chi);
    }

    {
     {// In classical limit, if one scales (say, makes it $ m $ times larger) simultaneously the
      // mass of the emitting particle (hence the curvature radius) and the wavelength of the
      // emitted photon, the emission process would be similar to the non-scaled case in the sense
      // in which two triangles can be similar. Namely, both "transverse" and "longitudinal" time
      // scales of the emission process in the second case $ m $ times larger than in the first
      // one. This allows to find the relation which is tested here: $ W(m, omega / m) = m^3 W(1,
      // omega) $, where $ W $ is the emission probability with first argument the particle mass
      // and the second the photon frequency.
        double b       = 1e-6;
        double gamma_e = 1e2;
        // just for the record, theta = 0 and omega = 0.42 omega_c is approximate point of the
        // maximum of the synchrotron energy spectrum
        double omega   = 0.42 * omega_c(gamma_e);
        double theta   = 0.5 / gamma_e;
        double m       = 10;
        double w0 = bks_synchrotron_emission_probability(1, m, b, gamma_e, theta, omega / m);
        double w1 = bks_synchrotron_emission_probability(1, 1, b, gamma_e, theta, omega)
                  * pow(m, 3);
        assert(fabs(w0 / w1 - 1) < 1e-3);
     }
     {// The trick described above does not hold for the quantum case (chi >~ 1) mostly because the
      // rising spin term in the emission probability. However, if the particle mass and the photon
      // frequency are multiplied by $ m $ (that preserves the recoil effect), and the magnetic
      // field by $ m^2 $ (that preserves the ratio of the trajectory curvature radius and the
      // photon wavelength), the emission probability depend on $ m $ as $ m^{-3} $ (see equations
      // in the comments to bks_emission_probability from ../src/radiation.hpp).
        double b       = 1e-2;
        double gamma_e = 3e2;
        double omega   = 0.5 * gamma_e / b;
        double theta   = 0.5 / gamma_e;
        double m = 10;
        // Note here that as we chenge the magnetic field, and the time and frequency normalization
        // depend on it, the scaling described above looks different in the normalized variables.
        // In particular, one should take into account that bks_synchrotron_emission_probability
        // returns the product of the virtual volume and the emission probability, and in the
        // normalized units the virtual volume scales as $ m^6 $.
        double w0
            = bks_synchrotron_emission_probability(1, m, m * m * b, gamma_e, theta, omega / m);
        double w1
            = bks_synchrotron_emission_probability(1, 1,         b, gamma_e, theta, omega)
            * pow(m, 3);
        assert(fabs(w0 / w1 - 1) < 1e-3);
     }
     {// Furthermore, one can slightly increase the refractive index, and compensate the resulting
      // dephasing by the decrease of the electron Lorentz factor. To conserve the curvature
      // radius, the particle mass should be respectively increased. Thus, after these
      // manipulations the emission probability should be unchanged.
      double b        = 1e-2;
      double gamma_0  = 1e2; // 1 - v \approx 1 / 2 gamma_e^2
      double omega    = 0.2 * gamma_0 / b;
      double theta    = 0.5 / gamma_0;
      double delta_ri = 2e-5; // excess of the refractive index over the unity

      double ri       = 1 + delta_ri;
      double gamma_1  = 1 / sqrt(2 * delta_ri + 1 / (gamma_0 * gamma_0)); // 84.5
      double m        = gamma_0 / gamma_1;

      double w0 = bks_synchrotron_emission_probability( 1, 1, b, gamma_0, theta, omega);
      double w_ = bks_synchrotron_emission_probability(ri, 1, b, gamma_0, theta, omega);
      double w1 = bks_synchrotron_emission_probability(ri, m, b, gamma_1, theta, omega);
      // w_ and w0 differs on about several percents, whereas w1 should be equal to w0
      assert(fabs(w1 - w0) < 1e-3 * (w_ - w0));
     }
    }

    { // Here we compare results of bks_synchrotron_emission_probability in "Cherenkov zone" with
      // the results obtained with Wolfram Cloud, see
      // https://www.wolframcloud.com/obj/4df105db-f7ae-4f5c-b700-1661857f91d7
      double eps_s_fraction = 0.5; // \hbar \omega' / mc^2 \gamma_e
      double chi = sqrt(12 * M_PI) * eps_s_fraction; // ...thus tau_parallel = tau_perp, for
                                                     // theta = 0 and delta_ri = 0
      double gamma_e = 1e4;
      double b = chi / gamma_e;
      double omega = gamma_e / b * eps_s_fraction / (1 + eps_s_fraction);

      // ri2 and ri3 lead to the chenge of the signum of the linear term in the phase; ri2 yields
      // the same value of tau_parallel as ri = 1, and ri3 yields 1/9 of that value for
      // tau_parallel
      double ri2 = 1 + 1 / pow(gamma_e, 2);
      double ri3 = 1 + 5 / pow(gamma_e, 2);
      double w1 = bks_synchrotron_emission_probability(  1, 1, b, gamma_e, 0, omega);
      double w2 = bks_synchrotron_emission_probability(ri2, 1, b, gamma_e, 0, omega);
      double w3 = bks_synchrotron_emission_probability(ri3, 1, b, gamma_e, 0, omega);

      double a1
          = 0.5 * (1 + pow(1 / eps_s_fraction, 2)); // (\eps^2 + \eps'^2) / 2 \eps^2, see
                                                    // equation for W_m in bks_emission_probability
                                                    // from ../src/radiation.hpp
      double a2 = 0.5 * pow(omega * b / gamma_e, 2);
      // function that computes value \propto W_m from values of the integrals I and J, see ref. to
      // wolframcloud above
      auto f = [=](double I, double J) { return a1 * J * J + a2 * I * I; };

      cout << w1 / w2 << '\t' << f(0.0464, 0.0236) / f(0.0373, 0.621) << '\n'
           << w1 / w3 << '\t' << f(0.0464, 0.0236) / f(0.0670, 1.07) << '\n';
    }

    cout << "assertions: \x1b[32mpassed\x1b[0m\n";
    return 0;
}
