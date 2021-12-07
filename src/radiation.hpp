/**
 * @file
 * @brief Formulas for the photon emission by an ultrarelativistic electron.
 * @details See \ref radiation and \ref synchrotron_radiation for further details.
 */

#include <complex>
#include <cmath>
#include <optional>
#include <range/v3/view/drop.hpp>
#include <range/v3/view/iota.hpp>
#include <range/v3/view/reverse.hpp>
#include <range/v3/view/transform.hpp>

#include <iostream> // qwe

using namespace std::literals::complex_literals;

/**
 * @defgroup radiation Radiation module
 * @brief General formulas for the photon emission by an ultrarelativistic electron.
 * @details Formulas describing the emission of photons by an ultrarelativistic electron in the
 * general case. I.e., the emission can occur either in the synchrotron or in the dipole regime.
 * However, we use here normalization natural for the synchrotron emission. Thus, some value of the
 * magnetic field @f$ B @f$ is used as follows. In functions provided by this module time is
 * normalized to the radiation formation time 
 * @f[
 * t_{rf} = m c / (e B),
 * @f]
 * coordinates are normalized to @f$ c t_{rf} @f$ and energy to the electron
 * rest energy @f$ m c^2 @f$, i.e. @f$ m @f$ the electron mass. Velocity is
 * normalized to the speed of light @f$ c @f$ and frequency is normalized to
 * @f$ 1 / t_{rf} @f$.  Note that the photon of the cyclic frequency @f$ \omega
 * @f$ in the normalized units has energy @f$ \omega \, b @f$ with @f$ b = B /
 * B_S @f$, and @f$ B_S = m^2 c^3 / e \hbar @f$ the Sauter--Schwinger field
 * (the critical field of QED). For the synchrotrom radiation formulas, see
 * @ref synchrotron_radiation.
 * @{
 */

/**
 * The fine structure constant, CODATA 2018 recommended value
 */
const double alpha = 1 / 137.035999084;

/**
 * <a href = "https://en.wikipedia.org/wiki/Bisection_method">Bisection method</a> to find the root
 * of the equation @f$ f(x) = 0 @f$ on the interval @f$ [x_a, x_b] @f$.
 * @param  f  @f$ f @f$, a continuous function
 * @param  xa @f$ x_a @f$
 * @param  xb @f$ x_b @f$
 * @param  n  the number of the iterations; the method error is less than @f$ (x_b - x_a) / 2^n @f$
 * @return    <tt> std::nullopt </tt> if the root is out of @f$ [x_a, x_b] @f$, (<tt> std::optional
 *            </tt> of) @f$ x_a @f$ if @f$ f(x_a) = 0 @f$, @f$ x_b @f$ if @f$ f(x_b) = 0 @f$, @f$
 *            (x_a + x_b) / 2 @f$ if @p n < 1, and approximate root value otherwise
 */
std::optional<double> bisection( std::function<double(const double&)> f
                               , double xa
                               , double xb
                               , int    n
                               ) {
    double f1 = f(xa);
    double f2 = f(xb);
    if (f1 * f2 > 0) {
        return std::nullopt;
    } else if (f1 == 0) {
        return std::optional(xa);
    } else if (f2 == 0) {
        return std::optional(xb);
    } else { // f1 * f2 < 0
        double x1 = xa;
        double x2 = xb;
        double x = 0.5 * (x1 + x2);
        double v = f(x);
        for (int i = 0; i < (n - 1) and v != 0; ++i) {
            if (v * f1 < 0) {
                x2 = x;
                f2 = v;
            } else { // v * f2 < 0
                x1 = x;
                f1 = v;
            }
            x = 0.5 * (x1 + x2);
            v = f(x);
        }
        return std::optional(x);
    }
}

/**
 * Trapezoidal rule of integration.
 * @param f       function to be integrated
 * @param t_nodes nodes in which the values of @p f are computed; this function can take @p t_nodes
 *                as @p range or @p view (from C++20 ranges) or as @p std::vector, with doubles
 *                within
 */
template <typename R, typename V>
R trap_rule(std::function<R(double)> f, V t_nodes) {
    R acc(0);
    double t_prev = *(t_nodes.begin());
    R      f_prev = f(t_prev);
    for (auto t: t_nodes | ranges::v3::views::drop(1)) {
        R f_curr = f(t);
        acc += 0.5 * (t - t_prev) * (f_curr + f_prev);
        t_prev = t;
        f_prev = f_curr;
    }
    return acc;
}

/**
 * The mode amplitude @f$ C_m @f$ in the emitted field (for a unit volume @f$ V = 1 @f$), computed
 * numerically, as follows.
 *
 * The emitted field can be decomposed by modes of virtual superconductive rectangular box of size
 * @f$ L_x \times L_y \times L_z @f$:
 * @f[
 *     \mathbf{E} = \sum_m C_m \mathbf{E}_m, \quad
 *     \mathbf{B} = \sum_m C_m \mathbf{B}_m,
 * @f]
 * with well-known sine/cosine spatial and @f$ e^{-i \omega_m t} @f$ temporal structure of complex
 * rectangular resonator modes. The modes are orthogonal
 * @f[
 *     \frac{1}{8 \pi} \int_V \left( \mathbf{E}_m \mathbf{E}_l^* + \mathbf{B}_m \mathbf{B}_l^*
 *     \right) \, dV = \hbar \omega_m \delta_{ml},
 * @f]
 * thus, the energy of the emitted field is
 * @f[
 *     I = \frac{1}{8 \pi} \int_V \left( \mathbf{EE}^* + \mathbf{BB}^* \right) \, dV = \sum_m \hbar
 *     \omega_m |C_m|^2,
 * @f]
 * and @f$ |C_m|^2 @f$ can be interpreted as the emission probability of the photon of mode @f$ m
 * @f$. It can be found from Maxwell's equations []:
 * @f[
 *     C_m = -\frac{1}{2 \hbar \omega_m} \int_t \int_V \mathbf{j E}_m^* \, dV \, dt,
 * @f]
 * where @f$ \mathbf{j} = -e \mathbf{v} \delta(\mathbf{r} - \mathbf{r}(t)) @f$ is the current of
 * the emitting electron which position is @f$ \mathbf{r}(t) @f$.
 *
 * For an ultrarelativistic particle, which emits mostly in the forward direction, the computation
 * of @f$ C_m @f$ can be further simplified. Each of the complex modes is formed by eight complex
 * plane waves @f$ \sim \exp(-i \omega_m t + i \mathbf{k}_m \mathbf{r}) @f$ (exept a few modes with
 * wave vector parallel to the box boundaries). This yields eight terms in the integral over @f$ t
 * @f$ in the expression for @f$ C_m @f$. One can note that one of the terms oscillates much slower
 * than the others which hence can be dropped (e.g., @f$ \exp( i \omega t - i k_x x(t)) @f$ cannot
 * be dropped, whereas @f$ \exp( i \omega t + i k_x x(t)) @f$ can be, if @f$ k_x \approx \omega_m /
 * c @f$ and @f$ x(t) \approx c t @f$). The remaining term corresponds to the wave polarization
 * direction @f$ \mathbf{p_m} @f$ (with @f$ p_m^2 = 1 @f$). The amplitude of the remaining wave,
 * @f$ a_m @f$, can be found from the normalization: the wave energy is @f$ \hbar \omega_m / 8 @f$
 * hence @f$ a_m = (2 \pi \hbar \omega_m / V)^{1/2} @f$.  Therefore,
 * @f[
 *     C_m = - \frac{e}{2} \sqrt{ \frac{\pi}{\hbar \omega_m V} } \int_t \mathbf{v p}_m \exp \left[
 *     i \omega_m t - i \mathbf{k}_m \mathbf{r}(t) \right] \, dt.
 * @f]
 * This formula is used here to compute @f$ C_m @f$ with #trap_rule. To compute the emission
 * probability, use <tt> <a href = "https://en.cppreference.com/w/cpp/numeric/complex/norm"> norm
 * </a>(c_m) </tt>.
 * @param vp      the product of the mode polarization direction @f$ \mathbf{p} @f$ and the
 *                electron velocity @f$ \mathbf{v} @f$, as the function of @f$ t @f$.
 * @param nr      @f$ \mathbf{ n \cdot r }(t) @f$, with @f$ \mathbf{n} = c \mathbf{k}_m / \omega_m
 *                @f$ and @f$ \mathbf{r}(t) @f$ the electron position
 * @param t_nodes values of @f$ t_0, t_1, t_2, ... @f$ used for the integral computation; this
 *                function can take @p t_nodes as @p range or @p view (from C++20 ranges) or as @p
 *                std::vector, with doubles within
 * @param omega   > 0, normalized frequency of the emitted photon, @f$ \omega_m @f$
 */
template <typename V>
std::complex<double> c_m( std::function<double(double)> vp
                        , std::function<double(double)> nr
                        , V                             t_nodes
                        , double                        omega
                   ) {
    auto f = std::function<std::complex<double>(double)>(
        [=](double t) {
            return vp(t) *  exp( 1i * omega * ( t - nr(t) ) );
        }
    );
    return sqrt(alpha * M_PI / (4 * omega)) * trap_rule(f, t_nodes);
}

/**
 * @}
 */

/**
 * @defgroup synchrotron_radiation Synchrotron radiation module
 * @brief Formulas for classical synchrotron radiation.
 * @details If the electron trajectory is bended wider than the angle of @f$ 1 / \gamma_e @f$, the
 * photon emission occurs in the synchrotron regime.  To describe this, the same normalization as
 * in @ref radiation is used here, and it is assumed that the electron trajectory is a circle which
 * in the normalized units is
 * @f[
 * x = r \sin( t / \gamma_e ), \\
 * @f]
 * @f[
 * y = r [\cos( t / \gamma_e ) - 1],
 * @f]
 * with @f$ r = \sqrt{\gamma_e^2 - 1} @f$ the curvature radius and @f$ \gamma_e @f$ the electron
 * Lorentz factor.
 * @{
 */

/**
 * The normalized curvature radius of the electron trajectory
 * @param gamma_e the electron Lorentz factor
 */
double r(double gamma_e) {
    return sqrt(gamma_e * gamma_e - 1);
}

/**
 * The normalized critical frequency of the classical synchrotron emission, see Eq. (14.85) from
 * [Jackson J.D., Classical Electrodynamics, Wiley, 1962].
 * @param gamma_e the electron Lorentz factor
 */
double omega_c(double gamma_e) {
    return 3 * pow(gamma_e, 3) / r(gamma_e);
}

/**
 * The (normalized) energy radiated in unit frequency interval per unit solid angle, as defined in
 * Eq. (14.83) from [Jackson J.D., Classical Electrodynamics, Wiley, 1962].
 * @param b       the normalized magnetic field strength
 * @param gamma_e the electron Lorentz factor
 * @param theta   angle in the plane perpendicular to the normal vector of the trajectory; @p theta
 *                = 0 is the direction of the tangent to the trajectory
 * @param omega   > 0, normalized frequency of the emitted photon
 */
double jackson1483(double b, double gamma_e, double theta, double omega) {
    double a = 1 / (gamma_e * gamma_e) + theta * theta;
    double xi = omega * r(gamma_e) / 3 * pow(a, 3.0/2);
    double k13 = std::cyl_bessel_k(1.0/3, xi);
    double k23 = std::cyl_bessel_k(2.0/3, xi);
    return alpha * b / (3 * M_PI * M_PI)
        * pow(omega * r(gamma_e) * a, 2) * (k23 * k23 + pow(theta * k13, 2) / a);
}

/**
 * Probability of the photon emission in classical synchrotron process, computed numerically with
 * #c_m. The integration parameters (timestep etc.) are set to yield the spectrum #jackson1483_num
 * whose relative difference from the theoretical one is less than 1% up to frequences/angles where
 * photon emission is negligibly small.
 * @param b       the magnetic field strength normalized to the Sauter-Schwinger field
 * @param gamma_e the electron Lorentz factor
 * @param p       the polarization of the emitted wave: <tt> p = true </tt> means polarization
 *                parallel to the normal vector of the trajectory, and <tt> p = false </tt> means
 *                polarization in the plane perpendicular to the normal vector
 * @param theta   angle in the plane perpendicular to the normal vector of the trajectory; @p theta
 *                = 0 is the direction of the tangent to the trajectory
 * @param omega   > 0, normalized frequency of the emitted photon
 */
double synchrotron_emission_probability( double b
                                       , double gamma_e
                                       , bool   p
                                       , double theta
                                       , double omega) {
    // the exponential in the integral is approximately equal to exp(i phi(t))
    double tau_parallel = 4 * M_PI / (omega * (theta * theta + 1 / (gamma_e * gamma_e)));
    double tau_perp     = gamma_e * pow(12 * M_PI / (omega * r(gamma_e)), 1/3.0);

    // we integrate approximately on the interval of +-l scales of the descent of the integrals
    double l = 3;
    double period_fraction = 1/2.0; // fraction of the minimal period which determines the timestep

    // a point where linear and cubic terms in the phase yield the same oscillation period; if
    // varsigma = -1, then d\phi / dt = 0 at t = +-t_s
    double ts = sqrt(pow(tau_perp, 3) / (tau_parallel * 3));
    // tb is the upper limit of the integration
    double tb = l * std::max(tau_perp, ts);

    // the estimate of the period of the exponent oscillations
    // (T = 2 \pi / (d\phi/dt)) at t = tb
    double osc_period = 1 / (1 / tau_parallel + 3 * pow(tb / tau_perp, 2) / tau_perp);
    // lower limit of the integration
    double ta = -tb;
    // number of points for the exponent integration
    long long int nt = llround((tb - ta) / (osc_period * period_fraction));
    // step of the integration
    double dt = (tb - ta) / static_cast<double>(nt - 1);
    auto t_nodes = ranges::v3::iota_view(0, nt)
                 | ranges::v3::views::transform(
                       [=](long long i){ return ta + dt * static_cast<double>(i); }
                   );

    // for c_m function
    auto nr = std::function<double(double)>(
        [=](double t) {
            return r(gamma_e) * sin(t / gamma_e) * cos(theta);
        }
    );

    std::complex<double> c;
    if (p == true) {
        auto vp = std::function<double(double)>(
            [=](double t) {
                return r(gamma_e) / gamma_e * sin(t / gamma_e)
                     * 0.25 * (1 - tanh(8 * (t / tb - 0.7)))
                     *        (1 + tanh(8 * (t / tb + 0.7)));
                     // note that artificial attenuation is added here; the idea behind is that
                     // with this attenuation neighboring bumps quench each other earlier hence
                     // smaller integration period can be used
            }
        );
        c = c_m(vp, nr, t_nodes, omega);
    } else {
        auto vp = std::function<double(double)>(
            [=](double t) {
                return r(gamma_e) / gamma_e * sin(theta) * cos(t / gamma_e)
                     * 0.25 * (1 - tanh(8 * (t / tb - 0.7)))
                     *        (1 + tanh(8 * (t / tb + 0.7)));
                     // note that artificial attenuation is added here; the idea behind is that
                     // with this attenuation neighboring bumps quench each other earlier hence
                     // smaller integration period can be used
            }
        );
        c = c_m(vp, nr, t_nodes, omega);
    }

    return norm(c);
}

/**
 * The (normalized) energy radiated per unit frequency interval per unit solid angle, computed
 * numerically from #synchrotron_emission_probability.
 *
 * An ultrarelativistic particle emit photons in a narrow cone around the direction of the particle
 * velocity. Thus the energy radiated in a certain direction can be readily computed from the
 * energy of the modes of the superconductive virtual box (see #c_m). The density of the modes
 * which has a plane-wave component in some certain unit solid angle and unit frequency interval
 * can be easily found from the boundary conditions, that gives
 * @f[
 *     I = \sum_m \hbar \omega_m |C_m|^2
 *     = \frac{\hbar V }{\pi^3 c^3} \int \int \omega^3 \sum_\mathbf{p} |C|^2 \, d\omega \, d\Omega
 *     = \frac{e^2}{4 \pi^2 c^3} \int \int \omega^2 \sum_\mathbf{p} \left| \int \mathbf{v}(t)
 *     \mathbf{p} \exp \left[ i \omega t - i \mathbf{k r}(t) \right] \, dt \right|^2 \, d\omega \,
 *     d\Omega,
 * @f]
 * with @f$ \mathbf{p} @f$ one of two polarization directions.  This formula is used here in the
 * computations. Note that it coincides with Eq. (14.78) from [Jackson J.D., Classical
 * Electrodynamics, Wiley, 1962].
 * @param b       the normalized magnetic field strength
 * @param gamma_e the electron Lorentz factor
 * @param theta   angle in the plane perpendicular to the normal vector of the trajectory; @p theta
 *                = 0 is the direction of the tangent to the trajectory
 * @param omega   > 0, normalized frequency of the emitted photon
 */
double jackson1483_num(double b, double gamma_e, double theta, double omega) {
    return pow(omega / M_PI, 3) * b
        * ( synchrotron_emission_probability(b, gamma_e, true,  theta, omega)
          + synchrotron_emission_probability(b, gamma_e, false, theta, omega));
}

/**
 * @}
 */

/**
 * @defgroup bks_radiation Quasiclassical radiation module
 * @brief General formulas from BKS theory with vacuum polarization taken into
 * account. The same normalization as in @ref radiation is used here
 * @{
 */

/**
 * Unity <tt>std::function</tt> of any number of arguments.
 * @return 1
 */
template<typename... Types>
auto unity = std::function<double(Types...)>( [](Types... _) { return 1; } );

/**
 * Photon emission probability by an electron, @f$ W_m @f$, found in the framework of
 * Baier-Katkov-Strakhovenko (BKS) quasiclassical theory, namely a sum of the emission
 * probabilities for both photon polarizations, which are then averaged over the polarizations of
 * the emitting electron. Here @f$ m @f$ designates the generalized mode number, and unit virtual
 * box volume is assumed (@f$ V = 1 @f$, see #c_m for details).
 *
 *
 * BKS theory gives the following expression for the energy radiated in unit frequency interval per
 * unit solid angle, see Eq. () from [V.N. Baier and V.M. Katkov and V.M. Strakhovenko,
 * Electromagnetic processes at high energies in oriented single crystals, World Scientific,
 * Singapore, 1998], or see also Eq. (6) in the supplementary material of <a href =
 * "https://doi.org/10.1038/s41467-018-03165-4">this paper</a>:
 *
 * @f[
 *     \frac{d^2 I}{d\omega d\Omega} = \frac{e^2}{4 \pi^2} \left\{
 *         \frac{ \varepsilon^2 + \varepsilon'^2 }{ 2 \varepsilon^2 }
 *         \left| \int dt\, \frac{ \mathbf{n \times [(n - v) \times \dot v]} }
 *                               { (1 - \mathbf{n v})^2 }
 *                \exp [ i \omega' ( t - \mathbf{n r} ) ]
 *         \right|^2
 *       + \frac{ \omega^2 m^2 }{ 2 \varepsilon^4 }
 *         \left| \int dt\, \frac{ \mathbf{n \dot v} }
 *                               { (1 - \mathbf{n v})^2 }
 *                \exp [ i \omega' ( t - \mathbf{n r} ) ]
 *         \right|^2
 *     \right\}
 * @f]
 * with @f$ \hbar = 1 @f$ and @f$ c = 1 @f$, @f$ \mathbf{n} = \mathbf{k} / \omega @f$, and @f$
 * \varepsilon' = \varepsilon - \omega @f$, @f$ \omega' = \omega \varepsilon / \varepsilon' @f$.
 * One can note that
 * @f[
 *     \frac{d}{dt} \frac{1}{1 - \mathbf{n v}}
 *       = \frac{ \mathbf{n \dot v} }{ (1 - \mathbf{n v})^2 }, \quad
 *     \frac{d}{dt} \frac{ \mathbf{ n \times [ n \times v ] } }{ 1 - \mathbf{n v} }
 *       = \frac{ \mathbf{ n \times [(n - v) \times \dot v] } }{ (1 - \mathbf{n v})^2 },
 *       \quad
 *     \frac{d}{dt} \exp [ i \omega' ( t - \mathbf{n r} ) ]
 *       = i \omega' (1 - \mathbf{n v}) \exp [ i \omega' ( t - \mathbf{n r} ) ].
 * @f]
 * Hence, integrating by parts one get
 * @f[
 *     \frac{d^2 I}{d\omega d\Omega} = \frac{e^2 \omega'^2}{4 \pi^2} \left\{
 *         \frac{ \varepsilon^2 + \varepsilon'^2 }{ 2 \varepsilon^2 }
 *         \sum_\mathbf{p}
 *         \left| \int dt\, \mathbf{v p}
 *                \exp [ i \omega' ( t - \mathbf{n r} ) ]
 *         \right|^2
 *       + \frac{ \omega^2 m^2 }{ 2 \varepsilon^4 }
 *         \left| \int dt\,
 *                \exp [ i \omega' ( t - \mathbf{n r} ) ]
 *         \right|^2
 *     \right\},
 * @f]
 * where the product @f$ \mathbf{n \times [n \times v]} @f$ is rewritten using @f$ \mathbf{v p}_1
 * @f$ and @f$ \mathbf{v p}_2 @f$. Then the relation between the emission probability and the
 * distribution of the emitted energy (see #jackson1483_num) yields
 * @f[
 *     W_m = \frac{e^2 \pi \omega'^2}{4 \omega^3 V} \left\{
 *         \frac{ (\varepsilon^2 + \varepsilon'^2) }{ 2 \varepsilon^2 }
 *         \sum_\mathbf{p}
 *         \left| \int dt\, \mathbf{v p}
 *                \exp [ i \omega' ( t - \mathbf{n r} ) ]
 *         \right|^2
 *       + \frac{ \omega^2 m^2 }{ 2 \varepsilon^4 }
 *         \left| \int dt\,
 *                \exp [ i \omega' ( t - \mathbf{n r} ) ]
 *         \right|^2
 *     \right\}.
 * @f]
 * This formula is used here in the computations, together with #trap_rule.
 * @param vp1     the product of the mode polarization direction @f$ \mathbf{p}_1 @f$ and the
 *                electron velocity @f$ \mathbf{v} @f$, as the function of @f$ t @f$
 * @param vp2     the product @f$ \mathbf{v p}_2 @f$, as function of @f$ t @f$
 * @param nr      @f$ \mathbf{ n r }(t) @f$, with @f$ \mathbf{n} = \mathbf{k}_m / \omega_m @f$ and
 *                @f$ \mathbf{r}(t) @f$ the electron position
 * @param t_nodes values of @f$ t_0, t_1, t_2, ... @f$ used for the computation of the integrals;
 *                this function can take @p t_nodes as @p range or @p view (from C++20 ranges) or
 *                as @p std::vector, with doubles within
 * @param m       particle mass normalized to the electron mass
 * @param b       magnetic field used in the normalization, normalized to the critical field
 * @param gamma_p the Lorentz factor of the emitted particle
 * @param omega   @f$ 0 < \omega < \gamma_e / b @f$, frequency of the emitted photon, @f$ \omega_m
 *                @f$, normalized to the reverse radiation formation time
 */
template <typename V>
double bks_emission_probability( std::function<double(double)> vp1
                               , std::function<double(double)> vp2
                               , std::function<double(double)> nr
                               , V                             t_nodes
                               , double                        m
                               , double                        b
                               , double                        gamma_p
                               , double                        omega
                               ) {
    double epsilon = m * gamma_p / b; // the initial energy of the emitting particle in the
                                      // frequency units
    double epsilon_s = epsilon - omega;
    double omega_s   = omega * epsilon / epsilon_s;

    auto g  = std::function<std::complex<double>(double)>(
        [=](double t) {
            return exp( 1i * omega_s * ( t - nr(t) ) );
        }
    );
    auto f1 = std::function<std::complex<double>(double)>(
        [=](double t) {
            return vp1(t) * g(t);
        }
    );
    auto f2 = std::function<std::complex<double>(double)>(
        [=](double t) {
            return vp2(t) * g(t);
        }
    );

    std::complex<double> c1 = trap_rule(f1, t_nodes);
    std::complex<double> c2 = trap_rule(f2, t_nodes);
    std::complex<double> c3 = trap_rule( g, t_nodes);

    return alpha * M_PI / (8 * omega) * pow(omega_s / omega, 2)
               * ( (1 + pow(epsilon_s / epsilon, 2))   * (norm(c1) + norm(c2))
                 + pow(omega / (epsilon * gamma_p), 2) * norm(c3)
                 );
}

/**
 * An approximation of the vacuum refractive index in a strong magnetic field @f$
 * B @f$ for photons which polarization is @a perpendicular to the magnetic field (photon wave
 * vector is also assumed perpendicular to the magnetic field). The refractive index is
 * @f[
 *     n_\perp = 1 + \frac{\alpha b^2}{4 \pi} N_\perp(\chi),
 * @f]
 * with @f$ \chi @f$ the photon quantum parameter ("kappa").
 *
 * The exact expression for the refractive index can be get from [Narozhny N. B., Zh. Eksp. Teor.
 * Fiz. 55, 714 (1968)] and from (2.11) of [Ritus V. I., Sov. Phys. JETP 30, 1181 (1970)], see also
 * [Thomas Erber, High- Energy Electromagnetic Conversion Processes in Intense Magnetic Fields,
 * Reviews of modern physics, vol.  38, num. 4, oct. 1966] and Fig. 9 in [McDonald K. T. et al.,
 * Princeton U. preprint DOE ER, 3072-38 (1986)]. The asymptotics of the function @f$ \mathrm{Re}~N_\perp @f$
 * are
 * @f[
 *     \mathrm{Re}~N_\perp(\chi) = 
          \left\{
            \begin{array}{lr}
              14/45, \quad &\chi \ll 1, \\
              -1.65 \chi^{-4/3}, \quad &\chi \gg 1.
            \end{array}
          \right.
 * @f]
 * The approximation used here is based on numerical values of @f$ N @f$ obtained by Egor Sozinov
 * in the framework of QED (not published yet). The comparison between numerical data (red and blue
 * for perpendicular and parallel polarizations, respectively) and approximations used here (dark
 * green and brown) is shown below:
 * @image html "N_kappa_interpolation.png"
 * For the complex part...
 * @image html "N_kappa_interpolation.png"
 * @param b        the normalized magnetic field strength, @f$ B / B_{cr} @f$
 * @param omega    the normalized photon cyclic frequency, @f$ \omega t_{rf} @f$
*/
std::complex<double> vacuum_refractive_index_perp( double b
                                                 , double omega
                                                 ) {
    double k = omega * b * b; // kappa, i.e. chi of the photon
    double c1 = 2.6;
    double c2 = 60;
    double c3 = 0.07;
    double c4 = 0.9;
    double c5 = 1.0;
    double c6 = 0.008;
    double c7 = 0.045;
    double c8 = 5.2;
    double c9 = 15;
    double Re_N = 14.0/45 / (1 + pow(k / c1, 2))
                - 1.65 * pow(k + c9, -4.0/3) * k / (k + c2)
                + c3 * exp( -c4 * (k + c5 * c5 / k - 2 * c5))
                - c6 * exp( -c7 * (k + c8 * c8 / k - 2 * c8));
    double d1 = 1.12;
    double d2 = 1.3;
    double d3 = 0.26;
    double d4 = 128;
    double Im_N = 4 * M_PI * d1 * pow(3.0/8, 3.0/2) * exp( - 8.0 / ( 3 * k ) )
                / ( k * ( 1
                        + (d1 * pow(3.0/8, 3.0/2) * exp(0.0))
                          / 0.38 * (d2 + d3 * k / (k + d4)) * pow(k, 1.0/3)
                        )
                  );
    return 1.0 + alpha * b * b * (Re_N + 1i * Im_N) / (4 * M_PI);
}

/**
 * An approximation of the vacuum refractive index in a strong magnetic field @f$
 * B @f$ for photons which polarization is @a parallel to the magnetic field (photon wave vector is
 * assumed perpendicular to the magnetic field). The function @f$ N @f$ in this case has the
 * following asymptotics:
 * @f[
 *     \mathrm{Re} N_\parallel(\chi) = 
          \left\{
            \begin{array}{lr}
              8/45, \quad &\chi \ll 1, \\
              -1.1 \chi^{-4/3}, \quad &\chi \gg 1.
            \end{array}
          \right.
 * @f]
 * See #vacuum_refractive_index_perp for further details.
 * @param b        the normalized magnetic field strength, @f$ B / B_{cr} @f$
 * @param omega    the normalized photon cyclic frequency, @f$ \omega t_{rf} @f$
*/
std::complex<double> vacuum_refractive_index_parallel( double b
                                                     , double omega
                                                     ) {
    double k = omega * b * b; // kappa, i.e. chi of the photon
    double c1 = 3.0;
    double c2 = 25;
    double c3 = 0.046;
    double c4 = 1.9;
    double c5 = 1.0;
    double c6 = 0.0062;
    double c7 = 0.063;
    double c8 = 3.5;
    double c9 = 35;
    double Re_N = 14.0/45 / (1 + pow(k / c1, 2))
                - 1.65 * pow(k + c9, -4.0/3) * k / (k + c2)
                + c3 * exp( -c4 * (k + c5 * c5 / k - 2 * c5))
                - c6 * exp( -c7 * (k + c8 * c8 / k - 2 * c8));
    double d1 = 0.43;
    double d2 = 1.4;
    double d3 = 0.79;
    double d4 = 79;
    double Im_N = 4 * M_PI * d1 * pow(3.0/8, 3.0/2) * exp( - 8.0 / ( 3 * k ) )
                / ( k * ( 1
                        + (d1 * pow(3.0/8, 3.0/2) * exp(0.0))
                          / 0.38 * (d2 + d3 * k / (k + d4)) * pow(k, 1.0/3)
                        )
                  );
    return 1.0 + alpha * b * b * (Re_N + 1i * Im_N) / (4 * M_PI);
}

 /**
  * Probability of synchrotron emission in a medium with a given refractive index. The emission
  * phobability is computed with #bks_emission_probability, where the refractive index is added
  * only to the phase in the exponent.
  * @param ri      the medium refractive index, 1 for vacuum
  * @param m       the rest mass of the emitting particle (in electron masses), thus set to 1 for
  *                an electron
  * @param b       the magnetic field strength normalized to the Sauter-Schwinger field
  * @param gamma_p the Lorentz factor of the emitting particle
  * @param theta   angle in the plane perpendicular to the normal vector of the trajectory; @p
  *                theta = 0 is the direction of the tangent to the trajectory
  * @param omega   @f$ 0 < \omega < \gamma_p / b @f$, frequency of the emitted photon normalized to
  *                the reverse radiation formation time
 */
double bks_synchrotron_emission_probability( double ri
                                           , double m
                                           , double b
                                           , double gamma_p
                                           , double theta
                                           , double omega
                                           ) {
    double epsilon   = m * gamma_p / b; // m c^2 \gamma_p / (\hbar / t_rf)
    double epsilon_s = epsilon - omega;
    double omega_s   = omega * epsilon / epsilon_s;

    // the exponential in the integrals is approximately equal to exp(i phi(t)), with phi = 2 *
    // M_PI * (t * signed_reverse_tau_parallel + pow(t / tau_perp, 3));
    double delta_ri = ri - 1;
    // +-1 / tau_parallel, to avoid the point tau_parallel = \infty
    double signed_reverse_tau_parallel
        = omega_s / (4 * M_PI) * ( theta * theta + 1 / (gamma_p * gamma_p) - 2 * delta_ri);
    double reverse_tau_parallel = fabs(signed_reverse_tau_parallel);
    double varsigma = (signed_reverse_tau_parallel >= 0) ? (+1) : (-1);
    double tau_perp = gamma_p * pow(12 * M_PI * m * m / (omega_s * r(gamma_p)), 1/3.0);

    // we integrate approximately on the interval of +-l scales of the descent of the integrals
    double l = 3;
    // fraction of the minimal period which determines the timestep; irrational to avoid resonances
    double period_fraction = 1 / (exp(1) - 0.5);

    // a point where linear and cubic terms in the phase yield the same oscillation period; if
    // varsigma = -1, then d\phi / dt = 0 at t = +-t_s
    double ts = sqrt(pow(tau_perp, 3) * reverse_tau_parallel / 3);
    // tw is a width of a leading bump (one of two) in the case varsigma = -1 and tau_perp >>
    // tau_parallel, i.e. tw = T(t_s) in this case
    double tw = sqrt(pow(tau_perp, 3) / ts);

    if ( varsigma < 0 and l * tw < ts ) {
        // parabolic approximation of phi(t) at t_s works, and the integration should be performed
        // over two intervals (one around +t_s and the other around -t_s);
        // tb is the upper limit of the integration
        double tb = ts + l * tw;
        // ta is the lower limit of the integration
        double ta = ts - l * tw;
        // the estimate of the period of the exponent oscillations
        // (T = 2 \pi / (d\phi/dt)) at t = tb
        double osc_period = 1 / (-reverse_tau_parallel + 3 * pow(tb / tau_perp, 2) / tau_perp);
        // number of points for the exponent integration
        long long int nt = llround((tb - ta) / (osc_period * period_fraction));
        // step of the integration
        double dt = (tb - ta) / static_cast<double>(nt - 1);
        auto t_nodes  = ranges::v3::iota_view(0, nt)
                      | ranges::v3::views::transform(
                          [=](long long i){ return  ta + dt * static_cast<double>(i); }
                        );
        auto t_nodes_ = ranges::v3::iota_view(0, nt)
                      | ranges::v3::views::transform(
                          [=](long long i){ return -tb + dt * static_cast<double>(i); }
                        );

        auto nr = std::function<double(double)>(
            [=](double t) {
                return ri * m * r(gamma_p) * sin(t / (m * gamma_p)) * cos(theta);
            }
        );

        auto vp1 = std::function<double(double)>(
            [=](double t) {
                return r(gamma_p) / gamma_p * sin(t / (m * gamma_p))
                     * 0.25 * (1 - tanh(8 * ((t - ts) / (l * tw) - 0.7)))
                     *        (1 + tanh(8 * ((t - ts) / (l * tw) + 0.7)));
                /* qwe
                     * 0.25 * (1 - tanh(8 * (t / tb - 0.7)))
                     *        (1 + tanh(8 * (t / tb + 0.7)));
                     */
                     // attenuation; see comments in #synchrotron_emission_probability
            }
        );
        auto vp1_ = std::function<double(double)>(
            [=](double t) {
                return r(gamma_p) / gamma_p * sin(t / (m * gamma_p))
                     * 0.25 * (1 - tanh(8 * ((t + ts) / (l * tw) - 0.7)))
                     *        (1 + tanh(8 * ((t + ts) / (l * tw) + 0.7)));
                /* qwe
                     * 0.25 * (1 - tanh(8 * (t / tb - 0.7)))
                     *        (1 + tanh(8 * (t / tb + 0.7)));
                     */
                     // attenuation; see comments in #synchrotron_emission_probability
            }
        );

        auto vp2 = std::function<double(double)>(
            [=](double t) {
                return r(gamma_p) / gamma_p * sin(theta) * cos(t / (m * gamma_p))
                     * 0.25 * (1 - tanh(8 * ((t - ts) / (l * tw) - 0.7)))
                     *        (1 + tanh(8 * ((t - ts) / (l * tw) + 0.7)));
                /* qwe
                     * 0.25 * (1 - tanh(8 * (t / tb - 0.7)))
                     *        (1 + tanh(8 * (t / tb + 0.7)));
                     */
                     // attenuation; see comments in #synchrotron_emission_probability
            }
        );
        auto vp2_ = std::function<double(double)>(
            [=](double t) {
                return r(gamma_p) / gamma_p * sin(theta) * cos(t / (m * gamma_p))
                     * 0.25 * (1 - tanh(8 * ((t + ts) / (l * tw) - 0.7)))
                     *        (1 + tanh(8 * ((t + ts) / (l * tw) + 0.7)));
                /* qwe
                     * 0.25 * (1 - tanh(8 * (t / tb - 0.7)))
                     *        (1 + tanh(8 * (t / tb + 0.7)));
                     */
                     // attenuation; see comments in #synchrotron_emission_probability
            }
        );
        return bks_emission_probability(vp1 , vp2 , nr, t_nodes , m, b, gamma_p, omega)
             + bks_emission_probability(vp1_, vp2_, nr, t_nodes_, m, b, gamma_p, omega);
    } else {
        // parabolic approximation of phi(t) at t_s doesn't work, and we should integrate over
        // single interval;
        // tb is the upper limit of the integration
        double tb = l * std::max(tau_perp, ts);
        // the estimate of the period of the exponent oscillations
        // (T = 2 \pi / (d\phi/dt)) at t = tb
        double osc_period = 1 / (reverse_tau_parallel + 3 * pow(tb / tau_perp, 2) / tau_perp);
        // lower limit of the integration
        double ta = -tb;
        // number of points for the exponent integration
        long long int nt = llround((tb - ta) / (osc_period * period_fraction));
        // step of the integration
        double dt = (tb - ta) / static_cast<double>(nt - 1);
        auto t_nodes = ranges::v3::iota_view(0, nt)
                     | ranges::v3::views::transform(
                           [=](long long i){ return ta + dt * static_cast<double>(i); }
                       );

        auto nr = std::function<double(double)>(
            [=](double t) {
                return ri * m * r(gamma_p) * sin(t / (m * gamma_p)) * cos(theta);
            }
        );
        auto vp1 = std::function<double(double)>(
            [=](double t) {
                return r(gamma_p) / gamma_p * sin(t / (m * gamma_p))
                     * 0.25 * (1 - tanh(8 * (t / tb - 0.7)))
                     *        (1 + tanh(8 * (t / tb + 0.7)));
                     // attenuation; see comments in #synchrotron_emission_probability
            }
        );
        auto vp2 = std::function<double(double)>(
            [=](double t) {
                return r(gamma_p) / gamma_p * sin(theta) * cos(t / (m * gamma_p))
                     * 0.25 * (1 - tanh(8 * (t / tb - 0.7)))
                     *        (1 + tanh(8 * (t / tb + 0.7)));
                     // attenuation; see comments in #synchrotron_emission_probability
            }
        );
        return bks_emission_probability(vp1, vp2, nr, t_nodes, m, b, gamma_p, omega);
    }
}

/**
 * Target distribution which can be used in the Metropolis algorithm in order to generate the
 * distributions of synchrotron photons in (@p theta, @p omega) space, i.e. this function gives
 * probability of photon emission per frequency interval and per angle interval. Here @p theta is
 * the angle in the plane perpendicular to the normal vector of the trajectory, and <tt> theta = 0
 * </tt> is the direction of the tangent to the trajectory; @p omega is the normalized photon
 * frequency.
 *
 * The normalization of this function (which does not matter for the Metropolis algorithm, but can
 * be important for other applications) is chosen such that @p bks_synchrotron_td integrated over
 * @p theta and @p omega yield the overall probability of photon emission during single full circle
 * path of the emitting particle.
 *
 * For parameters description, see #bks_synchrotron_emission_probability.
 */
std::function< double( std::tuple<double, double> ) >
bks_synchrotron_td( double ri
                  , double m
                  , double b
                  , double gamma_p
                  ) {
    // The normalization can be obtained analogous to jackson1483_num, see comments there on the
    // mode density etc. Note that the generated photons are evenly distributed around the electron
    // trajectory (which is a circle), i.e. they are evenly distributed in the angle $ \phi $ which
    // in turn corresponds to the direction of the photon wavevector in the plane of the electron
    // trajectory. Thus, the solid angle is $ d\Omega = \cos \theta \, d\phi d\theta $, and for the
    // emission probability we get $ d^2 W / d\theta d\omega = 2 W_m \cos \theta \omega^2 / \pi^2 $
    // with $ W_m $ the emission probability computed with bks_synchrotron_emission_probability
    // function.
    return std::function< double( std::tuple<double, double> ) >(
        [=]( std::tuple<double, double> theta_omega ) {
            double theta = std::get<0>(theta_omega);
            double omega = std::get<1>(theta_omega);
            if (omega > 0 and omega * b < m * gamma_p) {
                return 2 * pow(omega / M_PI, 2) * cos(theta) / m
                         * bks_synchrotron_emission_probability(ri, m, b, gamma_p, theta, omega);
            } else {
                return 0.0;
            }
        }
    );
}

/**
 * The (normalized) energy radiated per unit frequency interval per unit solid angle, computed
 * numerically from #bks_synchrotron_emission_probability. In contrast to #jackson1483_num,
 * #bks_jackson1483_num take into account the spin term and the radiation recoil effect.
 * @param b       the normalized magnetic field strength
 * @param gamma_e the electron Lorentz factor
 * @param theta   angle in the plane perpendicular to the normal vector of the trajectory; @p theta
 *                = 0 is the direction of the tangent to the trajectory
 * @param omega   > 0, normalized frequency of the emitted photon
 */
double bks_jackson1483_num(double b, double gamma_e, double theta, double omega) {
    std::function<double(std::tuple<double, double>)> probability_distribution
        = bks_synchrotron_td(1, 1, b, gamma_e);
    return b * omega / (2 * M_PI)
             * probability_distribution(std::make_tuple(theta, omega));
}

/**
 * @}
 */

/**
 * @defgroup synchrotron_complex_n Synchrotron emission for complex refractive index
 * * @brief Classical radiation reaction formulas for the synchrotron emission with imaginary part
 * of the refractive index taken into account. The same normalization as in @ref radiation is used
 * here.
 * @{
 */

 /**
  * (Positive) kernel of the @a synchrotron power of the radiation losses in a medium with a given @a complex
  * refractive index, in the framework of classical electrodynamics. If the imaginary part of the
  * refractive index is zero, this kernel is the spectal distribution of the emitted power,
  * multiplied by (-1). If the imaginary part of the refractive index is non-zero, the emitted
  * photons decay (also during the radiation process) and there is no obvious physical sense of the
  * spectrum. Anyway, the power of the losses can be given as an integral over real "wave vector" k
  * (from 0 to infinity in the classical case):
  * @f[
  *     \frac{dP}{dk} = -\frac{\alpha}{\pi} b k \int_0^\infty
  *     \frac{ {\mathbf v}(\tau/2) {\mathbf v}(-\tau/2) }
  *          { |{\mathbf r}(\tau/2) - {\mathbf r}(-\tau/2)| }
  *     \exp(-\psi \nu \tau) \\
  *     \times \left \{ \eta \sin [\psi \eta \tau - k | {\mathbf r}(\tau/2) - {\mathbf r}(-\tau/2)|]
  *                   + \nu  \cos [\psi \eta \tau - k | {\mathbf r}(\tau/2) - {\mathbf r}(-\tau/2)|]
  *            \right\} \, d\tau,
  * @f]
  * with @f$ \eta @f$ and @f$ \nu @f$ the real and imaginary parts of the refractive index @f$ n =
  * \eta + i \nu @f$, @f$ \psi = k / |n|^2 @f$, and we assume the pre-exponent multiplier @f$ \mu / \eta
  * \approx 1 @f$ (where @f$ \mu @f$ the magnetic permeability).
  *
  * Note that the formula above is not general and valid only for the circular motion of a
  * particle. For an ultrarelativistic particle the formula for the kernel can be further
  * simplified:
  * @f[
  * @f]
  * @param ri      complex refractive index
  * @param m       the rest mass of the emitting particle (in electron masses)
  * @param b       the magnetic field strength normalized to the Sauter-Schwinger field
  * @param gamma_p the Lorentz factor of the emitting particle
  * @param k       @f$ k > 0 @f$, the wavevector amplitude used in the kernel, normalized to the
  *                reverse radiation formation length @f$ 1 / c t_{rf} @f$
 */
double classical_synchrotron_power_kernel_deprecated( std::complex<double> ri
                                                    , double m
                                                    , double b
                                                    , double gamma_p
                                                    , double k
                                                    ) {
    /*double epsilon   = m * gamma_p / b; // m c^2 \gamma_p / (\hbar / t_rf)
    double epsilon_s = epsilon - omega;
    double omega_s   = omega * epsilon / epsilon_s;
    */
    double eta = real(ri);
    double nu  = imag(ri);

    // the phase in the sine and cosine functions is phi = 2 * M_PI
    // * (tau * signed_reverse_tau_parallel + pow(tau / tau_perp, 3));
    // +-1 / tau_parallel, to avoid the point tau_parallel = \infty
    double signed_reverse_tau_parallel
        = k / (2 * M_PI) * (eta / norm(ri) - 1 + 1 / (2 * gamma_p * gamma_p));
    double reverse_tau_parallel = fabs(signed_reverse_tau_parallel);
    double varsigma = (signed_reverse_tau_parallel >= 0) ? (+1) : (-1);
    double tau_perp = gamma_p * pow(48 * M_PI * m * m / (k * r(gamma_p)), 1/3.0);
    double reverse_tau_d = 1 / (k * nu);

    // we integrate approximately on the interval from 0 to l scales of the descent of the
    // integrals, and on +-l scales around the saddle point if the Cherenkov condition is
    // exceeded
    double l = 3;
    double period_fraction = 1/2.0; // fraction of the minimal period which determines the timestep

    // a point where linear and cubic terms in the phase yield the same oscillation period; if
    // varsigma = -1, then d\phi / dt = 0 at t = +-t_s (saddle points)
    double ts = sqrt(pow(tau_perp, 3) * reverse_tau_parallel / 3);
    // tw is a width of a leading bump (one of two) in the case varsigma = -1 and tau_perp >>
    // tau_parallel, i.e. tw = T(t_s) in this case
    double tw = sqrt(pow(tau_perp, 3) / ts);
    // common multiplier
    double a = alpha * b * k / M_PI;

    if ( varsigma < 0 and l * tw < ts ) {
        // parabolic approximation of phi(t) at t_s works, and the integration should be performed
        // around t_s
        // tb is the upper limit of the integration
        double tb = ts + l * tw;
        // ta is the lower limit of the integration, ta > 0
        double ta = ts - l * tw;
        // the estimate of the period of the exponent oscillations
        // (T = 2 \pi / (d\phi/dt)) at t = tb
        double osc_period = 1 / (-reverse_tau_parallel + 3 * pow(tb / tau_perp, 2) / tau_perp);
        // number of points for the exponent integration
        long long int nt = llround((tb - ta) / (osc_period * period_fraction));
        // step of the integration
        double dt = (tb - ta) / static_cast<double>(nt - 1);
        auto t_nodes  = ranges::v3::iota_view(0, nt)
                      | ranges::v3::views::transform(
                          [=](long long i){ return  ta + dt * static_cast<double>(i); }
                        );

        auto phi = std::function<double(double)>(
            [=](double t) {
                return 2 * M_PI * ( t * signed_reverse_tau_parallel + pow(t / tau_perp, 3) );
            }
        );

        auto f = std::function<double(double)>(
            [=](double t) {
                return exp(-t * reverse_tau_d) * sin(phi(t)) / t;
            }
        );

        return a * trap_rule(f, t_nodes);
    } else {
        // parabolic approximation of phi(t) at t_s doesn't work, and we should start to integrate
        // from 0
        // tb is the upper limit of the integration
        double tb = l * std::max(tau_perp, ts);
        // the estimate of the period of the exponent oscillations
        // (T = 2 \pi / (d\phi/dt)) at t = tb
        double osc_period = 1 / (reverse_tau_parallel + 3 * pow(tb / tau_perp, 2) / tau_perp);
        // lower limit of the integration
        double ta = osc_period * period_fraction * 1e-2; // to avoid sin(0) / 0
        // number of points for the exponent integration
        long long int nt = llround((tb - ta) / (osc_period * period_fraction));
        // step of the integration
        double dt = (tb - ta) / static_cast<double>(nt - 1);
        auto t_nodes = ranges::v3::iota_view(0, nt)
                     | ranges::v3::views::transform(
                           [=](long long i){ return ta + dt * static_cast<double>(i); }
                       );

        auto f = std::function<double(double)>(
            [=](double t) {
                double phi =  2 * M_PI
                           * ( t * signed_reverse_tau_parallel + pow(t / tau_perp, 3) );
                //double vv = 1 / pow(gamma_p, 2) + 0.5 * pow(t / (m * gamma_p), 2);
                //return vv * exp(-t * reverse_tau_d) * sin(phi) / t;
                double vv = t * t;
                return vv * exp(-t * reverse_tau_d) * sin(phi) / t;
            }
        );

        return a * trap_rule(f, t_nodes);
    }
}

 /**
  * WARNING: похоже, здесь была потеряна двойка, множетель 2 введён в формулу по сравнению с
  * предыдущей версией;
  * Kernel..., should be integrated over theta and k. Note that here the angle theta is not the
  * same in bks...
  * The magic consts are such that classical_syn... test
  * WARNING: we do not compute for high nt, nt > some -> int = 0
  * @param ri      the medium refractive index (complex)
  * @param m       the rest mass of the emitting particle (in electron masses)
  * @param b       the magnetic field strength normalized to the Sauter-Schwinger field
  * @param gamma_p the Lorentz factor of the emitting particle
  * @param theta   @f$ theta > 0 @f$, the angle between the tangent to the trajectory and the
  *                "wavevector" k. theta = 0 is the direction of the tangent to the trajectory. The
  *                kernel here is already integrated over the azimuthal angle
  * @param k       @f$ k > 0 @f$, the wavevector amplitude used in the kernel, normalized to the
  *                reverse radiation formation length @f$ 1 / c t_{rf} @f$
 */
double classical_synchrotron_power_kernel( std::complex<double> ri
                                         , double m
                                         , double b
                                         , double gamma_p
                                         , double theta
                                         , double k
                                         ) {
    double eta = real(ri);
    double nu  = imag(ri);

    double signed_reverse_tau_parallel
        = k / (2 * M_PI) * ( theta * theta / 2 + eta / norm(ri) - 1 + 1 / (2 * gamma_p * gamma_p));
    double reverse_tau_parallel = fabs(signed_reverse_tau_parallel);
    double varsigma = (signed_reverse_tau_parallel >= 0) ? (+1) : (-1);
    double tau_perp = gamma_p * pow(48 * M_PI * m * m / (k * r(gamma_p)), 1/3.0);
    double reverse_tau_d = k * nu;

    // we integrate approximately on +-l scales around the saddle point if the Cherenkov condition
    // is exceeded
    double l = 2;
    double period_fraction = 1/2.0; // fraction of the minimal period which determines the
                                    // timestep, in the Cherenkov case
    long long int nt_max = 80'000;

    // a point where linear and cubic terms in the phase yield the same oscillation period; if
    // varsigma = -1, then d\phi / dt = 0 at t = t_s
    double ts = sqrt(pow(tau_perp, 3) * reverse_tau_parallel / 3);
    // tw is a width of a leading bump in the case varsigma = -1 and tau_perp >>
    // tau_parallel, i.e. tw = T(t_s) in this case
    double tw = sqrt(pow(tau_perp, 3) / ts);

    // if true, the contribution from the bump ts +- l * tw is greater than the contribution from
    // [0, td].
    bool bump = reverse_tau_parallel * tw * exp(-ts * reverse_tau_d) > 1;

    if ( varsigma < 0 and l * tw < ts and bump) {
        // parabolic approximation of phi(t) at t_s works, and the integration should be performed
        // over the interval around t_s;
        // tb is the upper limit of the integration
        double tb = ts + l * tw;
        // ta is the lower limit of the integration
        double ta = ts - l * tw;
        // the estimate of the period of the exponent oscillations
        // (T = 2 \pi / (d\phi/dt)) at t = tb
        double osc_period = 1 / (-reverse_tau_parallel + 3 * pow(tb / tau_perp, 2) / tau_perp);
        // number of points for the exponent integration
        long long int nt = llround((tb - ta) / (osc_period * period_fraction));
        // step of the integration
        double dt = (tb - ta) / static_cast<double>(nt - 1);

        auto t_nodes  = ranges::v3::iota_view(0, nt)
                      | ranges::v3::views::transform(
                          [=](long long i){ return  ta + dt * static_cast<double>(i); }
                        );

        auto f = std::function<double(double)>(
            [=](double t) {
                double phi = 2 * M_PI * ( t * signed_reverse_tau_parallel + pow(t / tau_perp, 3) );
                return alpha * b / M_PI * k * k * theta
                     * (theta * theta - t * t / (4 * gamma_p * gamma_p * m * m))
                     * exp(-t * reverse_tau_d)
                     * ( eta * cos(phi) - nu * sin(phi) )
                     * 0.25 * (1 - tanh(8 * ((t - ts) / (l * tw) - 0.7)))
                     *        (1 + tanh(8 * ((t - ts) / (l * tw) + 0.7)));
//                     * 0.5 * (1 - tanh(8 * (t / tb - 0.7))); // qwe QWE QWE !!! wrong in this case
                     // note that artificial attenuation is added here; the idea behind is that
                     // with this attenuation neighboring bumps quench each other earlier hence
                     // smaller integration interval can be used
            }
        );

        if (nt < nt_max) {
            return trap_rule(f, t_nodes);
        } else {
            return 0;
        }

    } else if (reverse_tau_d * std::max(tau_perp, ts) < 1 ) {
        // parabolic approximation of phi(t) at t_s doesn't work; photon decay is not crucial

        // we integrate in this (more "synchrotron") case of the interval from 0 to l_s scales of
        // the descent of the integrals
        double l_s = 11;
        double period_fraction_s = 1 / exp(1); // fraction of the minimal period which
                                               // determines the timestep; irrational to
                                               // avoid resonances

        // tb is the upper limit of the integration
        /*double tb = l_s * ( reverse_tau_d == 0 ? // qwe we assume that nu > 0 !!!
                            std::max(tau_perp, ts) :
                            std::min( 1 / reverse_tau_d, std::max(tau_perp, ts) ) );
                            */
        double tb = l_s * std::max(tau_perp, ts);

        // the estimate of the period of the exponent oscillations
        // (T = 2 \pi / (d\phi/dt)) at t = tb
        double osc_period = 1 / (reverse_tau_parallel + 3 * pow(tb / tau_perp, 2) / tau_perp);
        // lower limit of the integration
        double ta = 0;
        // number of points for the exponent integration
        long long int nt = llround((tb - ta) / (osc_period * period_fraction_s));
        // step of the integration
        double dt = (tb - ta) / static_cast<double>(nt - 1);
        auto t_nodes = ranges::v3::iota_view(0, nt)
                     | ranges::v3::views::transform(
                           [=](long long i){ return ta + dt * static_cast<double>(i); }
                       );

    // qwe
    /*
    if(reverse_tau_d != 0)
    std::cout << k * b * b << '\t' << theta/1e-3 << '\t' << 1/reverse_tau_d << '\t' << sqrt(pow(tau_perp, 3) * reverse_tau_parallel) << '\t' << tb << '\t' << osc_period << std::endl;
    */

        auto f = std::function<double(double)>(
            [=](double t) {
                double phi = 2 * M_PI * ( t * signed_reverse_tau_parallel + pow(t / tau_perp, 3) );
                return alpha * b / M_PI * k * k * theta
                     * (theta * theta - t * t / (4 * gamma_p * gamma_p * m * m))
                     * exp(-t * reverse_tau_d)
                     * ( eta * cos(phi) - nu * sin(phi) )
                     * 0.5 * (1 - tanh(9 * (t / tb - 0.55)));
                     // note that artificial attenuation is added here; the idea behind is that
                     // with this attenuation neighboring bumps quench each other earlier hence
                     // smaller integration interval can be used
            }
        );
        auto g = std::function<double(double)>( // qwe
            [=](double t) {
                double phi = 2 * M_PI * ( t * signed_reverse_tau_parallel + pow(t / tau_perp, 3) );
                return alpha * b / M_PI * k * k * theta
                     * (theta * theta - t * t / (4 * gamma_p * gamma_p * m * m))
                     * ( eta * cos(phi) - nu * sin(phi) )
                     * 0.5 * (1 - tanh(9 * (t / tb - 0.55)));
            }
        );

        if (nt < nt_max) {
            return trap_rule(f, t_nodes);
        } else {
            return 0;
        }

    } else {
        //  photon decay is crucial; parabolic approximation of phi(t) at t_s doesn't work;

        // we integrate in this (more "synchrotron") case of the interval from 0 to l_s scales of
        // the descent of the integrals
        double l_s = 11;
        double period_fraction_s = 1 / (30 * exp(1)); // fraction of the minimal period which
                                               // determines the timestep; irrational to
                                               // avoid resonances

        // tb is the upper limit of the integration
        double tb = l_s / reverse_tau_d;

        // the estimate of the period of the exponent oscillations
        // (T = 2 \pi / (d\phi/dt)) at t = tb
        double osc_period = std::min( 1 / (reverse_tau_parallel + 3 * pow(tb / tau_perp, 2) / tau_perp)
                                    , 1 / reverse_tau_d );
        // lower limit of the integration
        double ta = 0;
        // number of points for the exponent integration
        long long int nt = llround((tb - ta) / (osc_period * period_fraction_s));
        // step of the integration
        double dt = (tb - ta) / static_cast<double>(nt - 1);
        auto t_nodes = ranges::v3::iota_view(0, nt)
                     | ranges::v3::views::transform(
                           [=](long long i){ return ta + dt * static_cast<double>(i); }
                       );

        auto f = std::function<double(double)>(
            [=](double t) {
                double phi = 2 * M_PI * ( t * signed_reverse_tau_parallel + pow(t / tau_perp, 3) );
                return alpha * b / M_PI * k * k * theta
                     * (theta * theta - t * t / (4 * gamma_p * gamma_p * m * m))
                     * exp(-t * reverse_tau_d)
                     * ( eta * cos(phi) - nu * sin(phi) );
            }
        );

        if (nt < nt_max) {
            return trap_rule(f, t_nodes);
        } else {
            return 0;
        }

    }
}

/**
 * @}
 */

