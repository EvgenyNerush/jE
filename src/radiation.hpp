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
 * coordinates are normalized to @f$ c t_{rf} @f$ and energy to @f$ m c^2 @f$. Velocity is
 * normalized to the speed of light @f$ c @f$ and frequency is normalized to @f$ 1 / t_{rf} @f$.
 * Note that the photon of the cyclic frequency @f$ \omega @f$ in the normalized units has energy
 * @f$ \omega \, b @f$ with @f$ b = B / B_S @f$, and @f$ B_S = m^2 c^3 / e \hbar @f$ the
 * Sauter--Schwinger field (the critical field of QED). For the synchrotrom radiation formulas, see
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
 * @param pv      the product of the mode polarization direction @f$ \mathbf{p} @f$ and the
 *                electron velocity @f$ \mathbf{v} @f$, as the function of @f$ t @f$.
 * @param nr      @f$ \mathbf{ n \cdot r }(t) @f$, with @f$ \mathbf{n} = \mathbf{k}_m / k_m @f$ and
 *                @f$ \mathbf{r}(t) @f$ the electron position
 * @param t_nodes values of @f$ t_0, t_1, t_2, ... @f$ used for the integral computation; this
 *                function can take @p t_nodes as @p range or @p view (from C++20 ranges) or as @p
 *                std::vector, with doubles within
 * @param omega   > 0, normalized frequency of the emitted photon, @f$ \omega_m @f$
 */
template <typename V>
std::complex<double> c_m( std::function<double(double)> pv
                   , std::function<double(double)> nr
                   , V t_nodes
                   , double omega
                   ) {
    auto f = std::function<std::complex<double>(double)>(
        [=](double t) {
            return pv(t) *  exp( 1i * omega * ( t - nr(t) ) );
        }
    );
    return sqrt(alpha * M_PI / (4 * omega)) * trap_rule(f, t_nodes);
}

/*
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
    // the exponent in the integral is approximately equal to exp(i phi(t))
    double tau_parallel = 4 * M_PI / (omega * (theta * theta + 1 / (gamma_e * gamma_e)));
    double tau_perp     = gamma_e * pow(12 * M_PI / (omega * r(gamma_e)), 1/3.0);
    auto phi = [=](double t){ return 2 * M_PI * (t / tau_parallel + pow(t / tau_perp, 3)); };
    // we integrate approximately on +-l periods of the exponent
    int l = 5;
    auto f = std::function<double(const double&)>(
                 [=](const double& t) {
                      return phi(t) - 2 * M_PI * static_cast<double>(l);
                 }
             );

    double tb;
    if (tau_perp < tau_parallel) {
        tb = tau_perp * pow(static_cast<double>(l), 1/3.0);
    } else {
        tb = tau_parallel * static_cast<double>(l);
    } //... thus the root of f(t) = 0 is in [0, tb]
    int n_bisections = 10; // number of the iterations in the bisection method; 1 / 2^10 ~ 10^{-3}
    tb = bisection(f, 0, tb, n_bisections).value();
    double ta = -tb;

    // the estimate of the period of the exponent oscillations (T = 2 \pi / (d\phi/dt)) at t = tb
    double osc_period = 1 / (1 / tau_parallel + 3 * tb * pow(tb / tau_perp, 2));
    // number of point for the exponent integration
    long long int nt = llround(3 * (tb - ta) / osc_period);
    // step of the integration
    double dt = (tb - ta) / static_cast<double>(nt - 1);
    auto t_nodes  = ranges::v3::iota_view(0, nt)
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
        auto pv = std::function<double(double)>(
            [=](double t) {
                return r(gamma_e) / gamma_e * sin(t / gamma_e)
                    * exp( -8 * pow(t / tb, 8) ); // note that artificial attenuation is added here
                                                  // to avoid emission caused by the current on and
                                                  // off
            }
        );
        c = c_m(pv, nr, t_nodes, omega);
    } else {
        auto pv = std::function<double(double)>(
            [=](double t) {
                return r(gamma_e) / gamma_e * sin(theta) * cos(t / gamma_e)
                    * exp( -8 * pow(t / tb, 8) ); // note that artificial attenuation is added here
                                                  // to avoid emission caused by the current on and
                                                  // off
            }
        );
        c = c_m(pv, nr, t_nodes, omega);
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
 *     = \frac{e^2}{4 \pi^2 c^3} \int \int \omega^2 \sum_\mathbf{p} \left| \int_t \mathbf{v}(t)
 *     \mathbf{p} \exp \left[ i \omega t - i \mathbf{k r}(t) \right] \right|^2 \, d\omega \,
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
 * The approximation for the real part of the index of refraction the exact expression for which was taken from [Thomas Erber, High-
 * Energy Electromagnetic Conversion Processes in Intense Magnetic Fields, reviews of modern physics, vol 38, num 4, oct. 1966]
 * The expression depends on @f$ \chi_{ph} = \omega * H/H_{ph} @f$, where H - magnetic field and @f$ \omega @f$ is the photon
 * angular frequency.
 * @param h        the normalized magnetic field strength
 * @param omega    the photon angular frequency
*/
long double index_n(double h, double omega) {
    double chi = h * omega;
    return 1.0 + (chi <= 2.0) * 0.22 * exp(1.3636 * chi*chi - 1.26 * pow(chi, 2.22)) + (chi > 2.0) * 0.365756 / pow(chi, 4.0/3.0);
}

/**
 * @}
 */
