/**
 * @file
 * @brief Formulas for synchrotron emission
 * @details Formulas describing the emission of photons by ultrarelativistic electrons rotating in
 * magnetic field @f$ B @f$. In functions provided by this file time is normalized to the radiation
 * formation time 
 * @f[
 * t_{rf} = m c / (e B),
 * @f]
 * coordinates to @f$ c t_{rf} @f$ and energy to @f$ m c^2 @f$. Thus, velocity is normalized to the
 * speed of light @f$ c @f$ and frequency is normalized to @f$ 1 / t_{rf} @f$. The electron
 * trajectory is a circle which in the normalized units is
 * @f[
 * x = r \sin( t / \gamma_e ), \\
 * @f]
 * @f[
 * y = r \cos( t / \gamma_e ),
 * @f]
 * with @f$ r = \sqrt{\gamma_e^2 - 1} @f$ the curvature radius and @f$ \gamma_e @f$ the electron
 * Lorentz factor.  Note that the photon of the cyclic frequency @f$ \omega @f$ in the normalized
 * units has energy @f$ \omega \, b @f$ with @f$ b = B / B_S @f$, and @f$ B_S = m^2 c^3 / e \hbar
 * @f$ the Sauter--Schwinger field (the critical field of QED).
 */

#include <complex>
#include <cmath>

using namespace std::literals::complex_literals;

/**
 * The fine structure constant, CODATA 2018 recommended value
 */
const double alpha = 1 / 137.035999084;

/**
 * The normalized curvature radius of the electron trajectory
 * @param gamma_e the electron Lorentz factor
 */
double r(double gamma_e) {
    return sqrt(gamma_e * gamma_e - 1);
}

/**
 * The normalized critical frequency of the classical synchrotron emission, see Eq.  (14.85) from
 * [Jackson J.D., Classical Electrodynamics, Wiley, 1962].
 * @param gamma_e the electron Lorentz factor
 */
double omega_c(double gamma_e) {
    return 3 * pow(gamma_e, 3) / r(gamma_e);
}

/**
 * The (normalized) energy radiated per unit frequency interval per unit solid angle, as defined in
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
 * The product of the probability of the photon emission and the volume of the virtual box, @f$ V
 * |C_m|^2 @f$, computed numerically, as follows.
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
 * the emitting electron.
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
 * This formula is used here to compute @f$ C_m @f$ with midpoint rule of integration (very similar to the
 * trapezoidal rule). Note that the product @f$ V |C_m|^2 @f$ do not contain the volume of the
 * virtual box.
 * @param a       value which shrinks/enlarges the limits of integration, i.e. the ratio of the
 *                upper (lower) limit of integration to the "transverse" scale of the exponent to
 *                integrate, @f$ \tau_\perp = (m c \gamma_e / e B) (6 \pi c / r \omega)^{1/3} @f$.
 * @param n       the number of points for numerical integration; use #emission_probability_s if
 *                not know what value to choose
 * @param b       the normalized magnetic field strength
 * @param gamma_e the electron Lorentz factor
 * @param theta   angle in the plane perpendicular to the normal vector of the trajectory; @p theta
 *                = 0 is the direction of the tangent to the trajectory
 * @param p       the polarization of the emitted wave: <tt> p = true </tt> means polarization
 *                parallel to the normal vector of the trajectory, and <tt> p = false </tt> means
 *                polarization in the plane perpendicular to the normal vector
 * @param omega   > 0, normalized frequency of the emitted photon
 */
double emission_probability(double a, long long n,
        double b, double gamma_e, double theta, bool p, double omega) {
    // one of two scales of the exponent which we integrate
    double tau_perp = gamma_e * pow(6 * M_PI / (omega * r(gamma_e)), 1/3.0);
    // the limits of the numerical integration are [-tau_m, tau_m]
    double tau_m = a * tau_perp;
    // step of the integration
    double dt = 2 * tau_m / static_cast<double>(n);
    // c_m
    std::complex<double> c = 0i;
    if (p == true) {
        for (long long int j = 0; j < n; ++j) {
            double t = -tau_m + dt * (0.5 + static_cast<double>(j));
            c -= dt * r(gamma_e) / gamma_e
                * sin(t / gamma_e) * exp( 1i * omega * t
                                        - 1i * omega * r(gamma_e) * sin(t / gamma_e) * cos(theta));
        }
    } else {
        for (long long int j = 0; j < n; ++j) {
            double t = -tau_m + dt * (0.5 + static_cast<double>(j));
            c += dt * r(gamma_e) / gamma_e * sin(theta)
                * cos(t / gamma_e) * exp( 1i * omega * t
                                        - 1i * omega * r(gamma_e) * sin(t / gamma_e) * cos(theta));
        }
    }
    return M_PI * alpha / (4 * omega) * norm(c);
}

/**
 * Simplified version of #emission_probability, where @p a and @p n are set to ensure the accuracy
 * of the emitted spectrum #jackson1483_num of about 1%.
 */
double emission_probability_s(double b, double gamma_e, double theta, bool p, double omega) {
    return emission_probability(5, 10'000, b, gamma_e, theta, p, omega);
}

/**
 * The (normalized) energy radiated per unit frequency interval per unit solid angle, computed
 * numerically from #emission_probability_s.
 *
 * An ultrarelativistic particle emit photons in a narrow cone around the direction of the particle
 * velocity. Thus the energy radiated in a certain direction can be readily computed from the
 * energy of the modes of the superconductive virtual box (see #emission_probability). The density
 * of the modes which has a plane-wave component in some certain unit solid angle and unit
 * frequency interval can be easily found from the boundary conditions, that gives
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
        * ( emission_probability_s(b, gamma_e, theta, true,  omega)
          + emission_probability_s(b, gamma_e, theta, false, omega));
}
