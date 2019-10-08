/**
 * @file
 * @brief Formulas for synchrotron emission
 * @details Formulas describing the emission of photons by ultrarelativistic electrons rotating in
 * magnetic field @f$ B @f$. Here time is normalized to the radiation formation time 
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
 * The curvature radius of the electron trajectory
 * @param gamma_e the electron Lorentz factor
 */
double r(double gamma_e) {
    return sqrt(gamma_e * gamma_e - 1);
}

/**
 * The critical frequency of the classical synchrotron emission, as defined in Eq. (14.85) from
 * [Jackson J.D., Classical Electrodynamics, Wiley, 1962].
 * @param gamma_e the electron Lorentz factor
 */
double omega_c(double gamma_e) {
    return 3 * pow(gamma_e, 3) / r(gamma_e);
}

/**
 * The energy radiated per unit frequency interval per unit solid angle, as defined in Eq. (14.83)
 * from [Jackson J.D., Classical Electrodynamics, Wiley, 1962].
 * @param b       the normalized magnetic field strength
 * @param gamma_e the electron Lorentz factor
 * @param theta   angle in the plane perpendicular to the normal vector of the trajectory; @p theta
 *                = 0 is the direction of the tangent to the trajectory
 * @param omega   > 0, normalized frequency of the emitted photon
 */
double jackson1483(double b, double gamma_e, double theta, double omega) {
    double a = 1 / (gamma_e * gamma_e) + theta * theta;
    double xi = omega * r(gamma_e) * pow(a, 3.0/2);
    double k13 = std::cyl_bessel_k(1.0/3, xi);
    double k23 = std::cyl_bessel_k(2.0/3, xi);
    return alpha * b * pow(omega * r(gamma_e) * a, 2) * (k23 * k23 + pow(theta * k13, 2) / a);
}

/**
 * The probability of the photon emission, computed numerically.
 *
 * The emitted field in the general
 * case can be decomposed by modes:
 * @f[
 *     \mathbf{E} = \sum_m C_m \mathbf{E}_m, \quad \mathbf{B} = \sum_m C_m \mathbf{B}_m,
 * @f]
 * where (complex) modes are orthogonal:
 * @f[
 *     \int_V \mathbf{E}_m \mathbf{E}_l^* \, dV = \int_V \mathbf{B}_m \mathbf{B}_l^* \, dV = 4 \pi
 *     \hbar \omega_m \delta_{ml}.
 * @f]
 * Hence, the energy of the emitted field is
 * @f[
 *     I = \sum_m \hbar \omega_m |C_m|^2,
 * @f]
 * and @f$ |C_m|^2 @f$ can be interpreted as the probability of emission of the photon of mode @f$
 * m @f$. It can be found from Maxwell's equations:
 * @f[
 *     C_m = \frac{1}{\hbar \omega_m} \int_t \int_V \mathbf{j E}_m \, dV \, dt,
 * @f]
 * where @f$ \mathbf{j} = e \mathbf{v} \delta(\mathbf{r} - \mathbf{r}(t)) @f$ is the current of the
 * emitting electron.
 *
 * Here @f$ C_m @f$ is computed with midpoint rule of integration (very similar to the trapezoidal
 * rule) for orthogonal modes of a virtual rectangular box of size @f$ L_x \times L_y \times L_z
 * @f$:
 * @f[
 *     \mathbf{E}_m = \sqrt{ \frac{4 \pi \hbar \omega_m}{V} } \mathbf{e}_m \exp \left[ -i \omega_m
 *     t + i \mathbf{k}_m \mathbf{r}(t) \right],
 * @f]
 * with @f$ \mathbf{e}_m @f$ the polarization vector, @f$ e_m^2 = 1 @f$, the wave vector @f$ (k_x,
 * k_y, k_z) = (2 \pi m_x / L_x, 2 \pi m_y / L_y, 2 \pi m_z / L_z) @f$, and the generalized index
 * @f$ m = (p, m_x, m_y, m_z) @f$, @f$ p @f$ is the polarization index and @f$ m_j = \pm 1, \pm 2,
 * ... @f$. The volume @f$ V = L_x L_y L_z @f$ is omitted in the computation because it anyway
 * vanishes from the expression for @f$ I @f$ if one goes from @f$ \sum_m @f$ to @f$ \int \, d^3
 * \mathbf{k} @f$.
 * @param a       value which shrinks/enlarges the limits of integration, i.e. the ratio of the
 *                upper (lower) limit of integration to the "transverse" scale of the exponent to
 *                integrate
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
    double tau_perp = gamma_e * pow(6 / (omega * r(gamma_e)), 1/3.0);
    // the limits of the numerical integration are [-tau_m, tau_m]
    double tau_m = a * tau_perp;
    // step of the integration
    double dt = 2 * tau_m / static_cast<double>(n);
    // C_m
    std::complex<double> c = 0i;
    if (p == true) {
        for (long long int j = 0; j < n; ++j) {
            double t = -tau_m + dt * (0.5 + static_cast<double>(j));
            c -= dt * r(gamma_e) / gamma_e
                * sin(t / gamma_e) * exp(-1i * omega * t
                                        + 1i * omega * r(gamma_e) * sin(t / gamma_e) * cos(theta));
        }
    } else {
        for (long long int j = 0; j < n; ++j) {
            double t = -tau_m + dt * (0.5 + static_cast<double>(j));
            c += dt * r(gamma_e) / gamma_e * sin(theta)
                * cos(t / gamma_e) * exp(-1i * omega * t
                                        + 1i * omega * r(gamma_e) * sin(t / gamma_e) * cos(theta));
        }
    }
    return 4 * M_PI * alpha / omega * norm(c);
}

/**
 * Simplified version of #emission_probability, where @p a and @p n are set to ensure the accuracy
 * of about 1%.
 */
double emission_probability_s(double b, double gamma_e, double theta, bool p, double omega) {
    return emission_probability(4, 10'000, b, gamma_e, theta, p, omega);
}
