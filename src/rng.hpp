/**
 * @file
 * @brief Random number generators in form appropriate for #metropolis and @p cast_to_01 functions
 * for them.
 */

/**
 * RAND_MAX for Park-Miller RNG; approximately @f$ 2 \times 10^9 @f$.
 */
const uint64_t pm_randmax = 0x7FFFFFFFull; // 

/**
 * Park--Miller random number generator
 */
void pm_rng(uint64_t& r) {
    r = 48271ull * r % pm_randmax;
}

/**
 * Function to compute evenly distributed values in (0, 1) from values returned by @p pm_rng.
 */
double pm_cast_to_01(uint64_t r) {
    return static_cast<double>(r) / static_cast<double>(pm_randmax);
}

/**
 * Function to compute values evenly distributed in <tt>(-a, a)</tt>, from values returned by
 * @pm_rng.
 */
double pm_cast_with_amplitude(uint64_t r, double a) {
    return a * (2 * pm_cast_to_01(r) - 1);
}
