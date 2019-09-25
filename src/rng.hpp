/**
 * @file
 * @brief Random number generators in form appropriate for @metropolis and @p cast_to_01 functions
 * for them.
 */

/**
 * RAND_MAX for Park-Miller RNG; approximately @f$ 2 \times 10^9 @f$.
 */
const uint64_t pm_randmax = 0x7FFFFFFFull; // 

/** Park--Miller RNG
 */
void pm_rng(uint64_t& r) {
    r = 48271ull * r % pm_randmax;
}

/** Function to compute evenly distributed values in (0, 1) from value returned by @p pm_rng.
 */
double pm_cast_to_01(uint64_t r) {
    return static_cast<double>(r) / static_cast<double>(pm_randmax);
}
