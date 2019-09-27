/**
 * @file
 * @brief Helps to plot particle distributions correctly.
 * @details Particle distribution can have singularity (e.g., @f$ 1 / \sqrt{x} @f$, @f$ x > 0 @f$),
 * and drawing the distribution by counting particles in the bins of the same size is not good in
 * this case. This file proposes @p histogram_1D function which finds the boundaries of bins (of
 * different sizes) which are the most appropriate for the drawing of the particle distribution:
 * near the singularity the bins are thinner than in the regions of small values of the
 * distribution function, and the relative error of the counted distribution function is the same
 * for all the bins.
 */

#include <vector>
#include <algorithm>

/**
 * This function places the bin boundaries such that the number of particles is @p n in every bin,
 * (thus <tt> N % n </tt> particles are dropped from the consideration, where <tt> N = xs.size() </tt>
 * the overall number of particles). The positions of the leftmost and rightmost boundaries are slightly shifted from minimum and
 * maximum of considered x values such that for all the bins the distribution function @f$ f @f$ can be
 * computed from the position of the left and rignt bin boundaries @f$ b_l @f$ and @f$ b_r @f$ as
 * @f\[
 *     f = \frac{n}{b_r - b_l}.
 * @f\]
 * The complexity of this function is @f$ O(N \log N) @f$ (operations).
 * @param n  desired number of particles per bin
 * @param xs coordinates of particles
 * @return position of the bins boundaries
 */
std::vector<double> histogram_1D(size_t n, const std::vector<double>& xs){
    std::vector<double> bs; // boundaries of the bins

    if (xs.size() > 1) { // for a single point the bin width cann't be found
        std::vector<double> ys(xs.size()  - xs.size() %  n);
        std::copy(xs.cbegin(), xs.cend() - xs.size() %  n, ys.begin());

        std::sort(ys.begin(), ys.end());
        // the leftmost boundary
        bs.push_back(ys[0]);

        size_t in = 0;
        for (size_t i = 0; i < ys.size() - 1; ++i) {
            in += 1;
            if (in == n) {
                in = 0;
                bs.push_back(0.5 * (ys[i] + ys[i + 1]));
            }
        }
        // the rightmost boundary
        bs.push_back(ys[ys.size() - 1]);

        /*
         * correction of the leftmost and the rightmost boundaries
         */
        if (bs.size() == 2) { // left and right boundaries exactly at the particle positions
            double d = (bs[1] - bs[0]) / static_cast<double>(n - 1);
            bs[0] -= 0.5 * d;
            bs[1] += 0.5 * d;
        } else if (bs.size() > 2) {
            double d1 = (bs[1] - bs[0]) / (static_cast<double>(n) - 0.5);
            double d2 = (bs[bs.size() - 1] - bs[bs.size() - 2]) / (static_cast<double>(n) - 0.5);
            bs[0]             -= 0.5 * d1;
            bs[bs.size() - 1] += 0.5 * d2;
        }
    }
    return bs;
}
