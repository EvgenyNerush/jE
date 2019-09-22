/**
 * @file
 * @brief Abstract <a href =
 * "https://en.wikipedia.org/wiki/Metropolis%E2%80%93Hastings_algorithm">Metropolis algorithm</a>
 * for symmetric proposal distribution.
*/

#include <iostream>
#include <stdexcept>
#include <cmath>
#include <vector>
#include <tuple>
#include <functional>

/**
 * A std::vector of random points distributed accordingly to the target distribution @p target_f.
 */
template<typename R, typename X>
std::vector<X> metropolis
            ( std::function<void(R&)> rng     ///< RNG used to accept/decline new points
            , R seed                          ///< seed for this RNG
            , std::function<double(R)> to01   ///< function to convert RNG values to doubles evenly
                                              ///  distributed in (0, 1)
            , std::function<void(std::pair<R, X>&)>
              proposal_density                ///< proposal distribution function (jumping
                                              ///  distribution) which updates the point and its
                                              //   own RNG state
                                              //   */
            , std::pair<R, X> rx0             ///< pair (seed for RNG of @p proposal_density,
                                              ///  initial point for it)
            , std::function
              <double(X)> target_distribution ///< distribution to generate
            , size_t n                        ///< number of points to generate
            ) { 
    std::vector<X> x(n);
    x[0] = rx0.second;
    std::pair<R, X> rx = rx0;
    R s = seed;
    for (size_t i = 1; i < n; ++i) {
        proposal_dansity(rx);
        // qwe; std::cout << std::get<0>(rx.second) << '\n';
        double acceptance = target_f(rx.second) / target_f(x[i - 1]);
        if (acceptance < 1) {
            rng(s);
            if (to_01(s) >= acceptance) {
                rx.second = x[i - 1];
            }
        }
        x[i] = rx.second;
    }
    // qwe?
    return std::move(x);
}
