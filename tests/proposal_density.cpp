/* Testing that make_proposal_density with pm_rng yield function which then generates evenly
 * distributed numbers.
 */

#include <iostream>
#include <cmath>
#include "../src/rng.hpp"
#include "../src/proposal_density.hpp"

using namespace std;

int main() {
    auto pd = make_proposal_density
                  ( pm_rng
                  , std::make_tuple<double, double>(1, 1)
                  , pm_cast_with_amplitude
                  );
    uint64_t rng_state = 123;
    size_t   n         = 1'000'000;

    // We generate n pairs (a, b) with a and b evenly distributed in (-1, 1), then we compute the
    // second momentum (sigma) for a and b; the theoretical value for evenly distributed numbers in
    // (-1, 1) is 1/sqrt(3)
    double sigma_x = 0;
    double sigma_y = 0;
    std::tuple<double, double> x;
    auto state = std::make_pair(rng_state, x);
    for (size_t i = 0; i < n; ++i) {
        get<0>(state.second) = 0;
        get<1>(state.second) = 0;
        pd(state);
        sigma_x += pow(get<0>(state.second), 2);
        sigma_y += pow(get<1>(state.second), 2);
    }
    sigma_x = sqrt(sigma_x / static_cast<double>(n));
    sigma_y = sqrt(sigma_y / static_cast<double>(n));

    double theor_sigma = 1 / sqrt(3.0);
    double err     = 0.001;

    if (  fabs(sigma_x - theor_sigma) < err
       && fabs(sigma_y - theor_sigma) < err) {
        cout << "proposal_density test: \x1b[32mpassed\x1b[0m\n";
    } else {
        cout << "proposal_density test: \x1b[1;31mfailed\x1b[0m\n";
        cout << "theor sigma = "      << theor_sigma << '\n'
             << "computed sigma_x = " << sigma_x     << '\n'
             << "computed sigma_y = " << sigma_y     << '\n';
    }

    return 0;
}
