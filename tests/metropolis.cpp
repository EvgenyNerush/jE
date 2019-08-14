/* In this test *n* x-coordinates with *target_f* distribution are generated, then the average
 * value of x, x^2 and x^3 are computed and compared with theoretical values.
 */

#include <iostream>
#include <vector>
#include <stdexcept>
#include <cmath>

double target_f(double x){
    if (x > 0 and x <= 1) {
        return 1 / sqrt(x);
    } else {
        return 0;
    }
}

// Park--Miller RNG
const uint64_t pm_randmax = 0x7FFFFFFFull;
uint64_t pm_rng(uint64_t r) {
    return 48271ull * r % pm_randmax;
}
// to_Type( p / pm_randmax)...
// walk function of (rng_state, x)
std::pair<uint64_t, double> walk_f(std::pair<uint64_t, double> rx) {
    uint64_t r = pm_rng(rx.first);
    double d = static_cast<double>(r) / static_cast<double>(pm_randmax);
    //std::cout << rx.second << '\t' << 0.5 * (2 * d - 1) << '\n';
    return std::make_pair(r, rx.second + 0.5 * (2 * d - 1));
}

double mean_x_pow(std::vector<double> xs, int pw){
    double acc = 0;
    for (auto x: xs) {
        acc += pow(x, pw);
    }
    return acc / xs.size();
}

// int_0^1 x^(p - 1/2) dx / \int_0^1 x^(-1/2) dx
double theor_mean_x_pow(int pw){
    if (pw > -0.5) {
        return 0.5 / (pw + 0.5);
    } else {
        std::invalid_argument("pw shoud be > -1/2");
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////
//
//uint32_t pm_state = 42;
//
// random numbers generated with metropolis algorithm
std::vector<double> metropolis( double target_f(double)
                              , size_t n
                              , std::pair<uint64_t, double> walk_f(std::pair<uint64_t, double>)
                              , std::pair<uint64_t, double> rx0 // ...
                              , uint64_t s0 // ...
                              ) { 
    std::vector<double> x(n);
    x[0] = rx0.second;
    std::pair<uint64_t, double> rx = rx0;
    uint64_t s = s0;
    for (size_t i = 1; i < n; ++i) {
        rx = walk_f(rx);
        double acceptance = target_f(rx.second) / target_f(x[i - 1]);
        // избежать срыва конвейера?
        if (acceptance < 1) {
            s = pm_rng(s);
            if (static_cast<double>(s) / static_cast<double>(pm_randmax) >= acceptance) {
                rx.second = x[i - 1];
            }
        }
        x[i] = rx.second;
    }
    return x;
}
///////////////////////////////////////////////////////////////////////////////////////////////////

int main() {
    size_t n = 1000000;
    double acceptable_accuracy = 5 / sqrt(static_cast<double>(n));
    std::vector<double> xs = metropolis(target_f, n, walk_f, std::make_pair(12345, 0.5), 123);
    /*for (auto x: xs) {
        std::cout << x << '\n';
    }
    */
    std::vector<int> pws {1, 2, 3};
    std::vector<double> relative_errs;
    bool acc = true;
    for (auto pw: pws) {
        double relative_err
            = fabs(mean_x_pow(xs, pw) - theor_mean_x_pow(pw)) / theor_mean_x_pow(pw);
        // std::cout << mean_x_pow(xs, pw) << '\t' << theor_mean_x_pow(pw) << '\n';
        relative_errs.push_back(relative_err);
        acc = acc and relative_err < acceptable_accuracy;
    }
    if (acc) {
        std::cout << "metropolis test: \x1b[32mpassed\x1b[0m\n";
    } else {
        std::cout << "metropolis test: \x1b[1;31mfailed\x1b[0m\n";
        std::cout << "acceptable_accuracy = " << acceptable_accuracy << '\n';
        std::cout << "power" << '\t' << "relative error" << '\n';
        for (size_t i = 0; i < pws.size(); ++i) {
            std::cout << pws[i] << '\t' << relative_errs[i] << '\n';
        }
    }
    return 0;
}
