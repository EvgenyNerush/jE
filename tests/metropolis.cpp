/* In this test *n* x-coordinates with *target_f* distribution are generated, then the average
 * value of x, x^2 and x^3 are computed and compared with theoretical values.
 * This test takes about 2.2 seconds with Xeon X5550 processor.
 */

#include <iostream>
#include <stdexcept>
#include <cmath>
#include <vector>
#include <tuple>
#include <functional>

// Park--Miller RNG
const uint64_t pm_randmax = 0x7FFFFFFFull; // approx 2*10^9
uint64_t pm_rng(uint64_t r) {
    return 48271ull * r % pm_randmax;
}

// Float (or double) in (0, 1), computed from value returned by pm_rng
template<typename F>
F to_floating_01(uint64_t r) {
    return static_cast<F>(r) / static_cast<F>(pm_randmax);
}

// target function used in this test
double target_f(double x){
    if (x > 0 and x <= 1) {
        return 1 / sqrt(x);
    } else {
        return 0;
    }
}

// target function used in this test
double alt_target_f(std::tuple<double> x){
    if (std::get<0>(x) > 0 and std::get<0>(x) <= 1) {
        return 1 / sqrt(std::get<0>(x));
    } else {
        return 0;
    }
}

// Element-wise operation on two tuples up to the given index, e.g. zipWith<1, (*)>((1, 2, 3), (4,
// 5, 6)) changes the first tuple to (4, 10, 3). One need partial specialization to end-up the
// recursion, but partial specialization is not allowed for functions. Hence zipWith is a C++
// functor.  (See, in Russian, http://artlang.net/post/c++11-obkhod-elementov-kortezhe-std-tuple/)
template<size_t i, typename... Ts>
struct zipWith {
    void operator()( auto f
             ,       std::tuple<Ts...>& x
             , const std::tuple<Ts...>& y) {
        std::get<i>(x) = f( std::get<i>(x), std::get<i>(y) );
        zipWith<i - 1, Ts...> z;
        z(f, x, y);
    }
};

// Variant for i = 0, to end-up the recursion
template<typename... Ts>
struct zipWith<0, Ts...> {
    void operator()( auto f
             ,       std::tuple<Ts...>& x
             , const std::tuple<Ts...>& y) {
        std::get<0>(x) = f( std::get<0>(x), std::get<0>(y) );
    }
};

// This function gets a random-number generator *rng* (e.g, void function(uint64_t& rng_state)), an
// amplitute *a* and a method *m* (which produces *dx* from *rng_state* and *a*), and returns a
// walk function, which for a given pair (rng_state, x) advances the *rng_state* and yields x = x +
// dx. Here (+), *m* and advance of *rng_state* are called for every element of *x* if *x* and *a*
// are tuples. Note that the walk function should be symmetric in the metropolis algorithm that
// constrains the choose of *m*.
template<typename R, typename... Ts>
std::function< void( std::pair<R, std::tuple<Ts...>>& ) >
make_walk_f( std::function<void (R&)> rng
           , std::tuple<Ts...> a
           , auto m) {
    return [=](std::pair<R, std::tuple<Ts...>>& rx){
        // qwe; auto f = [&rx, rng, m](auto elem_a, auto elem_x){
        auto f = [&rx, rng, m](auto elem_x, auto elem_a){
            rng(rx.first);
            // qwe;//std::cout << elem_a << '\t';
            return elem_x + m(rx.first, elem_a);
        };
        zipWith<std::tuple_size<std::tuple<Ts...>>::value - 1, double> z;
        z(f, rx.second, a);
    };
}

// Simple walk function of (rng_state, x)
std::pair<uint64_t, double> walk_f(std::pair<uint64_t, double> rx) {
    uint64_t r = pm_rng(rx.first);
    double d = to_floating_01<double>(r);
    //std::cout << rx.second << '\t' << 0.5 * (2 * d - 1) << '\n';
    return std::make_pair(r, rx.second + 0.5 * (2 * d - 1));
}

double m(uint64_t r, double a) {
    return a * (2 * to_floating_01<double>(r) - 1);
}
void myrng(uint64_t& r) {
    r = pm_rng(r);
}
//std::function<void(std::pair<uint64_t, double>&)> alt_walk_f
std::function<void(std::pair<uint64_t, std::tuple<double>>&)> alt_walk_f
    = make_walk_f( //[](uint64_t& r){ pm_rng(r); }
            std::function<void(uint64_t&)>(myrng)
                 , std::make_tuple<double>(0.1)
                 //, [](uint64_t r, double a){ return a * (2 * to_floating_01<double>(r) - 1); } );
                 , m );

double mean_x_pow(std::vector<double> xs, int pw){
    double acc = 0;
    for (auto x: xs) {
        acc += pow(x, pw);
    }
    return acc / xs.size();
}

double alt_mean_x_pow(std::vector<std::tuple<double>> xs, int pw){
    double acc = 0;
    for (auto x: xs) {
        acc += pow(std::get<0>(x), pw);
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

// Random points generated with metropolis algorithm
template<typename R, typename X>
std::vector<X> metropolis
            ( std::function<void(R&)> rng     // RNG used to accept/decline new points
            , uint64_t s0                     // seed for this RNG
            , std::function<double(R)> to_01  // converts RNG values to doubles evenly distributed
                                              // in (0, 1)
            , std::function
              <void(std::pair<R, X>&)> walk_f // walk function
            , std::pair<R, X> rx0             // (seed for RNG of walk_f, initial point)
            , std::function
              <double(X)> target_f            // target function
            , size_t n                        // number of points to generate
            ) { 
    std::vector<X> x(n);
    x[0] = rx0.second;
    std::pair<R, X> rx = rx0;
    R s = s0;
    for (size_t i = 1; i < n; ++i) {
        alt_walk_f(rx);
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

int mult(int a, int b) {
    return a * b;
}

int main() {
    /*
    std::tuple<int, int, int> x(1,2,3);
    std::tuple<int, int, int> y(4,5,6);
    zipWith<1, int, int, int> z;
    z(mult, x, y);
    std::cout << std::get<0>(x) << '\n';
    std::cout << std::get<1>(x) << '\n';
    std::cout << std::get<2>(x) << '\n';
    */
    size_t n = 16'000'000; // should be much less than pm_randmax
    //double acceptable_accuracy = 10 / sqrt(static_cast<double>(n));
    // !!!
    //qwe; QWE !!!!! absolute, not relative!
    double acceptable_accuracy = 4 / sqrt(static_cast<double>(n));
    //std::vector<double> xs = metropolis(target_f, n, walk_f, std::make_pair(12345, 0.5), 123);
    std::vector<std::tuple<double>> xs_t
        = metropolis<uint64_t, std::tuple<double>>( std::function<void(uint64_t&)>(myrng)
                    , 12345ul
                    , std::function<double(uint64_t)>(to_floating_01<double>)
                    , alt_walk_f
                    , std::make_pair(54321ul, std::tuple<double>(0.5))
                    , std::function<double(std::tuple<double>)>(alt_target_f)
                    , n);
    /*std::vector<double> xs(n);
    for (size_t i = 0; i < n; ++i) {
        xs[i] = std::get<0>(xs_t[i]);
    }
    */
    std::vector<int> pws {1, 2, 3};
    std::vector<double> theor;
    std::vector<double> num;
    std::vector<double> relative_errs;
    bool acc = true;
    for (auto pw: pws) {
        double theor_v = theor_mean_x_pow(pw);
        double num_v = alt_mean_x_pow(xs_t, pw);
        theor.push_back(theor_v);
        num.push_back(num_v);
        double relative_err = fabs(num_v - theor_v);// qwe; !!! / theor_v;
        relative_errs.push_back(relative_err);
        acc = acc and relative_err < acceptable_accuracy;
    }
    if (acc) {
        std::cout << "metropolis test: \x1b[32mpassed\x1b[0m\n";
    } else {
        std::cout << "metropolis test: \x1b[1;31mfailed\x1b[0m\n";
        std::cout << "acceptable_accuracy = " << acceptable_accuracy << '\n';
        std::cout << "power" << '\t'
                  << "theor value" << '\t'
                  << "num value" << '\t'
                  << "relative error" << '\n';
        for (size_t i = 0; i < pws.size(); ++i) {
            std::cout << pws[i] << '\t'
                      << theor[i] << '\t'
                      << num[i] << '\t'
                      // qwe; << (relative_errs[i] < acceptable_accuracy) << '\t'
                      << relative_errs[i] << '\n';
        }
    }
    return 0;
}
