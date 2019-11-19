/**
 * @file
 * @brief Simple proposal densities for spaces of arbitrary dimension.
 */

#include <tuple>
#include <functional>

/**
 * @defgroup proposal_density Proposal density
 * @brief Simple proposal densities for spaces of arbitrary dimension.
 * @{
 */

/**
 * A C++ functor which provides @ref zip "zipWith::operator()" to zip two tuples element-wise up to
 * a given index with a given function. The first argument is used to return the result.
 */
template<size_t i, typename... Ts>
struct zipWith {
    /**
     * @anchor zip
     * Examples:
     *
     * <tt> zipWith<2, f>((1, 2, 3), (1, 1, 1)) </tt> for <tt>auto f = ()[int x, int y]{return x + 
     * y;}</tt> changes the first tuple to (2, 3, 4). 
     *
     * <tt> zipWith<1, g>((1, 2, 3), (4, 5, 6)) </tt> for <tt>auto g = ()[int x, int y]{return x *
     * y;}</tt> changes the first tuple to (4, 10, 3). 
     */
    void operator()(       auto               f
                   ,       std::tuple<Ts...>& x
                   , const std::tuple<Ts...>& y
                   ) {
        std::get<i>(x) = f( std::get<i>(x), std::get<i>(y) );
        zipWith<i - 1, Ts...> z;
        z(f, x, y);
    }
};

/** Specialization of zipWith for i = 0, to end-up the recursion used in the code. Partial
 * specialization is not allowed for functions, hence #zipWith is a C++ functor, not a function.
 * See <a href = "http://artlang.net/post/c++11-obkhod-elementov-kortezhe-std-tuple/">this
 * post</a> for some details (in Russian).
 */
template<typename... Ts>
struct zipWith<0, Ts...> {
    void operator()(       auto               f
                   ,       std::tuple<Ts...>& x
                   , const std::tuple<Ts...>& y) {
        std::get<0>(x) = f( std::get<0>(x), std::get<0>(y) );
    }
};

/**
 * Returns proposal density,
which for a pair <tt>(rng_state, x)</tt>, advances @p rng_state and @p x such that <tt>x = x +
dx</tt> with @p dx found as described below. Note that if @p x is a tuple, the equation <tt>x = x +
dx</tt> should be considered element-wise, and advance of @p rng_state* is doing for every element
of the tuple @p x.
* @param rng a random number generator which takes and advances its state
* @param a   an amplitute to be used by method @p m
* @param m   a methot to produce @p dx from @p rng_state and @p a. To get @dx evenly distributed in
*            (-a, a) use appropriate <tt>*_cast_with_amplitude</tt> function from rng.hpp (i.e.
*            pm_cast_with_amplitude() if pm_rng() is used). Note also that the proposal
*            distribution should be symmetric in the metropolis algorithm that constrains the
*            choose of @p m.
*/
template<typename R, typename... Ts>
std::function< void( std::pair<R, std::tuple<Ts...>>& ) >
make_proposal_density( void              (rng)(R&)
                     , std::tuple<Ts...> a
                     , auto              m
                     ) {
    return [=](std::pair<R, std::tuple<Ts...>>& rx) {
        auto f = [&rx, rng, m](auto elem_x, auto elem_a){
            rng(rx.first);
            return elem_x + m(rx.first, elem_a);
        };
        zipWith<std::tuple_size<std::tuple<Ts...>>::value - 1, Ts...> z;
        z(f, rx.second, a);
    };
}

/**
 * @}
 */
