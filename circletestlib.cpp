#include <vector>
#include <cmath>
#include <complex>
using namespace std;

// for description, see circletest.py

const double radius = 30;
const double g = 30; // Lorentz factor, gamma
const double omega_c = g * g * g / radius;

const double v = sqrt(1 - 1 / (g * g));

// число точек для вычисления интеграла
//const long n = 5000;

const double a = 0.2 * omega_c; // amplitude of the walk function
double walk_f(double x) {
    double r = (double) rand() / (double) RAND_MAX;
    double y = x + a * (2 * r - 1);
    if (y < 0) { // это не нарушает требование симметрии of walk function: вероятность генерации x'
                 // из x раына вероятности генерации x из x'
        y = - y;
    }
    return y;
}

// Park--Miller RNG
const uint64_t pm_randmax = 0x7FFFFFFFull;
uint32_t pm_rng(uint32_t state) {
    return static_cast<uint32_t>(48271ull * static_cast<uint64_t>(state) % pm_randmax);
}
uint32_t pm_state = 42;

/*double phi_max(double omega) {
    // интегрируем в пределах порядка phi_2 (т.е. tau_\perp)
    //return 2 * pow(6 * M_PI / (omega * radius), 1/3);
    //return 2 / radius;
    //return 4 / g;
    double phi1 = omega_c / (omega * g);
    double phi2 = 1 / g * pow(omega_c / omega, 1/3);
    if (phi1 < phi2) {return 40 * phi1;} else {return 3 * phi2;}
}*/

// target function
double target_f(double omega) {
    // вычисляем интеграл с помощью метода Монте-Карло в интервале от -phi_max до phi_max
    complex<double> c1(0,0);
    double phim;
    // число точек для вычисления интеграла
    long n;
    double phi1 = omega_c / (omega * g);
    double phi2 = 1 / g * pow(omega_c / omega, 1/3);
    /*if (phi1 > phi2) {
        phim = 3 * phi2;
        n = 1000;
    } else {
        if (omega < 2 * omega_c) {
        phim = 3 * phi2;// 3 * phi2; // 40 * phi1;
        n = 10000;// * phi2 / phi1;
        } else {
            phim = 1e-9;
            n = 1;
        }
    }*/
    phim = 3 * phi2;
    /*if (phim > M_PI) {
        phim = M_PI;
    }*/
    n = 1000;
    // we use Park-Miller RNG here
    for (long j = 0; j < n; ++j) {
        pm_state = pm_rng(pm_state);
        double phi = phim * (2 * (double) pm_state / (double) pm_randmax - 1);
        complex<double> i(0,1);
        c1 += sin(phi) * exp(-i * omega * radius * (phi - v * sin(phi)));
    }
    double res;
    if (abs(c1) > 0 * phim * sqrt(n)) {
        res = norm(c1 * phim / (double) n);
    } else {
        res = 0;
    }
    return res;
}

extern "C" {
// random numbers generated with metropolis algorithm
double* rnd(size_t n) {
    srand(1235);
    double* x = new double[n];
    x[0] = 0.3 * omega_c;
    for (size_t i = 1; i < n; ++i) {
        double candidate = walk_f(x[i - 1]);
        double acceptance = target_f(candidate) / target_f(x[i - 1]);
        if (acceptance >= 1 || (double) rand() / (double) RAND_MAX < acceptance) {
            x[i] = candidate;
        } else {
            x[i] = x[i - 1];
        }
    }
    return x;
}
}
