#include <vector>
#include <cmath>
using namespace std;

// for target function testing
double f1(double x) {
    return exp(-100 * (x - 0.5) * (x - 0.5)) + pow((1 - x), 10) * (1 + tanh(50 * x));
}

double f2(double x) {
    return exp(-500 * (x + 0.5) * (x + 0.5)) + exp(-500 * (x - 0.5) * (x - 0.5));
}

const double a = 4; // 1.8
// random-walk function
double g(double x) {
    int r = rand();
    if (r < RAND_MAX / 2) {
        return x - a * (double) r / (double) RAND_MAX;
    } else {
        return x + a * ((double) r / (double) RAND_MAX - 0.5);
    }
}

// thin with wings, one-third of points in 0.5+-0.1 for b = 0.4. b = 1 corr to uniform from -a to
// a.
const double b = 0.5;
double g2(double x) {
    double r = (double) rand() / (double) RAND_MAX;
    r = b * r + (1 - b) * 0.5 * (pow((2 * r - 1), 3) + 1); // not uniform, in [0, 1) (?)
    return x + a * (r - 0.5);
}

extern "C" {
double* metropolis(size_t n) {
    srand(123);
    double* x = new double[n];
    x[0] = 0.25;
    for (size_t i = 1; i < n; ++i) {
        double y = g2(x[i - 1]);
        double alpha = f2(y) / f2(x[i - 1]);
        if (alpha >= 1) {
            x[i] = y;
        } else {
            int r = rand();
            if ((double) r / (double) RAND_MAX < alpha) {
                x[i] = y;
            } else {
                x[i] = x[i - 1];
            }
        }
    }
    return x;
}
}
