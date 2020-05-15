#include <vector>
#include <fstream>
#include <iostream>
#include <cmath>
#include "../src/radiation.hpp"
#include <thread>
using namespace std;

// linspace function NOT similar to one in python3
std::vector<double> linspace(double a, double b, size_t m) {
    double h = (b - a) / static_cast<double>(m);
    std::vector<double> xs(m);
	typename std::vector<double>::iterator x;
    double val;
    for (x = xs.begin(), val = a + 0.5 * h; x != xs.end(); ++x, val += h) {
        *x = val;
    }
    return xs;
}

// Note that all normalizations follow those of ../src/radiation.hpp
int main() {

    ofstream vacuum;
    vacuum.open("dipole_data/vacuum.txt"); // distribution d^2I/d\omega d\theta, computed with
                                          // delta_n DO NOT taken into account
    ofstream ref_ind;
    ref_ind.open("dipole_data/ref_ind.txt"); // distribution d^2I/d\omega d\theta, computed with
                                            // delta_n TAKEN into account

    ofstream parameters; // contains frequency and angle ranges, etc.
    parameters.open("sinus_data/parameters.txt");

    double gamma_e = 1e5;
    double b = 1;
	double chi = b * gamma_e;
	double mass = 1; // particle mass ratio to the electron mass
	double coeff = 2; // > 1
	double width_L = b * coeff; // normalized laser width
	double omega_L = 1.0 / b; // normalized laser frequency 

    // angle varies from 0 to theta_b
    const double theta_b = 0.5 * M_PI;
    size_t m = 87; // number of dots in angle;
    size_t n = 99; // number of dots in frequency ranges
    // angle and omega for the second plot

    vector<double> thetas( linspace(0, theta_b, m) );
    double omega_e = gamma_e / b;
    double omega_left = 0.0;
    double omega_right = omega_e;
    vector<double> omegas( linspace(omega_left, omega_right, n) );

    auto I = function<double(double, double, double, double)> (
        [=](double ri, double b, double theta, double omega) {
            auto f = bks_dipole_td(ri, mass, b, width_L, omega_L, gamma_e);
            return omega * b * f(make_tuple(theta, omega)) / (2 * M_PI * r(gamma_e));
        }
    );

    parameters << b           << '\n'
               << m           << '\n'
               << n           << '\n'
               << gamma_e     << '\n'
               << omega_left  << '\n'
               << omega_right << '\n'
               << omega_L     << '\n'
               << width_L     << '\n'
               ;
    for (auto omega: omegas) {
        for (auto theta: thetas) {
            vacuum  << I( 1,  b, theta, omega) << '\n';
//            ref_ind << I(ri,  b, theta, omega) << '\n';
        }
    }

    vacuum.close();
    ref_ind.close();
    parameters.close();

    return 0;
}
