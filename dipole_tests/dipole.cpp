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
	ofstream ch_angle;
	ch_angle.open("dipole_data/ch_angle.txt"); // cherenkov angle as a function of frequency

    ofstream parameters; // contains frequency and angle ranges, etc.
    parameters.open("dipole_data/parameters.txt");

    double gamma_e = 1.0E+8; // particle energy at the begining
	double beta_e = 1 - 0.5 / (gamma_e * gamma_e);
    double b = 0.1; // normalized field
	double mass = 1836; // particle mass ratio to the electron mass
	double omega_L = 0.9; // laser frequency or 1 / a0.
	double x_L = 100; 	// half size of the field region

    // angle varies from 0 to theta_b
	const double theta_right = 0.5 * M_PI;
	const double theta_left  = 0.0;
    size_t m = 27; // number of dots in angle;
    size_t n = 29; // number of dots in frequency ranges
    // angle and omega for the second plot

    vector<double> thetas( linspace( theta_left, theta_right, m) );
    double omega_e = gamma_e / b;
    double omega_left = 0.0;
    double omega_right = omega_e * 0.9;
    vector<double> omegas( linspace(omega_left, omega_right, n) );

    auto I = function<double(double, double, double, double)> (
        [=](double ri, double b, double theta, double omega) {
            auto f = bks_dipole_td(ri, mass, b, x_L, omega_L,gamma_e);
            return omega * b * f(make_tuple(theta, omega)) / (2 * M_PI * r(gamma_e));
        }
    );

    parameters << b           << '\n'
               << m           << '\n'
               << n           << '\n'
               << gamma_e     << '\n'
               << omega_left  << '\n'
               << omega_right << '\n'
               << x_L     	  << '\n'
               << omega_L     << '\n'
               << theta_right << '\n'
               << theta_left  << '\n'
               ;
    for (auto omega: omegas) {
		double ri = vacuum_refractive_index(b, omega);
		ch_angle << acos(1.0 / (ri * beta_e) ) << '\n';
        for (auto theta: thetas) {
            vacuum  << I(1.0, b, theta, omega) << '\n';
            ref_ind << I(ri,  b, theta, omega) << '\n';
        }
    }

    vacuum.close();
    ref_ind.close();
	ch_angle.close();
    parameters.close();

    return 0;
}
