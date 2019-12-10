#include <iostream>
#include <vector>
#include <iomanip>
#include <fstream>
#include "../src/radiation.hpp"
using namespace std;

int main() {
    ofstream vacuum_violation;
    ofstream vacuum_perturb;
    ofstream medium_violation;
    ofstream medium_perturb;
    ofstream other_data;
    vacuum_violation.open("vv.txt");
    vacuum_perturb.open("vp.txt");
    medium_violation.open("mv.txt");
    medium_perturb.open("mp.txt");
    other_data.open("od.txt");
    double gamma_e = 1000;
    double b = pow(1.0 / alpha, 1.5) / gamma_e;

    double oc = omega_c(gamma_e);
    const int n = 6;

    vector<double> thetas;
    for (int i = -n/2; i < n/2-1;  ++i) {
        thetas.push_back( static_cast<double>(i) / gamma_e );
        other_data << static_cast<double>(i) / gamma_e << "_";
    }
    thetas.push_back( static_cast<double>(n/2-1) / gamma_e );
    other_data << static_cast<double>(n/2-1) / gamma_e << "\n";

    vector<double> omegas;
    for (double x = 1; x < n; x += 1) {
        omegas.push_back(x * 0.1 * oc);
        other_data << x * 0.1 * oc << "_";
    }
    omegas.push_back( static_cast<double>(n) * 0.1 * oc);
    other_data << static_cast<double>(n) * 0.1 * oc;

    for (int i = 0; i < n ; i++) {
        for (int j = 0; j < n-1; j++) {

            vacuum_violation << bks_synchrotron_emission_probability(1, 1, b, gamma_e, thetas[j], omegas[i]) << "_";
            vacuum_perturb << bks_synchrotron_emission_probability(1, 1, 0.01 * b, gamma_e, thetas[j], omegas[i]) << "_";
            medium_violation << bks_synchrotron_emission_probability(vacuum_refractive_index(b,omegas[i]), 1, b, gamma_e, thetas[j], omegas[i]) << "_";
            medium_perturb << bks_synchrotron_emission_probability(vacuum_refractive_index(0.01 * b,omegas[i]), 1, 0.01 * b, gamma_e, thetas[j], omegas[i]) << "_";
        }
        vacuum_violation << bks_synchrotron_emission_probability(1, 1, b, gamma_e, thetas[n-1], omegas[i]) << "_";
        vacuum_perturb << bks_synchrotron_emission_probability(1, 1, 0.01 * b, gamma_e, thetas[n-1], omegas[i]) << "_";
        medium_violation << bks_synchrotron_emission_probability(vacuum_refractive_index(b,omegas[i]), 1, b, gamma_e, thetas[n-1], omegas[i]) << "_";
        medium_perturb << bks_synchrotron_emission_probability(vacuum_refractive_index(0.01 * b,omegas[i]), 1, 0.01 * b, gamma_e, thetas[n-1], omegas[i]) << "_";
    }
    vacuum_violation.close();
    vacuum_perturb.close();
    medium_violation.close();
    medium_perturb.close();
    other_data.close();
    return 0;
}
