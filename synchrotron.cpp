#include<vector>
#include<iostream>
#include<cmath>
#include<random>
#include<complex>

const double pi = M_PI;
// all quantities in Plank units (\hbar = 1, c = 1, G = 1)
double h = 0.0001;// h = H/H_cr - magnetic field Schwinger field(H_cr = m^2 c^3 /(e\hbar) ) ratio
const double m = 1.0; // lepton/proton mass in electron masses
double gamma0 = 100.0;// = energy / (m_{particle} c^2)
double beta = sqrt( 1 - 1/(gamma0 * gamma0) );
double omega0 = h/(m*m*gamma0);
double tcoh = 1.0*m*m / h;// coherence time /formation time
double omega_cr = 0.3*h * gamma0*gamma0;
using namespace std;
typedef complex<double> dcomp;
dcomp i(0,1);

double ind(double h1,double omega){
  double kappa = h1*omega;
  double dn = 0;
  if (kappa <= 2.01523 ){
    dn = 0.22;
  } else{
    dn = 0.56 / pow(kappa,4.0/3.0);
  }
  return 1 + 1.0 * h1 * h1 * dn /(pi * 137.0);
}

double tau1(double omega){
  return (gamma0 - omega/m) / (gamma0 * omega * abs(1 - beta ) );
}
double tau2(double omega){
  return pow( 6 * (gamma0 - omega/m) / (omega *gamma0 * omega0*omega0) , 1.0/3.0 );
}

dcomp sync_jE1(double omega,double theta,double t){
  //double phi = t/tau1(omega) + pow(t / tau2(omega),3.0);
  double phi = ( omega*t*(1 - beta) + omega*omega0*omega0*t*t*t/6 ) /(1.0 - omega /(m* gamma0) );
  //double phi =  omega*( t - beta * sin(omega0*t) * cos(theta) / omega0 ) /(1.0 - omega /(m* gamma0) );
  return exp( -i*phi )*sin(omega0 * t);
}

double c1s(double t1,double t2,double omega, double theta){
  // calculating |C_1|^2
  double step = 0;
  if (tau1(omega) > tau2(omega)){
    step = 0.1*tau2(omega);
  } else {
    step = 0.1*tau1(omega);
  }

  size_t n = (t2 - t1)/step;
  dcomp sum(0,0);
  for (unsigned int j = 0;j <= n;j++){
    sum += sync_jE1(omega, theta, t1 + step*j)*step;
  }
  return pow(abs(sum),2.0);
}

extern "C" {
  double* metropolis_spectrum_sy(size_t n){
    // spectrum for a given angle!
    default_random_engine gen{static_cast<long unsigned int>(time(NULL))};
    uniform_real_distribution<double> omega_dist(0.00001,40);//
    uniform_real_distribution<double> r(0,1);
    double* x = new double[n];
    x[0] = 1;//seed value
    double angle = 0.0;
    double form_time = tau2(x[0]);
    double tmp = c1s( -form_time / 2.0,form_time / 2.0,x[0],angle);// |C_m|^2
    for(unsigned int j = 1; j < n; ++j){
      double omega = omega_dist(gen);

      form_time = 10*tau2(omega);
      //form_time = tau2(omega) * pow( tau2(omega)/tau1(omega) , 0.5);
      double s = c1s( -form_time / 2.0,form_time / 2.0,omega,angle);
      if (s >= tmp){
        x[j] = omega;
        tmp = s;
      } else if (s > tmp * r(gen)){
        x[j] = omega;
        tmp = s;
      } else {
        x[j] = x[j-1];
      }
    }
    return x;
  }
}
