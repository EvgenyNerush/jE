#include<vector>
#include<cmath>
#include<random>
#include<complex>
const long double pi = M_PI;
using namespace std;
typedef complex<double> dcomp;
dcomp i(0,1);

dcomp dipole_jE(double omega, double theta,double t){

  double omega0 = 2 * pi;// 1 Hz dipol frequency
  double a = 0.01;// a = dipole_length / lambda0
  //k = omega / c
  double kd = a * 2 * pi * omega / omega0 ;//
  return sin(theta) * cos(omega0 * t) * exp( -i*omega * t + i * kd * cos(theta) * sin(omega0 * t) );
}

double c1d(double t1,double t2,double omega, double theta){
  const int m = 1000;
  dcomp sum(0,0);
  for (unsigned int j = 0;j <= m;j++){
    sum += dipole_jE(omega, theta, t1 + (t2 - t1) * j / m);
  }
  return pow(abs(sum),2.0);
}

dcomp sync_jE(double omega, double theta,double t){
  // all quantities in Plank units
  // m_e = 1
  double gamma = 100;
  double beta = sqrt(1 - 1/(gamma * gamma));
  double h = 0.001;// h = H/H_cr
  double omega0 = h / gamma;
  return cos(theta) * exp( (-i * omega * t + i * (omega / omega0) * beta * sin(omega0*t) * cos(theta) ) * gamma/(gamma - omega) );
}

double c1s(double t1,double t2,double omega, double theta){
  // calculating |C_1|^2
  const int n = 1000;
  dcomp sum(0,0);
  for (unsigned int j = 0;j <= n;j++){
    sum += sync_jE(omega, theta, t1 + (t2 - t1) * j / n);
  }
  return pow(abs(sum),2.0);
}

extern "C" {

  double* metropolis_integration(size_t n){
    // doesn't work now !
    default_random_engine gen{static_cast<long unsigned int>(time(NULL))};
    uniform_real_distribution<double> omega_dist(10,20000);//for dipole
    uniform_real_distribution<double> angle_dist(-pi/2,pi/2);
    uniform_real_distribution<double> r(0,1);
    double* x = new double[size_t (2*n)];
    x[0] = 100;// even - ang. freq. (omega)
    x[1] = pi / 2;// odd - angle
    double tmp = c1s(0,1000,x[0],x[1]);
    for(unsigned int j = 1; j < n; ++j){
      double omega = omega_dist(gen);
      double angle = angle_dist(gen);
      double s = c1s(0,1000,omega,angle);
      if (s >= tmp){
        x[2*j] = omega;
        x[2*j+1] = angle;
        tmp = s;
      } else if (s > tmp * r(gen)){
        x[2*j] = omega;
        x[2*j+1] = angle;
        tmp = s;
      } else {
        x[2*j] = x[2*(j-1)];
        x[2*j+1] = x[2*j-1];
      }
    }
    return x;
  }

  double* metropolis_spectrum(size_t n){
    // spectrum for a given angle
    default_random_engine gen{static_cast<long unsigned int>(time(NULL))};
    uniform_real_distribution<double> omega_dist(0.1, 10);// omega_cr = 1.5 h*gamma*gamma
    uniform_real_distribution<double> r(0,1);
    double* x = new double[n];
    x[0] = 2;//- ang. freq. (omega)
    double angle = 0.0;
    double t1 = -15000;
    double t2 = 15000;// t2 ~ 1/h
    double tmp = c1s(t1,t2,x[0],angle) * x[0];// |C_m|^2 * omega
    for(unsigned int j = 1; j < n; ++j){
      double omega = omega_dist(gen);
      double s = c1s(t1,t2,omega,angle) * omega;
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

  double* metropolis_radpattern(size_t n){
    // radiation pattern for a given frequency
    default_random_engine gen{static_cast<long unsigned int>(time(NULL))};
    uniform_real_distribution<double> distribution(0 ,pi);
    uniform_real_distribution<double> r(0,1);
    double* x = new double[n];
    x[0] = 0.1;
    double tmp = c1d(0,10,2*pi,x[0]);
    for(unsigned int j = 1; j<n; ++j){
      double tmp_x = distribution(gen);
      double s = c1d(0,10,2*pi,tmp_x);
      if (s >= tmp){
        x[j] = tmp_x;
        tmp = s;
      } else if (s > tmp * r(gen)){
        x[j] = tmp_x;
        tmp = s;
      } else {
        x[j] = x[j-1];
      }
    }
    return x;
  }
}
