#include<vector>
#include<cmath>
#include<random>
#include<complex>
const long double pi = 3.141592653589793238;
using namespace std;
typedef complex<double> dcomp;
dcomp i(0,1);
dcomp dipole_jE(double omega, double theta,double t){
  //double c = 3 * pow(10.0,10.10);
  double omega0 = 2 * pi;// 1 Hz
  //double lambda0 = 2 * pi * c / omega0 ;
  double a = 0.01;// length / lambda0
  //double e = 4.8 * pow(10.0,-10.0)
  //double dlength = a * lambda0 -- dipole length
  // 2 pi e omega0 * dlength = (2 pi)^2 * 4.8 * 3 = 568.5
  //double k = omega / c ;
  double kd = a * 2 * pi * omega / omega0 ;
  return - 568.5 * sin(theta) * cos(omega0 * t) * exp( -i*omega * t + i * kd * cos(theta) * sin(omega0 * t) );
}

double c1d(double t1,double t2,double omega, double theta){
  const int m = 1000;
  dcomp sum(0,0);
  for (unsigned int j = 0;j <= m;j++){
    sum += dipole_jE(omega, theta, t2 * j / m);
  }
  return pow(abs(sum),2.0);
}

dcomp sync_jE(double omega, double theta,double t){
  double gamma = 1000;
  double beta = 1 - 1 / (2 * gamma * gamma);
  double h = 0.01;// h = H/H_cr
  double omega0 = h / gamma;
  return -90.477 * beta * cos(theta) * exp( -i * omega * t + i * (omega / omega0) * beta * sin(omega0*t) * cos(theta) );
  //return -90.477 * cos(theta) * exp(-i * omega* t *(1 - beta*cos(theta) ) - i * omega* omega0*omega0*t*t*t / 6.0 );
}

double c1s(double t1,double t2,double omega, double theta){
  const int n = 1000;
  dcomp sum(0,0);
  for (unsigned int j = 0;j <= n;j++){
    sum += sync_jE(omega, theta, t2 * j / n);
  }
  return pow(abs(sum),2.0);
}

extern "C" {
  // for now dipol only
  double* metropolis_test(size_t n){
    default_random_engine gen{static_cast<long unsigned int>(time(NULL))};
    uniform_real_distribution<double> omega_dist(0,4*pi);//for dipole
    uniform_real_distribution<double> angle_dist(0,pi);
    uniform_real_distribution<double> r(0,1);
    double* x = new double[size_t (2*n)];
    x[0] = 0.1;// even - ang. freq. (omega)
    x[1] = 0.1;// odd - angle
    double tmp = c1d(0,10,x[0],x[1]);
    for(unsigned int j = 1; j < n; ++j){
      double omega = omega_dist(gen);
      double angle = angle_dist(gen);
      double s = c1d(0,10,omega,angle);
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
    default_random_engine gen{static_cast<long unsigned int>(time(NULL))};
    uniform_real_distribution<double> omega_dist(1000, 20000 );// omega_cr = h*gamma*gamma
    uniform_real_distribution<double> r(0,1);
    double* x = new double[n];
    x[0] = 2000;//- ang. freq. (omega)
    double angle = 0.0;
    double t2 = 1000;
    double tmp = c1s(0,t2,x[0],angle);// t2 ~ 1/h
    for(unsigned int j = 1; j < n; ++j){
      double omega = omega_dist(gen);
      double s = c1s(0,t2,omega,angle);
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
