#include<vector>
#include<cmath>
#include<random>
#include<complex>

const long double pi = M_PI;
// all quantities in Plank units (\hbar = 1, c = 1, G = 1)
double h = 0.001;// h = H/H_cr - magnetic field Schwinger field(H_cr = m^2 c^3 /(e\hbar) ) ratio
const double m = 1.0; // lepton/proton mass in electron masses
const double gamma0 = 100.0;// = energy / (m_{particle} c^2)
double omega_cr = 0.45*h * gamma0*gamma0;//
using namespace std;
typedef complex<double> dcomp;
dcomp i(0,1);

double ind(double h,double omega){
  double kappa = h*omega;
  if (kappa <= 2.01523 ){
    return 0.22;
  } else{
    return 0.56 / pow(kappa,4.0/3.0);
  }
}

dcomp dipole_jE(double omega, double theta,double t){

  double omega0 = 1;// \hbar / (m_e * c^2) units
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
  // rotating particle
  double beta = sqrt(1 - 1/(gamma0 * gamma0));
  double omega0 =  h /(gamma0 * m * m);
  double n = 1 + 1.0 * h * h * ind(h,omega) / ( pi * 137.0);//index of refraction (Re part)

  return cos(theta) * exp( (-i * omega* t + i * n * (omega / omega0) * beta * sin(omega0*t) * cos(theta) ) * gamma0/(gamma0 - omega / m) );
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

    default_random_engine gen{static_cast<long unsigned int>(time(NULL))};
    uniform_real_distribution<double> omega_dist(0.01,gamma0);
    uniform_real_distribution<double> angle_dist(-20.0/gamma0,20.0/gamma0);
    uniform_real_distribution<double> r(0,1);
    double* x = new double[size_t (2*n)];
    x[0] = 100;// even - ang. freq. (omega)
    x[1] = 0.0;// odd - angle
    double form_time = 1.0*m*m / h;// coherence time /formation time
    double t1 = - form_time / 2.0;
    double t2 = form_time / 2.0;
    double tmp = c1s(t1,t2,x[0],x[1]);
    for(unsigned int j = 1; j < n; ++j){
      double omega = omega_dist(gen);
      double angle = angle_dist(gen);
      double s = c1s(t1,t2,omega,angle);
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

  double* metropolis_spectrum_sy(size_t n){
    // spectrum for a given angle!
    default_random_engine gen{static_cast<long unsigned int>(time(NULL))};
    uniform_real_distribution<double> omega_dist(0.1,gamma0);//
    uniform_real_distribution<double> r(0,1);
    double* x = new double[n];
    x[0] = 1.0;//seed value
    double angle = 0.0;
    double form_time = 1.0*m*m / h;// coherence time /formation time
    double t1 = - form_time / 2.0;
    double t2 = form_time / 2.0;
    double tmp = c1s(t1,t2,x[0],angle) * x[0];// omega * |C_m|^2
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

  double* metropolis_radpattern_sy(size_t n){
    // radiation pattern for a given frequency
    default_random_engine gen{static_cast<long unsigned int>(time(NULL))};
    uniform_real_distribution<double> distribution(-20.0/gamma0,20.0/gamma0);
    uniform_real_distribution<double> r(0,1);
    double* x = new double[n];
    x[0] = 1.0;//seed value
    double form_time = 1.0*m*m / h;// coherence time /formation time
    double t1 = - form_time / 2.0;
    double t2 = form_time / 2.0;
    double tmp = c1s(t1,t2,omega_cr, x[0]);
    for(unsigned int j = 1; j<n; ++j){
      double tmp_x = distribution(gen);
      double s = c1s(t1,t2,omega_cr,tmp_x);
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
