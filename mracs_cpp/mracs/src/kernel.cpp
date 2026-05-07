#include"csmain.h"
#include<cmath>

#ifndef TWOPI
#define TWOPI M_PI*2
#endif



// shell window function: sin(2Pi*kR)/(2Pi*kR)
double WindowFunction_Shell(double R, double ki, double kj, double kk)
{
    double k = sqrt(ki * ki + kj * kj + kk * kk);
    if(k == 0) return 1.;
    double Phase = TWOPI * k * R;
    return sin(Phase)/Phase;
}

// sphere window function: 4Pi*[sin(2Pi*kR)-2Pi*kRcos(2Pi*kR)]/(2Pi*k)^3
double WindowFunction_Sphere(double R, double ki, double kj, double kk)
{
    double k = sqrt(ki * ki + kj * kj + kk * kk);
    if(k == 0) return 4.*M_PI*pow(R,3)/3.;
    double Phase = TWOPI * k * R;
    return (sin(Phase)-Phase*cos(Phase))/(2*M_PI*M_PI*pow(k,3));
}

// Gaussian window function: e^(-(1/2)(kR)^2)
double WindowFunction_Gaussian(double R, double ki, double kj, double kk)
{
    double k = sqrt(ki * ki + kj * kj + kk * kk);
    double Phase = TWOPI * k * R;
    return pow(1/M_E,Phase*Phase/2);
}

// dual-ring window function: J_0(2Pi*kRsin(a)sin(b))*cos(2Pi*kRcos(a)cos(b))
double WindowFunction_Dual_Ring(double R, double theta, double ki, double kj, double kk)
{
    return (std::cyl_bessel_j(0,TWOPI*sin(theta)*R*sqrt(ki*ki+kj*kj)))*cos(TWOPI * R * cos(theta) * kk);
}