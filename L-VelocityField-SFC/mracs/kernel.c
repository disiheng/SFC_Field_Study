#include <math.h>
#include "kernel.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

double window_shell(double R, double kx, double ky, double kz)
{
    double k = sqrt(kx * kx + ky * ky + kz * kz);
    if (k == 0.0) return 4.0 * M_PI * R * R;
    return 4.0 * M_PI * R * R * sin(k * R) / (k * R);
}

double window_sphere(double R, double kx, double ky, double kz)
{
    double k = sqrt(kx * kx + ky * ky + kz * kz);
    if (k == 0.0) return 4.0 * M_PI * R * R * R / 3.0;
    double kR = k * R;
    return 4.0 * M_PI * (sin(kR) - kR * cos(kR)) / (k * k * k);
}

double window_gaussian(double R, double kx, double ky, double kz)
{
    double k2 = kx * kx + ky * ky + kz * kz;
    double R2 = R * R;
    return pow(2.0 * M_PI * R2, 1.5) * exp(-k2 * R2 / 2.0);
}

double window_dual_ring(double R, double theta, double kx, double ky, double kz)
{
    double k = sqrt(kx * kx + ky * ky + kz * kz);
    if (k == 0.0) return 0.0;
    double cost = (kx * sin(theta) + ky * cos(theta)) / k;
    double kR = k * R;
    return 2.0 * cos(kR * cost) * 4.0 * M_PI * R * R * sin(kR) / (kR);
}

double kernel_eval(KernelType type, double R, double theta, double kx, double ky, double kz)
{
    switch (type) {
        case KERNEL_SHELL:      return window_shell(R, kx, ky, kz);
        case KERNEL_SPHERE:     return window_sphere(R, kx, ky, kz);
        case KERNEL_GAUSSIAN:   return window_gaussian(R, kx, ky, kz);
        case KERNEL_DUAL_RING:  return window_dual_ring(R, theta, kx, ky, kz);
        default:                return 0.0;
    }
}
