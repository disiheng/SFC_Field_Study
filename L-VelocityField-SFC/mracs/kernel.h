#ifndef KERNEL_H
#define KERNEL_H

typedef enum {
    KERNEL_SHELL,
    KERNEL_SPHERE,
    KERNEL_GAUSSIAN,
    KERNEL_DUAL_RING
} KernelType;

double window_shell(double R, double kx, double ky, double kz);
double window_sphere(double R, double kx, double ky, double kz);
double window_gaussian(double R, double kx, double ky, double kz);
double window_dual_ring(double R, double theta, double kx, double ky, double kz);
double kernel_eval(KernelType type, double R, double theta, double kx, double ky, double kz);

#endif
