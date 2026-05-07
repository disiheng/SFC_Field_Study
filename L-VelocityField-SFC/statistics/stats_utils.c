#include <mpi.h>
#include <stdio.h>
#include <math.h>

#include "../allvars.h"
#include "../proto.h"
#include "stats_utils.h"

double density_statistics(fftw_real *grid)
{
    long long i;
    double delta;
    double local_total = 0.0, total = 0.0;
    double local_var = 0.0, var = 0.0;
    double local_m3 = 0.0, m3 = 0.0, m1 = 0.0, m2 = 0.0;
    double local_min = 1e20, min = 1e20;
    double local_max = -1e20, max = -1e20;
    double local_mean = 0.0, mean = 0.0;
    double local_dm3 = 0.0, dm3 = 0.0;
    double local_dm2 = 0.0, dm2 = 0.0;

    for (i = 0; i < SQR(Ndim) * (nslab_x); i++)
        local_mean += (double) grid[i];
    MPI_Allreduce(&local_mean, &mean, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    mean = mean / CUB(Ndim);

    for (i = 0; i < SQR(Ndim) * (nslab_x); i++) {
        if (local_min > grid[i]) local_min = grid[i];
        if (local_max < grid[i]) local_max = grid[i];
        local_total += (double) grid[i];
        local_var += SQR(grid[i]);
        local_m3 += CUB(grid[i]);
    }

    MPI_Allreduce(&local_m3, &m3, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&local_var, &var, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&local_total, &total, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&local_min, &min, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
    MPI_Allreduce(&local_max, &max, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

    mean = total / CUB(Ndim);
    m1 = mean;
    m2 = var / CUB(Ndim);
    m3 = m3 / CUB(Ndim);

    for (i = 0; i < SQR(Ndim) * (nslab_x); i++) {
        delta = grid[i];
        local_dm2 += SQR(delta);
        local_dm3 += CUB(delta);
    }

    MPI_Allreduce(&local_dm2, &dm2, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&local_dm3, &dm3, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

    dm2 = dm2 / CUB(Ndim);
    dm3 = dm3 / CUB(Ndim);

    if (ThisTask == 0) {
        printf(" ------------------------------------------------\n");
        printf(" Grid statistics\n");
        printf(" Size: %ld\n", Ndim);
        printf(" Local Size: %d starting at %d\n", nslab_x, slabstart_x);
        printf(" Minimun Value: %f\n", min);
        printf(" Maximun Value: %f\n", max);
        printf(" Mean: %f\n", m1);
        printf(" <S^2>: %f\n", m2);
        printf(" <S^3>: %f\n", m3);
        printf(" <(S/<S>-1)^2>: %f\n", dm2);
        printf(" <(S/<S>-1)^3>: %f\n", dm3);
        printf(" Total: %f\n", total);
        printf(" Variance: %f\n", sqrt(var / Ndim - mean * mean));
        printf(" ------------------------------------------------\n");
    }

    return total;
}

double log_density(float *grid)
{
    long long i;
    double local_total = 0.0, total = 0.0;
    double local_log_total = 0.0, log_total = 0.0;
    double log_mean, mean, true_mean;

    true_mean = (double) TotNumHalo / CUB(Ndim);

    for (i = 0; i < SQR(Ndim) * (nslab_x); i++) {
        local_log_total += (double) log(grid[i] / true_mean);
        local_total += (double) (grid[i]);
    }

    MPI_Allreduce(&local_total, &total, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&local_log_total, &log_total, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

    mean = total / CUB(Ndim);
    log_mean = log_total / CUB(Ndim);

    if (ThisTask == 0) {
        printf(" ------------------------------------------------\n");
        printf(" Grid statistics\n");
        printf(" Size: %ld\n", Ndim);
        printf(" Local Size: %d starting at %d\n", nslab_x, slabstart_x);
        printf(" Mean: %f\n", mean);
        printf(" Log-Mean: %f\n", log_mean);
        printf(" True Mean: %f\n", true_mean);
        printf(" ------------------------------------------------\n");
    }

    return log_mean;
}
