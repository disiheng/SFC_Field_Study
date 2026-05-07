
#ifndef HISTOGRAM_H
#define HISTOGRAM_H

#include <stdlib.h>
#include <gsl/gsl_histogram.h>
#include <gsl/gsl_histogram2d.h>

void save_histogram(gsl_histogram *h, char *filename);
void save_histogram2d(gsl_histogram2d *h, char *filename);

#endif
