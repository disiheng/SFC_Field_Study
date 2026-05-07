
#ifndef CORREL_H
#define CORREL_H

#include "density.h"

void correl( char fname[], Grid grid, int foldnum, float boxsize); 

float periodic (float x, float currBox);
int wrapn (int ii, int ngrid);

#define TAG_GRID      20
#define TAG_POSGRID      340

#endif
