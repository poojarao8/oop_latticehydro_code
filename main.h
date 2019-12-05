#ifndef main_h
#define main_h

#include <iostream>
#include <cmath>
#include <mpi.h>
#include "HYPRE_struct_ls.h"
#include "Grid.h"
#include "Field.h"
#include "TimeIntegration.h"

#define X 0
#define Y 1
#define Z 2
#define YZ 0
#define ZX 1
#define XY 2 //XY, YZ and ZX are the oriented planes (i.e. 2-simplex) at each lattice point.
#define NDIMS 3
#define NGUARD 2

using namespace std;

const int nsteps = 1;
const double eta = 100.0;
const double coef = 2.5;
const double nu = 0.001;

#endif
