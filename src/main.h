#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <stdlib.h>
//#include <omp.h>
#include <algorithm>

#include "index_conversion.h"
#include "coarse_to_fine.h"
#include "fine_to_coarse.h"
#include "fluid_terms.h"
#include "lattice_calculus.h"
#include "solver.h"
#include "read_write.h"
#include "norm.h"
#include "visit_plot.h"
#include "initialize_var.h"
#include "parse_params.h"

#define RESET   "\033[0m"
#define RED     "\033[31m" //color code for red text
#define YZ 0
#define ZX 1
#define XY 2 //XY, YZ and ZX are the oriented planes (i.e. 2-simplex) at each lattice point.
#define X 0
#define Y 1
#define Z 2  //X, Y and Z are the oriented segments (i.e. 1-simplex) at each lattice point.
#define OUTNAME "Dfd_" //Output file for computed value. Python can read this by numpy package with  np.fromfile command.
#define OUTNAMEC "Coarse_"
#define L 64
#define W 64
#define H 64//Length, width and height of the fine grid. Coarse grid will be half.
#define MU 2.5

const int LL=L/2;
const int WW=W/2;
const int HH=H/2;
const int SCL_ARR_SIZE=L*W*H;
const int COARSE_SCL_ARR_SIZE=LL*WW*HH;
const int ARR_SIZE=3*L*W*H;
const int COARSE_ARR_SIZE=3*LL*WW*HH;
const int OCTALSIZE=LL*WW*HH;
const double factor=4.0*OCTALSIZE;
