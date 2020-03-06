#include <math.h>
#include <iostream>
#include <sstream>
#include <cstring>
#include <cstdlib>
#include <iomanip>
#include "visit_writer.h"
#include "visit_plot.h"


void visit_plot(const char* filename, int* nblocks, int* procGridLow, double* pressure, double* velocity, int* coords)
{
    // Set up storage
    int nCellsX = nblocks[0] + 2*(coords[0] % 2);
    int nCellsY = nblocks[1] + 2*(coords[1] % 2);
    int nCellsZ = nblocks[2] + 2*(coords[2] % 2);

    // Write out the mesh and the arrays.
    int dims[] = {nCellsX, nCellsY, nCellsZ};// cell centered
    int vardims[] = {1, 3}; // Two scalars
    int centering[] = {1, 1}; // node centered, cell centered
    const char * const varnames[] = {"pressure","velocity"};
    double *arrays[] = {pressure, velocity};

    write_regular_mesh(filename, 0, dims, procGridLow, 2, 
                       vardims, centering, varnames, arrays);
}

void output_file(int rank, int ts, char* outfile)
{
  const char* basename = "out";
  const char* extname = ".vtk";

  std::stringstream convert_ts;
  convert_ts << "_ts" << std::setfill('0') << std::setw(6) << ts << "_" << std::setw(3) << rank;

  strcpy(outfile,basename); //copy string one into the outfile
  strcat(outfile,convert_ts.str().c_str());
  strcat(outfile,extname); //append string two to the outfile
  std::cout << outfile << std::endl;
}
