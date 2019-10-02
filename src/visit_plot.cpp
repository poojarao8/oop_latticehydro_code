#include <math.h>
#include "visit_writer.h"
#include "visit_plot.h"
#include <sstream>
#include <iomanip>
#include <cstring>
#include <iostream>

using namespace std;

void visit_plot(const char* filename, int L, int W, int H, double* pressure, double* velocity, double phy_time, int smoothing_count)
{
    // Set up storage
    int nCellsX = L;
    int nCellsY = W;
    int nCellsZ = H;

    // Write out the mesh and the arrays.
    int dims[] = {nCellsX, nCellsY, nCellsZ};// cell centered
    int vardims[] = {1, 3}; // Two scalars
    int centering[] = {1, 1}; // node centered, cell centered
    const char * const varnames[] = { "L2_grad_velocity", "velocity" };
    double *arrays[] = { pressure, velocity };

    write_regular_mesh(filename, 0, dims, 2, vardims, centering,
                      varnames, arrays, phy_time, smoothing_count);

}


void output_file(int ts, double nu, char* outfile)
{
  const char* basename = "out";
  const char* extname = ".vtk";

  std::stringstream convert_nu;
  convert_nu << "_nu" << std::setfill('0') << std::setw(4) << nu;

  std::stringstream convert_ts;
  convert_ts << "_ts" << std::setfill('0') << std::setw(6) << ts;

  strcpy(outfile,basename); //copy string one into the outfile
  strcat(outfile,convert_nu.str().c_str());
  strcat(outfile,convert_ts.str().c_str());
  strcat(outfile,extname); //append string two to the outfile
  std::cout << outfile << std::endl;
}
