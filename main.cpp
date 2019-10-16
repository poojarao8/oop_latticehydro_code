#include "main.h"

using namespace std;

int main()
{
  int nsteps = 1;

  Grid coarse_grid(32, 32, 32);

  Field velocity(3, &coarse_grid);
  velocity.initialize();
  velocity.update_bdry('PERIODIC');

  Field pressure(1, &coarse_grid);
  pressure.initialize();
  pressure.update_bdry('PERIODIC');

  return 0;
}
