#include "main.h"

using namespace std;

int main()
{
  Grid coarse_grid(32, 32, 32);
  Grid fine_grid(64, 64, 64);

  Field velocity(3, &coarse_grid, 0);
  velocity.initialize();

  Field pressure(1, &coarse_grid, 0);
  pressure.initialize();

  TimeIntegration RK4;

  return 0;
}
