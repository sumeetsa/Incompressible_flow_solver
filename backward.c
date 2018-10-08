#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include "print2d.h"
#include "initialize.h"
#include "backward_grid.h"
#include "backward_transformation.h"
#include "backward_parameters.h"
#include "backward_initialize.h"
#include "backward_boundary.h"
#include "solver.h"



int main()
{
  int Nx = 100;
  int Ny = 100;
  backward_grid(Nx,Ny);
  backward_transformation(Nx,Ny);
  backward_parameters(Nx,Ny);
  backward_initialize(Nx,Ny);
  backward_boundary(Nx,Ny);
  solver(Nx,Ny);
  //print2d(Q_2,Nx,Ny);
  return 0;
}
