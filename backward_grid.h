double **x, **y;
void backward_grid(int Nx, int Ny)
{
  
  //int Nx = 5, Ny = 5;
  double dx =1, dy = 1;
  int Lx = 1;
  int Ly = 1;
  FILE *f;
  f = fopen("grid.dat","w");
  initialize(&x, Nx, Ny, 0);
  initialize(&y, Nx, Ny, 0);

  
  x[0][0] = 0;
  y[0][0] = 0;

  dx = (double)Lx/Nx;
  dy = (double)Ly/Ny;
  
  for(int i=0;i<=Ny;i++)
    {
      for(int j=0;j<=Nx;j++)
	{
	  x[i][j] = x[0][0] + j*dx;
	  y[i][j] = y[0][0] + i*dy;
	  fprintf(f,"%0.5f\t%0.5f\n",x[i][j],y[i][j]);
	}
    }
  fclose(f);
}
