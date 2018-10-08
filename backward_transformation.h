double **x_zeta, **y_eta, **zeta_x, **eta_y;
double **G, **J;
double **g11, **g22;
double **g11_half, **g22_half;
void backward_transformation(int Nx, int Ny)
{
  initialize(&x_zeta,Nx,Ny,0);
  initialize(&y_eta,Nx,Ny,0);
  initialize(&zeta_x,Nx,Ny,0);
  initialize(&eta_y,Nx,Ny,0);
  initialize(&G,Nx,Ny,0);
  initialize(&J,Nx,Ny,0);
  initialize(&g11,Nx,Ny,0);
  initialize(&g22,Nx,Ny,0);
  initialize(&g11_half,Nx,Ny,0);
  initialize(&g22_half,Nx,Ny,0);
  
  for(int i=1;i<Ny;i++)
    {
      for(int j=1;j<Nx;j++)
	{
	  x_zeta[i][j] = (x[i][j+1] - x[i][j-1])/2;
	  y_eta [i][j] = (y[i+1][j] - y[i-1][j])/2;
	  x_zeta[0][j] = (x[0][j+1]   - x[0][j-1])/2;
	  x_zeta[Ny][j] =(x[Ny][j+1]  - x[Ny][j-1])/2;
	  y_eta[i][0]   =(y[i+1][0]   - y[i-1][0])/2;
	  y_eta[i][Nx]  =(y[i+1][Nx]  - y[i-1][Nx])/2;
	}
    }
  

  for(int i=0;i<=Ny;i++)
    {
      x_zeta[i][0]  = x[i][1] - x[i][0];
      x_zeta[i][Nx] = x[i][Nx]- x[i][Nx-1];
    }

  for(int i=0;i<=Nx;i++)
    {
      y_eta[0][i] = y[1][i] - y[0][i];
      y_eta[Ny][i] =y[Ny][i]- y[Ny-1][i];
    }
  
  for(int i=0;i<=Ny;i++)
    {
      for(int j=0;j<=Nx;j++)
	{
	  G[i][j] = x_zeta[i][j] * y_eta[i][j];
	  J[i][j] = 1/G[i][j];
	}
    }
  
    for(int i=0;i<=Ny;i++)
    {
      for(int j=0;j<=Nx;j++)
	{
	  zeta_x[i][j] = J[i][j]*y_eta[i][j];
	  g11[i][j]    = zeta_x[i][j]*zeta_x[i][j];
	  eta_y[i][j]  = J[i][j]*x_zeta[i][j];
	  g22[i][j]    = eta_y[i][j]*eta_y[i][j];
	}
    }
  
  for(int i=1;i<Ny;i++)
    {
      for(int j=1;j<Nx;j++)
	{
	  g11_half[i][j] = 0.5*(g11[i][j+1] + g11[i][j]);
	  g22_half[i][j] = 0.5*(g22[i+1][j] + g22[i][j]);
	}
    }
 }
