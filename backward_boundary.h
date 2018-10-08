double Ub = 1;
double h = 1;

void backward_boundary(int Nx, int Ny)
{
  /*for(int i=(Ny/2 +1);i<=Ny;i++)
    {
      Q_2[i][0] = Ub*(1-(pow(((y[i][0] - 1.5*h)/(0.5*h)),2)));
      Q_2[Ny][0] = 0;
      ustar[i][0] = Ub*(1-(pow(((y[i][0] - 1.5*h)/(0.5*h)),2)));
      ustar[Ny][0] = 0;
      u_new[i][0] = Ub*(1-(pow(((y[i][0] - 1.5*h)/(0.5*h)),2)));
      u_new[Ny][0] = 0;
    }
  for(int i=0;i<=Ny;i++)
    {
      Q_1[i][Nx] = Q_1[i][Nx-1];
      Q_2[i][Nx] = Q_2[i][Nx-1];
      Q_3[i][Nx] = Q_3[i][Nx-1];
      ustar[i][Nx] = ustar[i][Nx-1];
      u_new[i][Nx] = u_new[i][Nx-1];
      vstar[i][Nx] = vstar[i][Nx-1];
      v_new[i][Nx] = v_new[i][Nx-1];
    }
    //print2d(u_new,Nx,Ny);*/

  for(int j=0;j<=Nx;j++)
    {
      Q_2[Ny][j] = 1;
      ustar[Ny][j] = 1;
      u_new[Ny][j] = 1;
      }//*/
}
