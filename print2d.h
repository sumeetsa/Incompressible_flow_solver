void print2d(double **X, int Nx, int Ny)
{
  //printf("hello");
  for(int i=0;i<=Ny;i++)
    {
      for(int j=0;j<=Nx;j++)
	{
	  printf("%0.4f\t",X[i][j]);
	}
      printf("\n");
    }
}
