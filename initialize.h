void initialize(double ***matrix, int Nx, int Ny,int p)
{
  //printf("Error from Initialize\n");
  *matrix = (double **)calloc(Ny,sizeof(double *));

  for(int i=0;i<=Ny;i++)
    {
      (*matrix)[i] = (double *)calloc(Nx, sizeof(double));
    }

  /*for(int i=0;i<=Ny;i++)
  {
    for(int j=0;j<=Nx;j++)
      {
	matrix[i][j] = 0;
      }
  }*/
  if(p==1)
    {
      print2d(*matrix,Nx,Ny);
    }
  
}
