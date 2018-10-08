#define MAX(x, y) (((x) > (y)) ? (x) : (y))
#define MIN(x, y) (((x) < (y)) ? (x) : (y))
void solver(int Nx, int Ny)
{
  FILE *vel_u, *vel_v, *vel_p, *itr, *vel_gridx, *vel_gridy;
  vel_u = fopen("vel_u.dat","w");
  vel_v = fopen("vel_v.dat","w");
  vel_p = fopen("pressure.dat","w");
  itr   = fopen("Iteration.dat","w");
  vel_gridx = fopen("xgrid.dat","w");
  vel_gridy = fopen("ygrid.dat","w");
  
  double a,b,c;
  int m,n,o,ans=0;
  double error_dp=1;
  int counter = 0,iteration=0;
  
  double error_p=1,\
    error_u=1,\
    error_v=1;
  
  double total_error_p = 1,\
    total_error_u = 1,\
    total_error_v = 1;
 
  //for(int m =0;m<5;m++)
  while (total_error_v > epsilon && total_error_v > epsilon)
    {
      //printf("hello total error %f\n",total_error_u);
      //velocity and pressure field updated
      for(int i=0;i<=Ny;i++)
	{
	  for(int j=0;j<=Nx;j++)
	    {
	      P_old[i][j] = Q_1[i][j];
	      U_old[i][j] = Q_2[i][j];
	      V_old[i][j] = Q_3[i][j];
	    }
	}
      error_p=1;
      error_u=1;
      error_v=1;
      //printf("Iteration = %d\t%0.12f\t%0.12f\n",iteration,total_error_u,total_error_v);
      m = 0,n=0,o=0;
      while (error_p > epsilon && error_u > epsilon && error_v > epsilon)
	{

	  //////////////////////////////////////////////////////
	  //Co-varient and contra-varient velocities calculated
	  for(int i=0;i<=Ny;i++)
	    {
	      for(int j=0;j<=Nx;j++)
		{
		  U[i][j] = Q_2[i][j] * zeta_x[i][j];
		  V[i][j] = Q_3[i][j] *  eta_y[i][j];
		}
	    }
	  //print2d(U, Nx, Ny);

	  //Calculating Spectral Radius
	  for(int i=1;i<Ny;i++)
	    {
	      for(int j=1;j<Nx;j++)
		{
		  rho_1[i][j] = (fabs(U[i][j]) + sqrt((U[i][j]*U[i][j])  + g11[i][j]))/J[i][j];
		  rho_2[i][j] = (fabs(V[i][j]) + sqrt((V[i][j]*V[i][j])  + g22[i][j]))/J[i][j];
		}
	    }

	  for(int i=1;i<Ny;i++)
	    {
	      for(int j=1;j<Nx;j++)
		{
		  rhox[i][j] = 0.5*(rho_1[i][j+1] + rho_1[i][j]);
		  rhoy[i][j] = 0.5*(rho_2[i+1][j] + rho_2[i][j]);
		}
	    }

	  for(int i=2;i<Ny-1;i++)
	    {
	      for(int j=2;j<Nx-1;j++)
		{
		  Q_half_x_1[i][j] = Q_1[i][j+2] - 3*Q_1[i][j+1] + 3*Q_1[i][j] - Q_1[i][j-1];
		  Q_half_x_2[i][j] = Q_2[i][j+2] - 3*Q_2[i][j+1] + 3*Q_2[i][j] - Q_2[i][j-1];
		  Q_half_x_3[i][j] = Q_3[i][j+2] - 3*Q_3[i][j+1] + 3*Q_3[i][j] - Q_3[i][j-1];

		  Q_half_y_1[i][j] = Q_1[i+2][j] - 3*Q_1[i+1][j] + 3*Q_1[i][j] - Q_1[i-1][j];
		  Q_half_y_2[i][j] = Q_2[i+2][j] - 3*Q_2[i+1][j] + 3*Q_2[i][j] - Q_2[i-1][j];
		  Q_half_y_3[i][j] = Q_3[i+2][j] - 3*Q_3[i+1][j] + 3*Q_3[i][j] - Q_3[i-1][j];
		}
	    }

	  for(int i=2;i<Ny-1;i++)
	    {
	      for(int j=2;j<Nx-1;j++)
		{
		  Diss_x_1[i][j] = epsi*((rhox[i][j]*Q_half_x_1[i][j]) - (rhox[i][j-1]*Q_half_x_1[i][j-1]));
		  Diss_x_2[i][j] = epsi*((rhox[i][j]*Q_half_x_2[i][j]) - (rhox[i][j-1]*Q_half_x_2[i][j-1]));
		  Diss_x_3[i][j] = epsi*((rhox[i][j]*Q_half_x_3[i][j]) - (rhox[i][j-1]*Q_half_x_3[i][j-1]));

		  Diss_y_1[i][j] = epsi*((rhoy[i][j]*Q_half_y_1[i][j]) - (rhoy[i-1][j]*Q_half_y_1[i-1][j]));
		  Diss_y_2[i][j] = epsi*((rhoy[i][j]*Q_half_y_2[i][j]) - (rhoy[i-1][j]*Q_half_y_2[i-1][j]));
		  Diss_y_3[i][j] = epsi*((rhoy[i][j]*Q_half_y_3[i][j]) - (rhoy[i-1][j]*Q_half_y_3[i-1][j]));
		}
	    }

	  for(int i=0;i<=Ny;i++)
	    {
	      for(int j=0;j<=Nx;j++)
		{
		  E_1_compute_x_1[i][j] = U[i][j]/J[i][j];
		  E_1_compute_x_2[i][j] = (Q_2[i][j]*U[i][j])/J[i][j];
		  E_1_compute_x_3[i][j] = (Q_3[i][j]*U[i][j])/J[i][j];

		  E_2_compute_y_1[i][j] = V[i][j]/J[i][j];
		  E_2_compute_y_2[i][j] = (Q_2[i][j]*V[i][j])/J[i][j];
		  E_2_compute_y_3[i][j] = (Q_3[i][j]*V[i][j])/J[i][j];
		}
	    }

	  for(int i=1;i<Ny;i++)
	    {
	      for(int j=1;j<Nx;j++)
		{
		  DE_1_compute_x_1[i][j] = 0.5*(E_1_compute_x_1[i][j+1] - E_1_compute_x_1[i][j-1]) + Diss_x_1[i][j];
		  DE_1_compute_x_2[i][j] = 0.5*(E_1_compute_x_2[i][j+1] - E_1_compute_x_2[i][j-1]) + Diss_x_2[i][j];
		  DE_1_compute_x_3[i][j] = 0.5*(E_1_compute_x_3[i][j+1] - E_1_compute_x_3[i][j-1]) + Diss_x_3[i][j];

		  DE_2_compute_y_1[i][j] = 0.5*(E_2_compute_y_1[i+1][j] - E_2_compute_y_1[i-1][j]) + Diss_y_1[i][j];
		  DE_2_compute_y_2[i][j] = 0.5*(E_2_compute_y_2[i+1][j] - E_2_compute_y_2[i-1][j]) + Diss_y_2[i][j];
		  DE_2_compute_y_3[i][j] = 0.5*(E_2_compute_y_3[i+1][j] - E_2_compute_y_3[i-1][j]) + Diss_y_3[i][j];
		}
	    }
	  //print2d(DE_1_compute_x_1,Nx,Ny);
	  for(int i=0;i<Ny;i++)
	    {
	      for(int j=0;j<Nx;j++)
		{
		  J_half_x[i][j] = 0.5*(J[i][j+1] + J[i][j]);
		  J_half_y[i][j] = 0.5*(J[i+1][j] + J[i][j]);

		  zeta_x_half[i][j] = pow(0.5*(zeta_x[i][j+1] + zeta_x[i][j]),2);
		  eta_y_half[i][j]  = pow(0.5*(eta_y[i+1][j]  + eta_y[i][j]),2);

		  u_11[i][j] = (Q_2[i][j+1] - Q_2[i][j]);
		  u_12[i][j] = (Q_3[i][j+1] - Q_3[i][j]);

		  Ev_1_compute_x_1[i][j] = 0;
		  Ev_1_compute_x_2[i][j] = (zeta_x_half[i][j]*u_11[i][j])/J_half_x[i][j];
		  Ev_1_compute_x_3[i][j] = (zeta_x_half[i][j]*u_12[i][j])/J_half_x[i][j];

		  u_21[i][j] = Q_2[i+1][j] - Q_2[i][j];
		  u_22[i][j] = Q_3[i+1][j] - Q_3[i][j];

		  Ev_2_compute_y_1[i][j] = 0;
		  Ev_2_compute_y_2[i][j] = eta_y_half[i][j]*u_21[i][j]/J_half_y[i][j];
		  Ev_2_compute_y_3[i][j] = eta_y_half[i][j]*u_22[i][j]/J_half_y[i][j];
		}
	    }

	  for(int i=1;i<Ny;i++)
	    {
	      for(int j=1;j<Nx;j++)
		{
		  DEv_1_compute_x_2[i][j] = (Ev_1_compute_x_2[i][j] - Ev_1_compute_x_2[i][j-1])/Re;
		  DEv_2_compute_y_2[i][j] = (Ev_2_compute_y_2[i][j] - Ev_2_compute_y_2[i-1][j])/Re;

		  DEv_1_compute_x_3[i][j] = (Ev_1_compute_x_3[i][j] - Ev_1_compute_x_3[i][j-1])/Re;
		  DEv_2_compute_y_3[i][j] = (Ev_2_compute_y_3[i][j] - Ev_2_compute_y_3[i-1][j])/Re;
		}
	    }

	  for(int i=1;i<Ny;i++)
	    {
	      for(int j=1;j<Nx;j++)
		{
		  grad_P_x[i][j] = J[i][j]*0.5*((Q_1[i][j+1]*zeta_x[i][j+1]/J[i][j+1])\
						- (Q_1[i][j-1]*zeta_x[i][j-1]/J[i][j-1]));
		  grad_P_y[i][j] = J[i][j]*0.5*((Q_1[i+1][j]*eta_y[i+1][j]/J[i+1][j])\
						- (Q_1[i-1][j]*eta_y[i-1][j]/J[i-1][j]));
		}
	    }

	  for(int i=1;i<Ny;i++)
	    {
	      for(int j=1;j<Nx;j++)
		{
		  RHS_compute_1[i][j] = J[i][j]*(-DE_1_compute_x_1[i][j] \
						 -DE_2_compute_y_1[i][j] \
						 +DEv_1_compute_x_1[i][j] \
						 +DEv_2_compute_y_1[i][j]);
		  
		  RHS_compute_2[i][j] = J[i][j]*(-DE_1_compute_x_2[i][j]\
						 -DE_2_compute_y_2[i][j]\
						 +DEv_1_compute_x_2[i][j]\
						 +DEv_2_compute_y_2[i][j]);

		  RHS_compute_3[i][j] = J[i][j]*(-DE_1_compute_x_3[i][j]\
						 -DE_2_compute_y_3[i][j]\
						 +DEv_1_compute_x_3[i][j]\
						 +DEv_2_compute_y_3[i][j]);
		   
		}
	    }
	  //print2d(RHS_compute_3,Nx,Ny);

	  //////////////////////////////////////////////////////////////////

	  for(int i=0;i<=Ny;i++)
	    {
	      for(int j=0;j<=Nx;j++)
		{
		  rho_1[i][j] = 0;
		  rho_2[i][j] = 0;
		}
	    }

	  for(int i=0;i<=Ny;i++)
	    {
	      for(int j=0;j<=Nx;j++)
		{
		  rho_1[i][j] = (fabs(U[i][j]) + sqrt(pow(U[i][j],2) + g11[i][j]))/J[i][j];
		  rho_2[i][j] = (fabs(V[i][j]) + sqrt(pow(V[i][j],2) + g22[i][j]))/J[i][j];

		  dtau_cfl[i][j] = cfl/(J[i][j]*MAX(rho_1[i][j],rho_2[i][j]));
		  dtau_VN[i][j]  = Re*VN/(MAX(g11[i][j],g22[i][j]));
		  dtau[i][j]     = MIN(dtau_cfl[i][j], dtau_VN[i][j]);
		  //printf("%f\t%f\n",MIN(dtau_cfl[i][j], dtau_VN[i][j]),MAX(g11[i][j],g22[i][j]));
		}
	    }

	  for(int i=1;i<Ny;i++)
	    {
	      for(int j=1;j<Nx;j++)
		{
		  ustar[i][j] = Q_2[i][j] + dtau[i][j]*(RHS_compute_2[i][j] - grad_P_x[i][j]);
		  vstar[i][j] = Q_3[i][j] + dtau[i][j]*(RHS_compute_3[i][j] - grad_P_y[i][j]);
		}
	    }

	  //Updating boundary conditions for ustar ?? make it an argument
  
	  /*for(int i=(Ny/2 + 1);i<=Ny;i++)
	    {
	    ustar[i][0] = Ub*(1-(pow(((y[i][0] - 1.5*h)/(0.5*h)),2)));
	    ustar[i][Nx] = ustar[i][Nx-1];
	    }

	    for(int i=0;i<=Ny;i++)
	    {
	    ustar[i][Nx] = ustar[i][Nx-1];
	    }
	    //print2d(ustar,Nx,Ny);*/
	  backward_boundary(Nx,Ny);

	  //Calculating Divergence of ustar
	  for(int i=1;i<Ny;i++)
	    {
	      for(int j=1;j<Nx;j++)
		{
		  diver[i][j] = 0.5*J[i][j]*((ustar[i][j+1]*zeta_x[i][j+1]/J[i][j+1])\
					     - (ustar[i][j-1]*zeta_x[i][j-1]/J[i][j-1])\
					     + (vstar[i+1][j]* eta_y[i+1][j]/J[i+1][j])\
					     - (vstar[i-1][j]* eta_y[i-1][j]/J[i-1][j]));
		}
	    }
	  //print2d(diver,Nx,Ny);

	  for(int i=0;i<=Ny;i++)
	    {
	      for(int j=0;j<=Nx;j++)
		{
		  dp[i][j] = 0;
		  dp1[i][j] = 0;
		}
	    }
	  counter = 0;
	  error_dp =1;

	  while (error_dp > epsilon/10)
	    {
	      //error_dp = 0;
	      for(int i=1;i<Ny;i++)
		{
		  for(int j=1;j<Nx;j++)
		    {
		      dp[i][j] = dp1[i][j];
		    }
		}
	      
	      for(int i=1;i<Ny;i++)
		{
		  for(int j=1;j<Nx;j++)
		    {
		      a = J[i][j]*(    (g11_half[i][j]*dp[i][j+1]/J_half_x[i][j])\
				       + (g11_half[i][j-1]*dp[i][j-1]/J_half_x[i][j-1])\
				       + (g22_half[i][j]*dp[i+1][j]/J_half_y[i][j])\
				       + (g22_half[i-1][j]*dp[i-1][j]/J_half_y[i-1][j]));
	  
		      b = -J[i][j]*(    (g11_half[i][j]/J_half_x[i][j])\
					+ (g11_half[i][j-1]/J_half_x[i][j-1])\
					+ (g22_half[i][j]/J_half_y[i][j])\
					+ (g22_half[i-1][j]/J_half_y[i-1][j]));

		      dp1[i][j] = ((diver[i][j]/dtau[i][j]) - a)/b;
		      
		      error_dp = error_dp + pow((dp1[i][j] - dp[i][j]),2);
		      //printf("error_dp = %0.10f\n",error_dp);

		      error_dp = sqrt(error_dp)/(Nx*Ny);
		    }
		}
	      counter = counter +1;
	      //printf("\nThe error from SIMPLE is %0.15f in iteration = %d\n",error_dp,counter);
	    }
	  /////////////////////////////////////////////////////////////////////////////////end of while loop
	  
	  //Calculating Gradient of delta Pressure
	  for(int i=1;i<Ny;i++)
	    {
	      for(int j=1;j<Nx;j++)
		{
		  dp[i][j] = dp1[i][j];
		  grad_dp_x[i][j] = 0.5*J[i][j]*((dp[i][j+1]*zeta_x[i][j+1]/J[i][j+1]) - (dp[i][j-1]*zeta_x[i][j-1]/J[i][j-1]));
		  grad_dp_y[i][j] = 0.5*J[i][j]*((dp[i+1][j]*eta_y[i+1][j]/J[i+1][j])  - (dp[i-1][j]*eta_y[i-1][j]/J[i-1][j]));
		}
	    }

	  for(int i=1;i<Ny;i++)
	    {
	      for(int j=1;j<Nx;j++)
		{
		  Pressure_new[i][j] = Q_1[i][j] + dp[i][j];
		  u_new[i][j] = ustar[i][j] - dtau[i][j]*grad_dp_x[i][j];
		  v_new[i][j] = vstar[i][j] - dtau[i][j]*grad_dp_y[i][j];
		}
	    }
	  //print2d(dp,Nx,Ny);
	  for(int i=0;i<=Ny;i++)
	    {
	      Pressure_new[i][0] = Pressure_new[i][1];
	      Pressure_new[i][Nx] = Pressure_new[i][Nx-1];
	      u_new[i][Nx] = u_new[i][Nx-1];
	      v_new[i][Nx] = v_new[i][Nx-1];
	    }

	  for(int j=0;j<=Nx;j++)
	    {
	      Pressure_new[0][j] = Pressure_new[1][j];
	      Pressure_new[Ny][j]= Pressure_new[Ny-1][j];
	    }

	  error_p = 0;
	  error_u = 0;
	  error_v = 0;
	  //print2d(Pressure_new,Nx,Ny);
	  for(int i=1;i<Ny;i++)
	    {
	      for(int j=1;j<Nx;j++)
		{
		  error_p = error_p + pow(((Pressure_new[i][j] - Q_1[i][j])),2);
		  error_u = error_u + pow(((u_new[i][j] - ustar[i][j])),2);
		  error_v = error_v + pow(((v_new[i][j] - vstar[i][j])),2);
		}
	    }

	  error_p = sqrt(error_p)/(Nx*Ny);
	  error_u = sqrt(error_u)/(Nx*Ny);
	  error_v = sqrt(error_v)/(Nx*Ny);

	  for(int i=0;i<=Ny;i++)
	    {
	      for(int j=0;j<=Nx;j++)
		{
		  Q_1[i][j] = Pressure_new[i][j];
		  Q_2[i][j] = u_new[i][j];
		  Q_3[i][j] = v_new[i][j];
		  Q_1[0][j] = Q_1[1][j];
		  Q_1[Ny][j] = Q_1[Ny-1][j];
		  Q_1[i][0] = Q_1[i][1];
		  Q_1[i][Nx] = Q_1[i][Nx-1];
		}
	    }
	  ans = m+n+o;
	}
	  
    
    
      ////////////////////////////////////////////////////////////////////////////////////////////////////end of while loop
      //printf("Error_p = %0.10f\tError_u = %0.10f\tError_v = %0.10f\n",error_p,error_u,error_v);
      total_error_p = 0;
      total_error_u = 0;
      total_error_v = 0;

      //printf("%d\t%0.12f\t%0.12f\n",iteration,total_error_u,total_error_v);
      for(int i=0;i<=Ny;i++)
	{
	  for(int j=0;j<=Nx;j++)
	    {
	      total_error_p = total_error_p + fabs((Q_1[i][j] - P_old[i][j]));
	      total_error_u = total_error_u + fabs((Q_2[i][j] - U_old[i][j]));
	      total_error_v = total_error_v + fabs((Q_3[i][j] - V_old[i][j]));
	      //printf("total error is %d\n",abs(Q_2[i][j] - U_old[i][j]));
	    }
	}
      //total_error_p = sqrt(total_error_p)/(Nx*Ny);
      //total_error_u = sqrt(total_error_u)/(Nx*Ny);
      //total_error_v = sqrt(total_error_v)/(Nx*Ny);

      iteration = iteration + 1;
      printf("%d\t%0.12f\t%0.12f\n",iteration,total_error_u,total_error_v);
      fprintf(itr,"%d\t%0.12f\t%0.12f\n",iteration,total_error_u,total_error_v);
    }
  //print2d(U_old,Nx,Ny);
  for(int i=0;i<=Ny;i++)
    {
      for(int j=0;j<=Nx;j++)
	{
	  fprintf(vel_u,"%0.10f\t",Q_2[i][j]);
	  fprintf(vel_v,"%0.10f\t",Q_3[i][j]);
	  fprintf(vel_p,"%0.10f\t",Q_1[i][j]);
	  fprintf(vel_gridx,"%0.10f\t",x[i][j]);
	  fprintf(vel_gridy,"%0.10f\t",y[i][j]);
	}
      fprintf(vel_u,"\n");
      fprintf(vel_v,"\n");
      fprintf(vel_p,"\n");
      fprintf(vel_gridx,"%\n");
	  fprintf(vel_gridy,"\n");
    }
  fclose(vel_u);
  fclose(vel_v);
  fclose(vel_p);
  fclose(itr);
  fclose(vel_gridx);
  fclose(vel_gridy);
  //////////////////////////////////////////////////////////////////////////////////////////////////////end of while loop  
}
