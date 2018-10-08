double **Q_1, **Q_2, **Q_3;
double **P_old, **U_old, **V_old;
double **dp;
double **ustar, **vstar;
double **u_new, **v_new;
double **vorticity;
double **rho_1, **rho_2;
double **Diss_x_1, **Diss_y_1;
double **Diss_x_2, **Diss_y_2;
double **Diss_x_3, **Diss_y_3;
double **diver;
double **rhox, **rhoy;
double **Q_half_x_1, **Q_half_y_1;
double **Q_half_x_2, **Q_half_y_2;
double **Q_half_x_3, **Q_half_y_3;
double **E_1_compute_x_1, **E_1_compute_x_2, **E_1_compute_x_3;
double **E_2_compute_y_1, **E_2_compute_y_2, **E_2_compute_y_3;
double **U, **V;
double **DE_1_compute_x_1, **DE_1_compute_x_2, **DE_1_compute_x_3;
double **DE_2_compute_y_1, **DE_2_compute_y_2, **DE_2_compute_y_3;
double **J_half_x, **J_half_y;
double **zeta_x_half, **eta_y_half;
double **u_11, **u_12, **u_21, **u_22;
double **Ev_1_compute_x_1, **Ev_1_compute_x_2, **Ev_1_compute_x_3;
double **Ev_2_compute_y_1, **Ev_2_compute_y_2, **Ev_2_compute_y_3;
double **DEv_1_compute_x_1,**DEv_1_compute_x_2, **DEv_1_compute_x_3;
double **DEv_2_compute_y_1,**DEv_2_compute_y_2, **DEv_2_compute_y_3;
double **grad_P_x, **grad_P_y;
double **grad_dp_x, **grad_dp_y;
double **dtau_cfl, **dtau_VN;
double **dtau;
double **RHS_compute_1,**RHS_compute_2, **RHS_compute_3;
double **dp, **dp1;
double **Pressure_new;


void backward_initialize(int Nx, int Ny)
{
  initialize(&Q_1, Nx, Ny,0);
  initialize(&Q_2, Nx, Ny,0);
  initialize(&Q_3, Nx, Ny,0);
  initialize(&P_old, Nx, Ny,0);
  initialize(&U_old, Nx, Ny,0);
  initialize(&V_old, Nx, Ny,0);
  initialize(&dp, Nx, Ny,0);
  initialize(&ustar, Nx, Ny,0);
  initialize(&vstar, Nx, Ny,0);
  initialize(&u_new, Nx, Ny,0);
  initialize(&v_new, Nx, Ny,0);
  initialize(&vorticity, Nx, Ny,0);
  initialize(&rho_1, Nx, Ny,0);
  initialize(&rho_2, Nx, Ny,0);
  initialize(&Diss_x_1, Nx, Ny,0);
  initialize(&Diss_y_1, Nx, Ny,0);
  initialize(&Diss_x_2, Nx, Ny,0);
  initialize(&Diss_y_2, Nx, Ny,0);
  initialize(&Diss_x_3, Nx, Ny,0);
  initialize(&Diss_y_3, Nx, Ny,0);
  initialize(&diver, Nx, Ny,0);
  initialize(&rhox, Nx, Ny, 0);
  initialize(&rhoy, Nx, Ny, 0);
  initialize(&Q_half_x_1, Nx, Ny, 0);
  initialize(&Q_half_y_1, Nx, Ny, 0);
  initialize(&Q_half_x_2, Nx, Ny, 0);
  initialize(&Q_half_y_2, Nx, Ny, 0);
  initialize(&Q_half_x_3, Nx, Ny, 0);
  initialize(&Q_half_y_3, Nx, Ny, 0);

  initialize(&E_1_compute_x_1, Nx, Ny, 0);
  initialize(&E_1_compute_x_2, Nx, Ny, 0);
  initialize(&E_1_compute_x_3, Nx, Ny, 0);

  initialize(&E_2_compute_y_1, Nx, Ny, 0);
  initialize(&E_2_compute_y_2, Nx, Ny, 0);
  initialize(&E_2_compute_y_3, Nx, Ny, 0);
  initialize(&U, Nx, Ny, 0);
  initialize(&V, Nx, Ny, 0);
  initialize(&DE_1_compute_x_1, Nx, Ny, 0);
  initialize(&DE_1_compute_x_2, Nx, Ny, 0);
  initialize(&DE_1_compute_x_3, Nx, Ny, 0);

  initialize(&DE_2_compute_y_1, Nx, Ny, 0);
  initialize(&DE_2_compute_y_2, Nx, Ny, 0);
  initialize(&DE_2_compute_y_3, Nx, Ny, 0);

  initialize(&J_half_x, Nx, Ny, 0);
  initialize(&J_half_y, Nx, Ny, 0);

  initialize(&zeta_x_half, Nx, Ny, 0);
  initialize(&eta_y_half, Nx, Ny, 0);

  initialize(&u_11, Nx, Ny, 0);
  initialize(&u_12, Nx, Ny, 0);
  initialize(&u_21, Nx, Ny, 0);
  initialize(&u_22, Nx, Ny, 0);

  initialize(&Ev_1_compute_x_1, Nx, Ny, 0);
  initialize(&Ev_1_compute_x_2, Nx, Ny, 0);
  initialize(&Ev_1_compute_x_3, Nx, Ny, 0);

  initialize(&Ev_2_compute_y_1, Nx, Ny, 0);
  initialize(&Ev_2_compute_y_2, Nx, Ny, 0);
  initialize(&Ev_2_compute_y_3, Nx, Ny, 0);

  initialize(&DEv_1_compute_x_1, Nx, Ny, 0);
  initialize(&DEv_1_compute_x_2, Nx, Ny, 0);
  initialize(&DEv_1_compute_x_3, Nx, Ny, 0);
  
  initialize(&DEv_2_compute_y_1, Nx, Ny, 0);
  initialize(&DEv_2_compute_y_2, Nx, Ny, 0);
  initialize(&DEv_2_compute_y_3, Nx, Ny, 0);
  

  initialize(&grad_P_x, Nx, Ny, 0);
  initialize(&grad_P_y, Nx, Ny, 0);

  initialize(&dtau_cfl, Nx, Ny, 0);
  initialize(&dtau_VN, Nx, Ny, 0);

  initialize(&dtau, Nx, Ny, 0);
  
  initialize(&RHS_compute_2, Nx, Ny, 0);
  initialize(&RHS_compute_3, Nx, Ny, 0);
  initialize(&RHS_compute_1, Nx, Ny, 0);

  initialize(&dp, Nx, Ny, 0);
  initialize(&dp1, Nx, Ny, 0);

  initialize(&grad_dp_x, Nx, Ny, 0);
  initialize(&grad_dp_y, Nx, Ny, 0);

  initialize(&Pressure_new, Nx, Ny, 0);
  
}
