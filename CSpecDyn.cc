#include "CSpecDyn.h"

CSpecDyn::CSpecDyn():
N(NUM), pdims(PDIMS), dt(DT), out_dir(OUT_DIR), out_interval(OUT_INTERVAL), end_simu(END_SIMU), L(LENGTH), nu(NU), eta(ETA), setup(SETUP)
{
  
  // init MPI
  MPI_Init(NULL, NULL);
 
  // init FFT  
  FFT = MikeFFT(N, pdims);
  
  // get sizes
  FFT.get_sizes(size_R, start_R, size_F, start_F);
  
  size_R_tot = size_R[0]*size_R[1]*size_R[2];
  size_F_tot = size_F[0]*size_F[1]*size_F[2];
  
  // get MPI info
  FFT.get_comm(comm);
  
  MPI_Comm_size(comm, &nprocs);
	MPI_Comm_rank(comm, &myRank);
  MPI_Cart_coords(comm, myRank, 2, mpi_coords);
  
  // drived quantities
  XB = -L*0.5;
  dx = L/N;
  dk = PI2/L;
  time = 0.;
  
  // allocate memory
  kx = new double[size_F[0]];
  ky = new double[size_F[1]];
  kz = new double[size_F[2]];
  k2 = new double[size_F_tot];
  
  Vx_R = FFT.malloc_R();
  Vy_R = FFT.malloc_R();
  Vz_R = FFT.malloc_R();
  
  Bx_R = FFT.malloc_R();
  By_R = FFT.malloc_R();
  Bz_R = FFT.malloc_R();
  
  Vx_F = FFT.malloc_C();
  Vy_F = FFT.malloc_C();
  Vz_F = FFT.malloc_C();
  Vx_F1 = FFT.malloc_C();
  Vy_F1 = FFT.malloc_C();
  Vz_F1 = FFT.malloc_C();
  Vx_F2 = FFT.malloc_C();
  Vy_F2 = FFT.malloc_C();
  Vz_F2 = FFT.malloc_C();
  Bx_F = FFT.malloc_C();
  By_F = FFT.malloc_C();
  Bz_F = FFT.malloc_C();
  Bx_F1 = FFT.malloc_C();
  By_F1 = FFT.malloc_C();
  Bz_F1 = FFT.malloc_C();
  Bx_F2 = FFT.malloc_C();
  By_F2 = FFT.malloc_C();
  Bz_F2 = FFT.malloc_C();
  RHS_Vx_F = FFT.malloc_C();
  RHS_Vy_F = FFT.malloc_C();
  RHS_Vz_F = FFT.malloc_C();
  RHS_Vx_F1 = FFT.malloc_C();
  RHS_Vy_F1 = FFT.malloc_C();
  RHS_Vz_F1 = FFT.malloc_C();
  RHS_Vx_F2 = FFT.malloc_C();
  RHS_Vy_F2 = FFT.malloc_C();
  RHS_Vz_F2 = FFT.malloc_C();
  RHS_Bx_F = FFT.malloc_C();
  RHS_By_F = FFT.malloc_C();
  RHS_Bz_F = FFT.malloc_C();
  RHS_Bx_F1 = FFT.malloc_C();
  RHS_By_F1 = FFT.malloc_C();
  RHS_Bz_F1 = FFT.malloc_C();
  RHS_Bx_F2 = FFT.malloc_C();
  RHS_By_F2 = FFT.malloc_C();
  RHS_Bz_F2 = FFT.malloc_C();
  
  RHS_Vx_R = FFT.malloc_R();
  RHS_Vy_R = FFT.malloc_R();
  RHS_Vz_R = FFT.malloc_R();
  RHS_Bx_R = FFT.malloc_R();
  RHS_By_R = FFT.malloc_R();
  RHS_Bz_R = FFT.malloc_R();
  
  Wx_R = FFT.malloc_R();
  Wy_R = FFT.malloc_R();
  Wz_R = FFT.malloc_R();
  Jx_R = FFT.malloc_R();
  Jy_R = FFT.malloc_R();
  Jz_R = FFT.malloc_R();
  Wx_F = FFT.malloc_C();
  Wy_F = FFT.malloc_C();
  Wz_F = FFT.malloc_C();
  Jx_F = FFT.malloc_C();
  Jy_F = FFT.malloc_C();
  Jz_F = FFT.malloc_C();
  
  float_array        = (float*) malloc(sizeof(float)*size_R_tot);
  float_array_vector = (float*) malloc(sizeof(float)*size_R_tot*3);
  
  // setup initial fields
  setup_k();
  setup_fields();
  
  // create output directory
  if(myRank == 0)
	{
    mkdir(out_dir.c_str(), 0777);
	}MPI_Barrier(comm);
  
  // create derived MPI_datatypes for vti output
  int size_total[3] = {N,N,N};
  MPI_Type_create_subarray(3, size_total, size_R, start_R, MPI_ORDER_C, MPI_FLOAT, &vti_subarray);
  MPI_Type_commit(&vti_subarray);
  
  MPI_Type_contiguous(3, MPI_FLOAT, &vti_float3);
  MPI_Type_commit(&vti_float3);
  
  MPI_Type_create_subarray(3, size_total, size_R, start_R, MPI_ORDER_C, vti_float3, &vti_subarray_vector);
  MPI_Type_commit(&vti_subarray_vector);

  // random
  normal_eng = std::mt19937(42);
  normal = std::normal_distribution<double>(0.,1.);
  
  // Energy Spectrum
  energySpectrum_V = new double[N_bin];
  energySpectrum_V_loc = new double[N_bin];
  bin_counter_V    = new int[N_bin];
  bin_counter_V_loc    = new int[N_bin];
  energySpectrum_B = new double[N_bin];
  energySpectrum_B_loc = new double[N_bin];
  bin_counter_B    = new int[N_bin];
  bin_counter_B_loc    = new int[N_bin];
  
}

void CSpecDyn::setup_k()
{
  // kx
  for(int ix = 0; ix < size_F[0]; ix++)
	{
    
    int ix_glob = start_F[0] + ix;
    
    if(ix_glob < N/2)
    {
      kx[ix] = ix_glob * dk;
    }
    else
    {
      kx[ix] = (ix_glob-N) * dk;
		}
	}
  
  // ky
  for(int iy = 0; iy < size_F[1]; iy++)
	{
    
    int iy_glob = start_F[1] + iy;
    
    if(iy_glob < N/2)
    {
      ky[iy] = iy_glob * dk;
    }
    else
    {
      ky[iy] = (iy_glob-N) * dk;
		}
	}
  
  // kx
  for(int iz = 0; iz < size_F[2]; iz++)
	{
		int iz_glob = start_F[2] + iz;
		
    kz[iz] = iz_glob * dk;
	}
  
  // k2
  for(int ix = 0; ix<size_F[0]; ix++){
  for(int iy = 0; iy<size_F[1]; iy++){
  for(int iz = 0; iz<size_F[2]; iz++){

    int id = ix * size_F[1]*size_F[2] + iy * size_F[2] + iz;
    
    k2[id] = kx[ix]*kx[ix] + ky[iy]*ky[iy] + kz[iz]*kz[iz];
    
    k2[id] = std::max(k2[id], 1.e-12); 
    
  }}}
  
}

void CSpecDyn::setup_fields()
{
  std::mt19937 eng(myRank);
  std::uniform_real_distribution<double> phi(0.,PI2);
  double norm_loc = 0.;
  double norm;
  double s = 11./6.;
  
  switch(setup)
  {
    case 0:
      /** u = 0, B = 0 **/
      for(int id = 0; id < size_R_tot; id++)
      {    
        Vx_R[id] = 0.;
        Vy_R[id] = 0.;
        Vz_R[id] = 0.;
        
        Bx_R[id] = 0.;
        By_R[id] = 0.;
        Bz_R[id] = 0.;
      }
      fFFT(Vx_R, Vy_R, Vz_R, Vx_F, Vy_F, Vz_F);
      fFFT(Bx_R, By_R, Bz_R, Bx_F, By_F, Bz_F);
      break;
      
    case 1:
      /** Orszag-Tang **/
      for(int ix = 0; ix < size_R[0]; ix++){
      for(int iy = 0; iy < size_R[1]; iy++){
      for(int iz = 0; iz < size_R[2]; iz++){
      
        int id = ix * size_R[1]*size_R[2] + iy * size_R[2] + iz;
        
        double x_val = (start_R[0]+ix)*dx+XB;
        double y_val = (start_R[1]+iy)*dx+XB;
        double z_val = (start_R[2]+iz)*dx+XB;
        
        Vx_R[id] = -2.*sin(y_val);
        Vy_R[id] =  2.*sin(x_val);
        Vz_R[id] =  0.;
        
        double beta = 0.8;
        
        Bx_R[id] = ( -2.*sin(y_val) + sin(z_val) )*beta;
        By_R[id] = (  2.*sin(x_val) + sin(z_val) )*beta;
        Bz_R[id] = (     sin(x_val) + sin(y_val) )*beta;
    
      }}}
      fFFT(Vx_R, Vy_R, Vz_R, Vx_F, Vy_F, Vz_F);
      fFFT(Bx_R, By_R, Bz_R, Bx_F, By_F, Bz_F);
      break;
    
    case 2:
      /** random with energy spectrum **/
      
      // get normalization
      for(int id = 0; id < size_F_tot; id++){
        if(k2[id] > 0.1*dk*dk)
        {
          norm_loc += 1./pow(1+k2[id],s);
        }
      }
      MPI_Allreduce(&norm_loc, &norm, 1, MPI_DOUBLE, MPI_SUM, comm);
      
      norm = 1./norm;
      
      // produce energy spectrum
      for(int ix = 0; ix < size_F[0]; ix++){
      for(int iy = 0; iy < size_F[1]; iy++){
      for(int iz = 0; iz < size_F[2]; iz++){
      
        int id = ix * size_F[1]*size_F[2] + iy * size_F[2] + iz;
        
        double A = sqrt( norm*1./pow(1+k2[id],s) );
        double k2_inv = 1./k2[id];
        double k_x = kx[ix];
        double k_y = ky[iy];
        double k_z = kz[iz];
        
        Vx_F[id] = A * ( + ( 1. - k_x*k_x*k2_inv )*exp(IM*phi(eng)) -        k_x*k_y*k2_inv  *exp(IM*phi(eng)) -        k_x*k_z*k2_inv  *exp(IM*phi(eng)) );
        Vy_F[id] = A * ( -        k_y*k_x*k2_inv  *exp(IM*phi(eng)) + ( 1. - k_y*k_y*k2_inv )*exp(IM*phi(eng)) -        k_y*k_z*k2_inv  *exp(IM*phi(eng)) );
        Vz_F[id] = A * ( -        k_z*k_x*k2_inv  *exp(IM*phi(eng)) -        k_z*k_y*k2_inv  *exp(IM*phi(eng)) + ( 1. - k_z*k_z*k2_inv )*exp(IM*phi(eng)) );
        
        Bx_F[id] = A * ( + ( 1. - k_x*k_x*k2_inv )*exp(IM*phi(eng)) -        k_x*k_y*k2_inv  *exp(IM*phi(eng)) -        k_x*k_z*k2_inv  *exp(IM*phi(eng)) );
        By_F[id] = A * ( -        k_y*k_x*k2_inv  *exp(IM*phi(eng)) + ( 1. - k_y*k_y*k2_inv )*exp(IM*phi(eng)) -        k_y*k_z*k2_inv  *exp(IM*phi(eng)) );
        Bz_F[id] = A * ( -        k_z*k_x*k2_inv  *exp(IM*phi(eng)) -        k_z*k_y*k2_inv  *exp(IM*phi(eng)) + ( 1. - k_z*k_z*k2_inv )*exp(IM*phi(eng)) );
        
      }}}
      
      // normalize to |V|_max=1.
      bFFT(Vx_F, Vy_F, Vz_F, Vx_R, Vy_R, Vz_R);
      bFFT(Bx_F, By_F, Bz_F, Bx_R, By_R, Bz_R);
      norm_loc = 0.;
      for(int id = 0; id < size_R_tot; id++){
        norm_loc = std::max(sqrt(Vx_R[id]*Vx_R[id]+Vy_R[id]*Vy_R[id]+Vz_R[id]*Vz_R[id]), norm_loc);
      }
      MPI_Allreduce(&norm_loc, &norm, 1, MPI_DOUBLE, MPI_MAX, comm);
      for(int id = 0; id < size_R_tot; id++){
        Vx_R[id] /= norm;
        Vy_R[id] /= norm;
        Vz_R[id] /= norm;
        Bx_R[id] /= norm;
        By_R[id] /= norm;
        Bz_R[id] /= norm;
      }
      fFFT(Vx_R, Vy_R, Vz_R, Vx_F, Vy_F, Vz_F);
      fFFT(Bx_R, By_R, Bz_R, Bx_F, By_F, Bz_F);
      
      break;
    
    default: 
      if(myRank==0){printf("No valid setup provided! setup = %d\n", setup);}
      MPI_Barrier(comm);
      MPI_Finalize();
      exit(EXIT_FAILURE);
  }
  
}

void CSpecDyn::execute()
{
  double start_time = MPI_Wtime();
  double out_time = time;
  
  calc_EnergySpectrum();
  print_vti();

  while(time < end_simu)
  {
    time_step();
    out_time += dt;
    
    if(out_time > out_interval)
    {
      calc_EnergySpectrum();
      print_vti();
      out_time -= out_interval;
    }
    
  }
  
  double print_time = MPI_Wtime() - start_time;
  if(myRank==0){printf("Print time = %f, pdims = [%d,%d], N = %d\n", print_time, pdims[0], pdims[1], N);}
}

void CSpecDyn::time_step()
{
  OrnsteinUhlenbeck();
  
  double del_t;
  
  // step 1
  del_t = 1.;
  
  calc_RHS(RHS_Vx_F , RHS_Vy_F , RHS_Vz_F , Vx_F , Vy_F , Vz_F
          ,RHS_Bx_F , RHS_By_F , RHS_Bz_F , Bx_F , By_F , Bz_F, del_t, 0);
  
  for(int id = 0; id < size_F_tot; id++)
  { 
    Vx_F1[id] = Vx_F[id] + dt * RHS_Vx_F[id];
    Vy_F1[id] = Vy_F[id] + dt * RHS_Vy_F[id];
    Vz_F1[id] = Vz_F[id] + dt * RHS_Vz_F[id];
    
    Bx_F1[id] = Bx_F[id] + dt * RHS_Bx_F[id];
    By_F1[id] = By_F[id] + dt * RHS_By_F[id];
    Bz_F1[id] = Bz_F[id] + dt * RHS_Bz_F[id];
  }
  
  diffusion_correction(Vx_F1, Vy_F1, Vz_F1, Bx_F1, By_F1, Bz_F1, del_t);
  
  projection(Vx_F1, Vy_F1, Vz_F1);
  projection(Bx_F1, By_F1, Bz_F1);
          
  // step 2
  del_t = 0.5;
  
  calc_RHS(RHS_Vx_F1, RHS_Vy_F1, RHS_Vz_F1, Vx_F1, Vy_F1, Vz_F1
          ,RHS_Bx_F1, RHS_By_F1, RHS_Bz_F1, Bx_F1, By_F1, Bz_F1, del_t, 1);
   
  double dt_025 = 0.25*dt; 
          
  for(int id = 0; id < size_F_tot; id++)
  { 
    Vx_F2[id] = Vx_F[id] + dt_025 * ( RHS_Vx_F[id] + RHS_Vx_F1[id]);
    Vy_F2[id] = Vy_F[id] + dt_025 * ( RHS_Vy_F[id] + RHS_Vy_F1[id]);
    Vz_F2[id] = Vz_F[id] + dt_025 * ( RHS_Vz_F[id] + RHS_Vz_F1[id]);
    
    Bx_F2[id] = Bx_F[id] + dt_025 * ( RHS_Bx_F[id] + RHS_Bx_F1[id]);
    By_F2[id] = By_F[id] + dt_025 * ( RHS_By_F[id] + RHS_By_F1[id]);
    Bz_F2[id] = Bz_F[id] + dt_025 * ( RHS_Bz_F[id] + RHS_Bz_F1[id]);
  }
  
  diffusion_correction(Vx_F2, Vy_F2, Vz_F2, Bx_F2, By_F2, Bz_F2, del_t);
  
  projection(Vx_F2, Vy_F2, Vz_F2);
  projection(Bx_F2, By_F2, Bz_F2);
          
  // step 3
  del_t = 1.;
  
  calc_RHS(RHS_Vx_F2, RHS_Vy_F2, RHS_Vz_F2, Vx_F2, Vy_F2, Vz_F2
          ,RHS_Bx_F2, RHS_By_F2, RHS_Bz_F2, Bx_F2, By_F2, Bz_F2, del_t, 2);
  
  double dt_6 = dt/6.;
  
  for(int id = 0; id < size_F_tot; id++)
  { 
    Vx_F[id] = Vx_F[id] + dt_6 * (RHS_Vx_F[id] + RHS_Vx_F1[id] + 4.*RHS_Vx_F2[id]);
    Vy_F[id] = Vy_F[id] + dt_6 * (RHS_Vy_F[id] + RHS_Vy_F1[id] + 4.*RHS_Vy_F2[id]);
    Vz_F[id] = Vz_F[id] + dt_6 * (RHS_Vz_F[id] + RHS_Vz_F1[id] + 4.*RHS_Vz_F2[id]);
    
    Bx_F[id] = Bx_F[id] + dt_6 * (RHS_Bx_F[id] + RHS_Bx_F1[id] + 4.*RHS_Bx_F2[id]);
    By_F[id] = By_F[id] + dt_6 * (RHS_By_F[id] + RHS_By_F1[id] + 4.*RHS_By_F2[id]);
    Bz_F[id] = Bz_F[id] + dt_6 * (RHS_Bz_F[id] + RHS_Bz_F1[id] + 4.*RHS_Bz_F2[id]);

  }
  
  diffusion_correction(Vx_F, Vy_F, Vz_F, Bx_F, By_F, Bz_F, del_t);
  
  projection(Vx_F , Vy_F , Vz_F );
  projection(Bx_F , By_F , Bz_F );
  
  // update time
  time += dt;
  
  if(myRank == 0)
	{
		printf("  time step: time = %f, dt = %f\n", time, dt);
	}MPI_Barrier(comm);
}

void CSpecDyn::calc_RHS(CX* RHSV_X, CX* RHSV_Y, CX* RHSV_Z, CX* V_X, CX* V_Y, CX* V_Z,
                        CX* RHSB_X, CX* RHSB_Y, CX* RHSB_Z, CX* B_X, CX* B_Y, CX* B_Z,
                        double del_t, int substep)
{
  // W = rot(V)
  for(int ix = 0; ix<size_F[0]; ix++){
  for(int iy = 0; iy<size_F[1]; iy++){
  for(int iz = 0; iz<size_F[2]; iz++){
    
    int id = ix * size_F[1]*size_F[2] + iy * size_F[2] + iz;
    
    Wx_F[id] = IM * ( ky[iy]*V_Z[id] - kz[iz]*V_Y[id] );
    Wy_F[id] = IM * ( kz[iz]*V_X[id] - kx[ix]*V_Z[id] );
    Wz_F[id] = IM * ( kx[ix]*V_Y[id] - ky[iy]*V_X[id] );
    
  }}}
  
  // J = rot(B)
  for(int ix = 0; ix<size_F[0]; ix++){
  for(int iy = 0; iy<size_F[1]; iy++){
  for(int iz = 0; iz<size_F[2]; iz++){
    
    int id = ix * size_F[1]*size_F[2] + iy * size_F[2] + iz;
    
    Jx_F[id] = IM * ( ky[iy]*B_Z[id] - kz[iz]*B_Y[id] );
    Jy_F[id] = IM * ( kz[iz]*B_X[id] - kx[ix]*B_Z[id] );
    Jz_F[id] = IM * ( kx[ix]*B_Y[id] - ky[iy]*B_X[id] );
    
  }}}
  
  // dealias
  dealias(V_X, V_Y, V_Z);
  dealias(B_X, B_Y, B_Z);
  dealias(Wx_F, Wy_F, Wz_F);
  dealias(Jx_F, Jy_F, Jz_F);
  
  // FFT F->R
  bFFT(V_X, V_Y, V_Z, Vx_R, Vy_R, Vz_R);
  bFFT(B_X, B_Y, B_Z, Bx_R, By_R, Bz_R);
  bFFT(Wx_F, Wy_F, Wz_F, Wx_R, Wy_R, Wz_R);  
  bFFT(Jx_F, Jy_F, Jz_F, Jx_R, Jy_R, Jz_R);
  
  // RHS_V = VxW + JxB
  for(int id = 0; id < size_R_tot; id++){
    
    RHS_Vx_R[id] = Vy_R[id]*Wz_R[id]-Vz_R[id]*Wy_R[id] + Jy_R[id]*Bz_R[id]-Jz_R[id]*By_R[id];
    RHS_Vy_R[id] = Vz_R[id]*Wx_R[id]-Vx_R[id]*Wz_R[id] + Jz_R[id]*Bx_R[id]-Jx_R[id]*Bz_R[id];
    RHS_Vz_R[id] = Vx_R[id]*Wy_R[id]-Vy_R[id]*Wx_R[id] + Jx_R[id]*By_R[id]-Jy_R[id]*Bx_R[id];
    
  }
  
  // Taylor-Green forcing
  if(setup==2)
  {
    for(int ix = 0; ix < size_R[0]; ix++){
    for(int iy = 0; iy < size_R[1]; iy++){
    for(int iz = 0; iz < size_R[2]; iz++){
    
      int id = ix * size_R[1]*size_R[2] + iy * size_R[2] + iz;
      
      double x_val = (start_R[0]+ix)*dx+XB;
      double y_val = (start_R[1]+iy)*dx+XB;
      double z_val = (start_R[2]+iz)*dx+XB;
      
      double kf = 2.;
      
      RHS_Vx_R[id] += 0.125 * ( sin(kf*x_val)*cos(kf*y_val)*cos(kf*z_val) );
      RHS_Vy_R[id] -= 0.125 * ( cos(kf*x_val)*sin(kf*y_val)*cos(kf*z_val) );
      
    }}}
  }
  
  // RHS_B = VxB
  for(int id = 0; id < size_R_tot; id++){
    
    RHS_Bx_R[id] = Vy_R[id]*Bz_R[id]-Vz_R[id]*By_R[id];
    RHS_By_R[id] = Vz_R[id]*Bx_R[id]-Vx_R[id]*Bz_R[id];
    RHS_Bz_R[id] = Vx_R[id]*By_R[id]-Vy_R[id]*Bx_R[id];
    
  }
  
  // FFT R->F
  fFFT(Vx_R, Vy_R, Vz_R, V_X, V_Y, V_Z);
  fFFT(Bx_R, By_R, Bz_R, B_X, B_Y, B_Z);
  fFFT(RHS_Vx_R, RHS_Vy_R, RHS_Vz_R, RHSV_X, RHSV_Y, RHSV_Z);
  fFFT(RHS_Bx_R, RHS_By_R, RHS_Bz_R, RHSB_X, RHSB_Y, RHSB_Z);
  
  // dealias
  dealias(V_X, V_Y, V_Z);
  dealias(B_X, B_Y, B_Z);
  dealias(RHSV_X, RHSV_Y, RHSV_Z);
  dealias(RHSB_X, RHSB_Y, RHSB_Z);
  
  // Ornstein-Uhlenbeck forcing
  //~ double dk2 = dk*dk;
  
  //~ for(int ix = 0; ix<size_F[0]; ix++){
  //~ for(int iy = 0; iy<size_F[1]; iy++){
  //~ for(int iz = 0; iz<size_F[2]; iz++){
    
    //~ int id = ix * size_F[1]*size_F[2] + iy * size_F[2] + iz;
    
    //~ if(0.1*dk2 < k2[id] && k2[id] < 9.1*dk2) // force first few modes 
    //~ {
      
      //~ double k2_inv = 1./k2[id];
      //~ double k_x = kx[ix];
      //~ double k_y = ky[iy];
      //~ double k_z = kz[iz];

      //~ RHSV_X[id] += f_OU_X[substep] * (+ ( 1. - k_x*k_x*k2_inv ) -        k_x*k_y*k2_inv   -        k_x*k_z*k2_inv  );
      //~ RHSV_Y[id] += f_OU_Y[substep] * (-        k_y*k_x*k2_inv   + ( 1. - k_y*k_y*k2_inv ) -        k_y*k_z*k2_inv  );
      //~ RHSV_Z[id] += f_OU_Z[substep] * (-        k_z*k_x*k2_inv   -        k_z*k_y*k2_inv   + ( 1. - k_z*k_z*k2_inv ));

    //~ }
    
  //~ }}}
  
  // RHS_B = rot(VxB)
  for(int ix = 0; ix<size_F[0]; ix++){
  for(int iy = 0; iy<size_F[1]; iy++){
  for(int iz = 0; iz<size_F[2]; iz++){
    
    int id = ix * size_F[1]*size_F[2] + iy * size_F[2] + iz;
    
    CX VxB_X = RHSB_X[id];
    CX VxB_Y = RHSB_Y[id];
    CX VxB_Z = RHSB_Z[id];
    
    RHSB_X[id] = IM * ( ky[iy]*VxB_Z - kz[iz]*VxB_Y );
    RHSB_Y[id] = IM * ( kz[iz]*VxB_X - kx[ix]*VxB_Z );
    RHSB_Z[id] = IM * ( kx[ix]*VxB_Y - ky[iy]*VxB_X );
    
  }}}
  
  // diffusion
  for(int id = 0; id < size_F_tot; id++)
  {
    double exp_diff_V = exp(nu *k2[id]*del_t*dt);
    double exp_diff_B = exp(eta*k2[id]*del_t*dt);
    
    RHSV_X[id] *= exp_diff_V;
    RHSV_X[id] *= exp_diff_V;
    RHSV_X[id] *= exp_diff_V;
    
    RHSB_X[id] *= exp_diff_B;
    RHSB_X[id] *= exp_diff_B;
    RHSB_X[id] *= exp_diff_B;
  }
}

void CSpecDyn::diffusion_correction(CX* Vx, CX* Vy, CX* Vz, CX* Bx, CX* By, CX* Bz, double del_t)
{

  double exp_diff_V;
  double exp_diff_B;

  for(int id = 0; id < size_F_tot; id++)
  {
    exp_diff_V = exp(- nu *k2[id]*del_t*dt);
    exp_diff_B = exp(- eta*k2[id]*del_t*dt);
    
    Vx[id] *= exp_diff_V;
    Vy[id] *= exp_diff_V;
    Vz[id] *= exp_diff_V;
    
    Bx[id] *= exp_diff_B;
    By[id] *= exp_diff_B;
    Bz[id] *= exp_diff_B;
  }
}

void CSpecDyn::finalize()
{
  MPI_Finalize();
}

void CSpecDyn::projection(CX* fieldX, CX* fieldY, CX* fieldZ)
{

  for(int ix = 0; ix<size_F[0]; ix++){
  for(int iy = 0; iy<size_F[1]; iy++){
  for(int iz = 0; iz<size_F[2]; iz++){
    
    int id = ix * size_F[1]*size_F[2] + iy * size_F[2] + iz;
    
    CX temp = (kx[ix] * fieldX[id] + ky[iy] * fieldY[id] + kz[iz] * fieldZ[id]) / k2[id];
    
    fieldX[id] -= kx[ix] * temp;
    fieldY[id] -= ky[iy] * temp;
    fieldZ[id] -= kz[iz] * temp;	
    
  }}}
  
}

void CSpecDyn::fFFT(double* IN_x, double* IN_y, double* IN_z, CX* OUT_x, CX* OUT_y, CX* OUT_z)
{
  
  FFT.r2c(IN_x, OUT_x);
  FFT.r2c(IN_y, OUT_y);
  FFT.r2c(IN_z, OUT_z);
  
}

void CSpecDyn::bFFT(CX* IN_x, CX* IN_y, CX* IN_z, double* OUT_x, double* OUT_y, double* OUT_z)
{
  
  FFT.c2r(IN_x, OUT_x);
  FFT.c2r(IN_y, OUT_y);
  FFT.c2r(IN_z, OUT_z);
  
}

void CSpecDyn::print_vti()
{
  bFFT(Vx_F, Vy_F, Vz_F, Vx_R, Vy_R, Vz_R);
  bFFT(Bx_F, By_F, Bz_F, Bx_R, By_R, Bz_R);
  
  std::string file_name  = out_dir + "/step_" + std::to_string(vti_count) + ".vti";
  std::ofstream os;
  
  int offset = 0;
	int N_tot = N*N*N;
	int N_bytes_scalar  =   N_tot * sizeof(float);
	int N_bytes_vector  = 3*N_tot * sizeof(float);
  int bin_size_scalar = N_bytes_scalar + sizeof(uint64_t);// 2nd term is the size of the the leading integer announcing the numbers n the data chunk
  int bin_size_vector = N_bytes_vector + sizeof(uint64_t);
  
  // header
  if(myRank==0)
  {
    printf("Printing vti # %d!\n", vti_count);
    
    os.open(file_name.c_str(), std::ios::out);
    if(!os){
      std::cout << "Cannot write header to file '" << file_name << "'!\n";
    }
    
    // write header	
		int extend_l[3]  = {0, 0, 0};
		int extend_r[3]  = {N-1, N-1, N-1};
		double origin[3] = {XB,XB,XB};
    
    os << "<VTKFile type=\"ImageData\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt64\">" << std::endl;	
    os << "  <ImageData WholeExtent=\"" << extend_l[0] << " " << extend_r[0] << " " 
                                        << extend_l[1] << " " << extend_r[1] << " " 
                                        << extend_l[2] << " " << extend_r[2] 
				 << "\" Origin=\""  << origin[0]  << " " << origin[1]  << " " << origin[2] 
				 << "\" Spacing=\"" << dx << " " << dx << " " << dx
         << "\" Direction=\"0 0 1 0 1 0 1 0 0\">" << std::endl; // C -> FORTRAN order
    
    os << "      <FieldData>" << std::endl;
    os << "        <DataArray type=\"Float32\" Name=\"TimeValue\" NumberOfTuples=\"1\" format=\"ascii\">" << std::endl;
    os << "        "<< float(time) << std::endl;
    os << "        </DataArray>" << std::endl;
    os << "      </FieldData>" << std::endl;
        
		os << "    <Piece Extent=\"" << extend_l[0] << " " << extend_r[0] << " " 
                                 << extend_l[1] << " " << extend_r[1] << " " 
                                 << extend_l[2] << " " << extend_r[2] << "\">" << std::endl;
    
    os << "      <PointData Scalars=\"P\" Vectors=\"V\">" << std::endl;
    
    os << "        <DataArray type=\"Float32\" Name=\"V\" NumberOfComponents=\"3\" format=\"appended\" offset=\"" << offset << "\">" << std::endl;
    os << "        </DataArray>" << std::endl;
    offset += bin_size_vector;
    os << "        <DataArray type=\"Float32\" Name=\"B\" NumberOfComponents=\"3\" format=\"appended\" offset=\"" << offset << "\">" << std::endl;
    os << "        </DataArray>" << std::endl;
    offset += bin_size_vector;
    
    os << "      </PointData>" << std::endl;
    os << "      <CellData>" << std::endl;
    os << "      </CellData>" << std::endl;
    os << "    </Piece>" << std::endl;
    os << "  </ImageData>" << std::endl;
    os << "  <AppendedData encoding=\"raw\">" << std::endl;
    os << "   _";
                                
    os.close();
  
  }MPI_Barrier(comm);
  
  // print fields in parallel  
  print_mpi_vector(Vx_R, Vy_R, Vz_R, N_bytes_vector, file_name.c_str()); 
  print_mpi_vector(Bx_R, By_R, Bz_R, N_bytes_vector, file_name.c_str()); 
  
  // footer
  if(myRank==0)
  {
    os.open(file_name.c_str(), std::ios::out | std::ios::app);
    if(!os){
      std::cout << "Cannot write footer to file '" << file_name << "'!\n";
      exit(3);
    }
		
		os << std::endl << "  </AppendedData>" << std::endl;
    os<< "</VTKFile>" << std::endl;
	
    os.close();
  }MPI_Barrier(comm);
  
  vti_count++;
  
  fFFT(Vx_R, Vy_R, Vz_R, Vx_F, Vy_F, Vz_F);
  fFFT(Bx_R, By_R, Bz_R, Bx_F, By_F, Bz_F);
}

void CSpecDyn::print_mpi_vector(double* field_X, double* field_Y, double* field_Z, int& N_bytes_vector, const char* file_name)
{
  if(myRank==0)
  {
    std::ofstream binary_os(file_name, std::ios::out | std::ios::app | std::ios::binary );
    binary_os.write(reinterpret_cast<const char*>(&N_bytes_vector),sizeof(uint64_t)); // size of following binary package
    binary_os.close();
  }MPI_Barrier(comm);
  
  // open file
  MPI_File mpi_file;
  MPI_File_open(comm, file_name, MPI_MODE_APPEND|MPI_MODE_WRONLY, MPI_INFO_NULL, &mpi_file);
  
  // offset to end of file
  MPI_Offset mpi_eof;
  MPI_File_get_position(mpi_file, &mpi_eof);
  MPI_Barrier(comm);
  
  // data to float array
  for(int id = 0; id < size_R_tot; id++)
  {
    float_array_vector[3*id+0] = float(field_X[id]);
    float_array_vector[3*id+1] = float(field_Y[id]);
    float_array_vector[3*id+2] = float(field_Z[id]);
  }
  
  // write data
  MPI_File_set_view(mpi_file, mpi_eof, vti_float3, vti_subarray_vector, "native", MPI_INFO_NULL);
  MPI_File_write_all(mpi_file, float_array_vector, size_R_tot, vti_float3, MPI_STATUS_IGNORE);
  
  // close file
  MPI_File_close(&mpi_file);  
}

void CSpecDyn::print_mpi_scalar(double* field, int& N_bytes_scalar, const char* file_name)
{
  if(myRank==0)
  {
    std::ofstream binary_os(file_name, std::ios::out | std::ios::app | std::ios::binary );
    binary_os.write(reinterpret_cast<const char*>(&N_bytes_scalar),sizeof(uint64_t)); // size of following binary package
    binary_os.close();
  }MPI_Barrier(comm);
  
  // open file
  MPI_File mpi_file;
  MPI_File_open(comm, file_name, MPI_MODE_APPEND|MPI_MODE_WRONLY, MPI_INFO_NULL, &mpi_file);
  
  // offset to end of file
  MPI_Offset mpi_eof;
  MPI_File_get_position(mpi_file, &mpi_eof);
  MPI_Barrier(comm);
  
  // data to float array
  for(int id = 0; id < size_R_tot; id++)
  {
    float_array[id] = float(field[id]);
  }
  
  // write data
  MPI_File_set_view(mpi_file, mpi_eof, MPI_FLOAT, vti_subarray, "native", MPI_INFO_NULL);
  MPI_File_write_all(mpi_file, float_array, size_R_tot, MPI_FLOAT, MPI_STATUS_IGNORE);
  
  // close file
  MPI_File_close(&mpi_file);  
}

void CSpecDyn::dealias(CX* fieldX, CX* fieldY, CX* fieldZ)
{
  // spherical truncation
  double kmax  = sqrt(2)/3.*N/2.*dk;
  double kmax2 = kmax*kmax;
  
  //~ for(int id = 0; id < size_F_tot; id++){
   
    //~ if(k2[id] > kmax)
    //~ {
      //~ fieldX[id] = 0.;
      //~ fieldY[id] = 0.;
      //~ fieldZ[id] = 0.;
    //~ }
    
  //~ }
  
  // 2/3 rule
  //~ double kmax = N/2.*dk;
  //~ double kmax_23 = N/2.*dk*2./3.;
  
  //~ for(int ix = 0; ix<size_F[0]; ix++){
  //~ for(int iy = 0; iy<size_F[1]; iy++){
  //~ for(int iz = 0; iz<size_F[2]; iz++){
    
    //~ int id = ix * size_F[1]*size_F[2] + iy * size_F[2] + iz;
    
    //~ if( fabs(kx[ix]) > kmax_23 ||  fabs(ky[iy]) > kmax_23 ||  fabs(kz[iz]) > kmax_23) // Cube
    //~ {
      //~ fieldX[id] = 0.;
      //~ fieldY[id] = 0.;
      //~ fieldZ[id] = 0.;
    //~ }
    
  //~ }}}
}

void CSpecDyn::OrnsteinUhlenbeck()
{
  double kf    = 2.;
  double T     = 0.1;
  double T_inv = 1./T;
  double sigma = 0.25*sqrt(2.*T_inv);
 
  CX rand[2];
  
  // X
  rand[0] = normal(normal_eng) + IM * normal(normal_eng);
  rand[1] = normal(normal_eng) + IM * normal(normal_eng);
  
  f_OU_X[0] = f_OU_X[2];
  f_OU_X[1] = f_OU_X[0] - 0.5*dt*T_inv*f_OU_X[0] + sqrt(0.5*dt)*sigma*rand[0];
  f_OU_X[2] = f_OU_X[1] - 0.5*dt*T_inv*f_OU_X[1] + sqrt(0.5*dt)*sigma*rand[1];
  
  // Y
  rand[0] = normal(normal_eng) + IM * normal(normal_eng);
  rand[1] = normal(normal_eng) + IM * normal(normal_eng);
  
  f_OU_Y[0] = f_OU_Y[2];
  f_OU_Y[1] = f_OU_Y[0] - 0.5*dt*T_inv*f_OU_Y[0] + sqrt(0.5*dt)*sigma*rand[0];
  f_OU_Y[2] = f_OU_Y[1] - 0.5*dt*T_inv*f_OU_Y[1] + sqrt(0.5*dt)*sigma*rand[1];
  
  // Z
  rand[0] = normal(normal_eng) + IM * normal(normal_eng);
  rand[1] = normal(normal_eng) + IM * normal(normal_eng);
  
  f_OU_Z[0] = f_OU_Z[2];
  f_OU_Z[1] = f_OU_Z[0] - 0.5*dt*T_inv*f_OU_Z[0] + sqrt(0.5*dt)*sigma*rand[0];
  f_OU_Z[2] = f_OU_Z[1] - 0.5*dt*T_inv*f_OU_Z[1] + sqrt(0.5*dt)*sigma*rand[1];
  
}

// inspired by https://bertvandenbroucke.netlify.app/2019/05/24/computing-a-power-spectrum-in-python/
//~ void CSpecDyn::calc_EnergySpectrum()
//~ {
  //~ int kmax = N/2+1;
  
  //~ // clean containers
  //~ for(int ik = 0; ik < kmax; ik++)
  //~ {
    //~ energySpectrum_V_loc[ik] = 0.;
    //~ bin_counter_V_loc[ik]    = 0;
    //~ energySpectrum_B_loc[ik] = 0.;
    //~ bin_counter_B_loc[ik]    = 0;
  //~ }
  
  //~ double Vx, Vy, Vz;
  //~ double Bx, By, Bz;
  
  //~ for(int id = 0; id < size_F_tot; id++){
    
    //~ int id_k = int(round(sqrt(k2[id])/dk));
  
    //~ if(id_k < kmax)
    //~ {
      //~ Vx = abs(Vx_F[id]);
      //~ Vy = abs(Vy_F[id]);
      //~ Vz = abs(Vz_F[id]);
      //~ Bx = abs(Bx_F[id]);
      //~ By = abs(By_F[id]);
      //~ Bz = abs(Bz_F[id]);
      
      //~ energySpectrum_V_loc[id_k] += Vx*Vx+Vy*Vy+Vz*Vz;
      //~ energySpectrum_B_loc[id_k] += Bx*Bx+By*By+Bz*Bz;
      //~ bin_counter_V_loc   [id_k] += 1;
      //~ bin_counter_B_loc   [id_k] += 1;
    //~ }
  //~ }
  
  //~ MPI_Reduce(energySpectrum_V_loc, energySpectrum_V, kmax, MPI_DOUBLE, MPI_SUM, 0, comm);
  //~ MPI_Reduce(bin_counter_V_loc   , bin_counter_V   , kmax, MPI_INT   , MPI_SUM, 0, comm);
  //~ MPI_Reduce(energySpectrum_B_loc, energySpectrum_B, kmax, MPI_DOUBLE, MPI_SUM, 0, comm);
  //~ MPI_Reduce(bin_counter_B_loc   , bin_counter_B   , kmax, MPI_INT   , MPI_SUM, 0, comm);
  
  //~ if(myRank == 0)
  //~ {
    //~ std::string file_name  = out_dir + "/spectrum_" + std::to_string(vti_count) + ".csv";
    //~ std::ofstream os;
    
    //~ os.open(file_name.c_str(), std::ios::out);
    //~ if(!os){
      //~ std::cout << "Cannot write header to file '" << file_name << "'!\n";
    //~ }
    
    //~ for(int ik = 1; ik < kmax; ik++)
    //~ {
      //~ energySpectrum_V[ik] /= double(bin_counter_V[ik]); // divide by number of elements in bin
      //~ energySpectrum_V[ik] *= 0.5*M_PI*dk*dk* ( (ik+0.5)*(ik+0.5) - (ik-0.5)*(ik-0.5) ); // multiply with Fourier surface elemets
      //~ energySpectrum_B[ik] /= double(bin_counter_B[ik]);
      //~ energySpectrum_B[ik] *= 0.5*M_PI*dk*dk* ( (ik+0.5)*(ik+0.5) - (ik-0.5)*(ik-0.5) );
      //~ os << ik*dk << ", " << energySpectrum_V[ik] << ", " << energySpectrum_B[ik] << std::endl;
    //~ }
    
    //~ double E_v = 0.;
    //~ for(int ik = 1; ik < kmax; ik++)
    //~ {
      //~ E_v += energySpectrum_V[ik];
    //~ }
    //~ E_v *= dk;
    
    //~ printf("Total Energy = %f\n", E_v);
    
    //~ os.close();
  //~ }MPI_Barrier(comm);
  
//~ }

void CSpecDyn::calc_EnergySpectrum()
{
  double kmax = N/2*dk;
  double del_k = kmax/N_bin;
  
  // clean containers
  for(int ik = 0; ik < N_bin; ik++)
  {
    energySpectrum_V_loc[ik] = 0.;
    bin_counter_V_loc[ik]    = 0;
    energySpectrum_B_loc[ik] = 0.;
    bin_counter_B_loc[ik]    = 0;
  }
  
  double Vx, Vy, Vz;
  double Bx, By, Bz;
  
  for(int id = 0; id < size_F_tot; id++){
    
    double k = sqrt(k2[id]);
    int id_k = int(k/del_k);
  
    if(id_k < N_bin)
    {
      Vx = abs(Vx_F[id]);
      Vy = abs(Vy_F[id]);
      Vz = abs(Vz_F[id]);
      Bx = abs(Bx_F[id]);
      By = abs(By_F[id]);
      Bz = abs(Bz_F[id]);
      
      energySpectrum_V_loc[id_k] += Vx*Vx;//+Vy*Vy+Vz*Vz; // UWAGA!
      energySpectrum_B_loc[id_k] += Bx*Bx+By*By+Bz*Bz;
      bin_counter_V_loc   [id_k] += 1;
      bin_counter_B_loc   [id_k] += 1;
    }
  }
  
  MPI_Reduce(energySpectrum_V_loc, energySpectrum_V, kmax, MPI_DOUBLE, MPI_SUM, 0, comm);
  MPI_Reduce(bin_counter_V_loc   , bin_counter_V   , kmax, MPI_INT   , MPI_SUM, 0, comm);
  MPI_Reduce(energySpectrum_B_loc, energySpectrum_B, kmax, MPI_DOUBLE, MPI_SUM, 0, comm);
  MPI_Reduce(bin_counter_B_loc   , bin_counter_B   , kmax, MPI_INT   , MPI_SUM, 0, comm);
  
  if(myRank == 0)
  {
    std::string file_name  = out_dir + "/spectrum_" + std::to_string(vti_count) + ".csv";
    std::ofstream os;
    
    os.open(file_name.c_str(), std::ios::out);
    if(!os){
      std::cout << "Cannot write header to file '" << file_name << "'!\n";
    }
    
    for(int ik = 1; ik < N_bin; ik++)
    {
      energySpectrum_V[ik] /= double(bin_counter_V[ik]); // divide by number of elements in bin to get mean values
      energySpectrum_V[ik] *= 4./3.*M_PI*del_k*del_k* ( (ik+1)*(ik+1)*(ik+1) - ik*ik*ik ); // get discrete Energy density
      
      
      energySpectrum_B[ik] /= double(bin_counter_B[ik]); // divide by number of elements in bin to get mean values
      energySpectrum_B[ik] *= 4./3.*M_PI*del_k*del_k* ( (ik+1)*(ik+1)*(ik+1) - ik*ik*ik ); // get discrete Energy density
      
      os << (ik+0.5)*del_k << ", " << energySpectrum_V[ik] << ", " << energySpectrum_B[ik] << std::endl;
    }
    
    double E_v = 0.;
    for(int ik = 1; ik < kmax; ik++)
    {
      E_v += energySpectrum_V[ik];
    }
    E_v *= dk;
    
    printf("Total Energy = %f\n", E_v);
    
    os.close();
  }MPI_Barrier(comm);
  
}
