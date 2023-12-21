#include "CSpecDyn.h"

CSpecDyn::CSpecDyn():
N(NUM), pdims(PDIMS), dt(DT), out_dir(OUT_DIR), out_interval(OUT_INTERVAL), end_simu(END_SIMU), L(LENGTH), nu(NU), eta(ETA), setup(SETUP)
{
  
  // init MPI
  MPI_Init(NULL, NULL);
 
  // init FFT  
  FFT = MikeFFT(N, pdims);
  
  // get sizes
  FFT.get_sizes_real(size_R, start_R, size_F, start_F);
  
  size_R_tot = size_R[0]*size_R[1]*size_R[2];
  size_F_tot = size_F[0]*size_F[1]*size_F[2];
  
  // get MPI info
  FFT.get_comm(comm);
  
  MPI_Comm_size(comm, &nprocs);
	MPI_Comm_rank(comm, &myRank);
  MPI_Cart_coords(comm, myRank, 2, mpi_coords);
  
  if(myRank==0)
  {
    printf("Real   : (%d, %d, %d)\n", size_R[0], size_R[1], size_R[2]);
    printf("Fourier: (%d, %d, %d)\n", size_F[0], size_F[1], size_F[2]);
  }
  
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
  
  Vx_F = FFT.malloc_F();
  Vy_F = FFT.malloc_F();
  Vz_F = FFT.malloc_F();
  Vx_F1 = FFT.malloc_F();
  Vy_F1 = FFT.malloc_F();
  Vz_F1 = FFT.malloc_F();
  Vx_F2 = FFT.malloc_F();
  Vy_F2 = FFT.malloc_F();
  Vz_F2 = FFT.malloc_F();
  Bx_F = FFT.malloc_F();
  By_F = FFT.malloc_F();
  Bz_F = FFT.malloc_F();
  Bx_F1 = FFT.malloc_F();
  By_F1 = FFT.malloc_F();
  Bz_F1 = FFT.malloc_F();
  Bx_F2 = FFT.malloc_F();
  By_F2 = FFT.malloc_F();
  Bz_F2 = FFT.malloc_F();
  RHS_Vx_F = FFT.malloc_F();
  RHS_Vy_F = FFT.malloc_F();
  RHS_Vz_F = FFT.malloc_F();
  RHS_Vx_F1 = FFT.malloc_F();
  RHS_Vy_F1 = FFT.malloc_F();
  RHS_Vz_F1 = FFT.malloc_F();
  RHS_Vx_F2 = FFT.malloc_F();
  RHS_Vy_F2 = FFT.malloc_F();
  RHS_Vz_F2 = FFT.malloc_F();
  RHS_Bx_F = FFT.malloc_F();
  RHS_By_F = FFT.malloc_F();
  RHS_Bz_F = FFT.malloc_F();
  RHS_Bx_F1 = FFT.malloc_F();
  RHS_By_F1 = FFT.malloc_F();
  RHS_Bz_F1 = FFT.malloc_F();
  RHS_Bx_F2 = FFT.malloc_F();
  RHS_By_F2 = FFT.malloc_F();
  RHS_Bz_F2 = FFT.malloc_F();
  
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
  Wx_F = FFT.malloc_F();
  Wy_F = FFT.malloc_F();
  Wz_F = FFT.malloc_F();
  Jx_F = FFT.malloc_F();
  Jy_F = FFT.malloc_F();
  Jz_F = FFT.malloc_F();
  
  Force_X = FFT.malloc_R();
  Force_Y = FFT.malloc_R();
  Force_Z = FFT.malloc_R();
  
  B0x = FFT.malloc_R();
  B0y = FFT.malloc_R();
  B0z = FFT.malloc_R();
  
  float_array        = (float*) malloc(sizeof(float)*size_R_tot);
  float_array_vector = (float*) malloc(sizeof(float)*size_R_tot*3);
  
  // create output directory
  if(myRank == 0)
	{
    std::string vti_dir = out_dir + "/vti";
    std::string scale_dir = out_dir + "/scales";
    std::string spectra_dir = out_dir + "/spectra";
    
    mkdir(out_dir.c_str(), 0777);
    mkdir(vti_dir.c_str(), 0777);
    mkdir(scale_dir.c_str(), 0777);
    mkdir(spectra_dir.c_str(), 0777);
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
  angle_eng = std::mt19937(myRank);
  angle = std::uniform_real_distribution<double>(0.,PI2);
  length_eng = std::mt19937(myRank+2);
  length = std::uniform_real_distribution<double>(-1.,1.);
  
  // Energy Spectrum
  energySpectrum_V = new double[N_bin];
  energySpectrum_V_loc = new double[N_bin];
  bin_counter_V    = new int[N_bin];
  bin_counter_V_loc    = new int[N_bin];
  energySpectrum_B = new double[N_bin];
  energySpectrum_B_loc = new double[N_bin];
  bin_counter_B    = new int[N_bin];
  bin_counter_B_loc    = new int[N_bin];
  
  // setup initial fields
  setup_k();
  setup_fields();
 
  if(RESTART_STEP > 0)
  {
    restart();
  }
  else
  {
    print_Energy();
    print();  
  }
  
}

void CSpecDyn::finalize()
{
  MPI_Finalize();
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
  
  double E0_V = 0.5;
  double E0_B = 0.25;
  double energy_b0_init = 5.;
  
  std::mt19937 eng(myRank);
  std::uniform_real_distribution<double> phi(0.,PI2);
  std::uniform_real_distribution<double> rand_real(1.,2.);
  double norm_loc = 0.;
  double norm;
  double s = 11./6.;
  
  double energy_V_loc = 0.;
  double energy_B_loc = 0.;
  double energy_B0_loc = 0.;
  double energy_V;
  double energy_B;
  double energy_B0;
  
  double diss_V_loc = 0.;
  double diss_B_loc = 0.;
  double diss_V;
  double diss_B;
  
  double Vx, Vy, Vz;
  double Bx, By, Bz;
  
  double norm_V, norm_B, norm_B0;
  
  double urms, urms_loc;
  
  switch(setup)
  {
    case 0:
      if(myRank==0){printf("0\n");}
      
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
      if(myRank==0){printf("1\n");}
    
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
        
        Bx_R[id] = ( -2.*sin(2*y_val) + sin(z_val) )*beta;
        By_R[id] = (  2.*sin(  x_val) + sin(z_val) )*beta;
        Bz_R[id] = (     sin(  x_val) + sin(y_val) )*beta;
    
      }}}
      
      fFFT(Vx_R, Vy_R, Vz_R, Vx_F, Vy_F, Vz_F);
      fFFT(Bx_R, By_R, Bz_R, Bx_F, By_F, Bz_F);
      
      break;
      
    case 2:
      if(myRank==0){printf("2\n");}
      
    /** random with energy spectrum normalized to desired kinetic energy **/
    
      // Zufallsfeld in R
      for(int id = 0; id < size_R_tot; id++)
      {    
        Vx_R[id] = rand_real(eng);
        Vy_R[id] = rand_real(eng);
        Vz_R[id] = rand_real(eng);
        
        Bx_R[id] = rand_real(eng);
        By_R[id] = rand_real(eng);
        Bz_R[id] = rand_real(eng);
      }
      fFFT(Vx_R, Vy_R, Vz_R, Vx_F, Vy_F, Vz_F);
      fFFT(Bx_R, By_R, Bz_R, Bx_F, By_F, Bz_F);
      
      // Amplituden für Energie-Spektrum
      for(int id = 0; id < size_F_tot; id++)
      {
        double A = sqrt( 1./pow(1+k2[id],s) );
        
        //~ //double A = sqrt( k2[id]*k2[id] *exp( -2. * k2[id]/(k0*k0) ) );
        
        if(int(round(k2[id])) != 0)
        {
          Vx_F[id] *= A/abs(Vx_F[id]);
          Vy_F[id] *= A/abs(Vy_F[id]);
          Vz_F[id] *= A/abs(Vz_F[id]);
          Bx_F[id] *= A/abs(Bx_F[id]);
          By_F[id] *= A/abs(By_F[id]);
          Bz_F[id] *= A/abs(Bz_F[id]);
        }
        else
        {
          Vx_F[id] = 0.;
          Vy_F[id] = 0.;
          Vz_F[id] = 0.;
          Bx_F[id] = 0.;
          By_F[id] = 0.;
          Bz_F[id] = 0.;
        }
      }
      
      projection(Vx_F, Vy_F, Vz_F);
      projection(Bx_F, By_F, Bz_F);
      
      // Energie normieren
      for(int ix = 0; ix < size_F[0]; ix++){
      for(int iy = 0; iy < size_F[1]; iy++){
      for(int iz = 0; iz < size_F[2]; iz++){
          
        int id = ix * size_F[1]*size_F[2] + iy * size_F[2] + iz;
        
        Vx = abs(Vx_F[id]);
        Vy = abs(Vy_F[id]);
        Vz = abs(Vz_F[id]);
        Bx = abs(Bx_F[id]);
        By = abs(By_F[id]);
        Bz = abs(Bz_F[id]);
        
        energy_V_loc += (Vx*Vx+Vy*Vy+Vz*Vz);
        energy_B_loc += (Bx*Bx+By*By+Bz*Bz);

      }}}
      
      energy_V_loc *= 0.5/double(N*N*N); // Ortsmittelung und 0.5 aus Definition der Energie/Definition Energy Spectrum?
      energy_V_loc *= 1. /double(N*N*N); // wg Fourier Space
      MPI_Allreduce(&energy_V_loc, &energy_V, 1, MPI_DOUBLE, MPI_SUM, comm);
      
      energy_B_loc *= 0.5/double(N*N*N);
      energy_B_loc *= 1. /double(N*N*N);
      MPI_Allreduce(&energy_B_loc, &energy_B, 1, MPI_DOUBLE, MPI_SUM, comm);
      
      norm_V = sqrt(E0_V/energy_V);    
      norm_B = sqrt(E0_B/energy_B);    
      
      for(int id = 0; id < size_F_tot; id++)
      {    
        Vx_F[id] *= norm_V;
        Vy_F[id] *= norm_V;
        Vz_F[id] *= norm_V;
        Bx_F[id] *= norm_B;
        By_F[id] *= norm_B;
        Bz_F[id] *= norm_B;
        
        #ifdef NS
        Bx_F[id] = 0.;
        By_F[id] = 0.;
        Bz_F[id] = 0.;
        #endif
      }
        
      break;
      
    case 3:
      if(myRank==0){printf("3\n");}
    
      // read .dat files
      read_binary();
      
      // normalize fields to u_rms = 1
      bFFT(Vx_F, Vy_F, Vz_F, Vx_R, Vy_R, Vz_R);
      bFFT(Bx_F, By_F, Bz_F, Bx_R, By_R, Bz_R);
      
      urms_loc = 0.;
      
      for(int id = 0; id < size_R_tot; id++)
      {
        
        Vx = Vx_R[id];
        Vy = Vy_R[id];
        Vz = Vz_R[id];
        
        urms_loc += (Vx*Vx+Vy*Vy+Vz*Vz);

      }
      
      MPI_Allreduce(&urms_loc, &urms, 1, MPI_DOUBLE, MPI_SUM, comm);
      
      urms *= 1. /double(N*N*N);
      urms = sqrt(1./3. * urms);
      urms *= 0.5; // norm to urms = 0.5-> T = 1 for L = 0.5
      
      norm = 1./urms;
      
      if(myRank==0){printf("urms = %f\n", urms);}
      
      for(int id = 0; id < size_R_tot; id++)
      {    
        Vx_R[id] *= norm;
        Vy_R[id] *= norm;
        Vz_R[id] *= norm;
        Bx_R[id] *= norm;
        By_R[id] *= norm;
        Bz_R[id] *= norm;
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
  
  /** Background Field **/
  if(BACKGROUND)
  {
   /** B0 = B0*ez **/
    for(int id = 0; id < size_R_tot; id++)
    {    
      B0x[id] = 0.;
      B0y[id] = 0.;
      B0z[id] = 1.;
    }
    
    // Energie in B0 auf vorgegeben Wert setzen!
    for(int id = 0; id < size_R_tot; id++){
        
        Bx = abs(B0x[id]);
        By = abs(B0y[id]);
        Bz = abs(B0z[id]);
        
        energy_B0_loc += (Bx*Bx+By*By+Bz*Bz);

    }
    
    energy_B0_loc *= 0.5/double(N*N*N); // sind in RealSpace!
    MPI_Allreduce(&energy_B0_loc, &energy_B0, 1, MPI_DOUBLE, MPI_SUM, comm);
    
    norm_B0 = sqrt(energy_b0_init/energy_B0);      
    
    for(int id = 0; id < size_R_tot; id++)
    {    
      B0x[id] *= norm_B0;
      B0y[id] *= norm_B0;
      B0z[id] *= norm_B0;
    }
    
  }
  else
  {
    /** B0 = 0 **/
    for(int id = 0; id < size_R_tot; id++)
    {    
      B0x[id] = 0.;
      B0y[id] = 0.;
      B0z[id] = 0.;
    }
  }
  
  // set initial forcing to zero
  for(size_t id = 0; id < size_R_tot; id++)
  {
    Force_X[id] = 0.;
    Force_Y[id] = 0.;
    Force_Z[id] = 0.;
  }
}

void CSpecDyn::execute()
{
  double start_time = MPI_Wtime();
  double out_time = fmod(time,out_interval);
  
  // geometric output series
  double r = 1.4919;
  double a = 0.00067;
  
  out_time = a * r; 
  
  while(time < end_simu)
  {
    
    time_step();
    out_time += dt;
    
    print_Energy();
    
    if(out_time > out_interval)
    // //if(time > out_time) // log output
    {
      print();
      out_time -= out_interval;
      // out_time = a * pow(r,print_count); // log output
    }
    
  }
  
  double end_time = MPI_Wtime() - start_time;
  if(myRank==0)
  {
    printf("Execution Time: %f [s]\n", end_time);
  }
}

void CSpecDyn::time_step()
{
  
  set_dt();
  
  if(FORCING==0)
  {
    // nothing
  }
  else if(FORCING==1)
  {
    Alvelius();
  }
  
  double del_t;
  
  // step 1
  del_t = 1.;
  
  calc_RHS(RHS_Vx_F , RHS_Vy_F , RHS_Vz_F , Vx_F , Vy_F , Vz_F
          ,RHS_Bx_F , RHS_By_F , RHS_Bz_F , Bx_F , By_F , Bz_F, del_t);
  
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
          ,RHS_Bx_F1, RHS_By_F1, RHS_Bz_F1, Bx_F1, By_F1, Bz_F1, del_t);
   
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
          ,RHS_Bx_F2, RHS_By_F2, RHS_Bz_F2, Bx_F2, By_F2, Bz_F2, del_t);
  
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

// calc Energy and Dissipation
void CSpecDyn::calc_Energy(double& energy_V, double& diss_V)
{
  
    double energy_V_loc = 0.;
    double diss_V_loc = 0.;
    double Vx, Vy, Vz;
  
    for(int ix = 0; ix < size_F[0]; ix++){
    for(int iy = 0; iy < size_F[1]; iy++){
    for(int iz = 0; iz < size_F[2]; iz++){
      
      // globale id  
      int id = ix * size_F[1]*size_F[2] + iy * size_F[2] + iz;
      
      Vx = abs(Vx_F[id]);
      Vy = abs(Vy_F[id]);
      Vz = abs(Vz_F[id]);
      
      energy_V_loc += (Vx*Vx+Vy*Vy+Vz*Vz);
      diss_V_loc   += (Vx*Vx+Vy*Vy+Vz*Vz) * k2[id];

    }}}
  
    // berechnen globale Größen
    energy_V_loc *= 0.5/double(N*N*N); // 1/N^3: Ortsmittelung, 0.5: aus Definition der Energie/Definition Energy Spectrum
    energy_V_loc *= 1. /double(N*N*N); // wg DFT Normierung
    MPI_Allreduce(&energy_V_loc, &energy_V, 1, MPI_DOUBLE, MPI_SUM, comm);

    diss_V_loc *= 0.5/double(N*N*N);
    diss_V_loc *= 1. /double(N*N*N);
    MPI_Allreduce(&diss_V_loc, &diss_V, 1, MPI_DOUBLE, MPI_SUM, comm);
    diss_V *= 2.*nu;
    
}

void CSpecDyn::calc_Energy(double& energy_V, double& diss_V, double& energy_B, double& diss_B)
{
  
    double energy_V_loc = 0.;
    double diss_V_loc = 0.;
    double Vx, Vy, Vz;
    double energy_B_loc = 0.;
    double diss_B_loc = 0.;
    double Bx, By, Bz;
    
    double hs; // factor because of hermitian symmetry in z
  
    for(int ix = 0; ix < size_F[0]; ix++){
    for(int iy = 0; iy < size_F[1]; iy++){
    for(int iz = 0; iz < size_F[2]; iz++){
      
      // globale id  
      int id = ix * size_F[1]*size_F[2] + iy * size_F[2] + iz;
      
      Vx = abs(Vx_F[id]);
      Vy = abs(Vy_F[id]);
      Vz = abs(Vz_F[id]);
      
      energy_V_loc += (Vx*Vx+Vy*Vy+Vz*Vz);
      diss_V_loc   += (Vx*Vx+Vy*Vy+Vz*Vz) * k2[id];

      Bx = abs(Bx_F[id]);
      By = abs(By_F[id]);
      Bz = abs(Bz_F[id]);
      
      energy_B_loc += (Bx*Bx+By*By+Bz*Bz);
      diss_B_loc   += (Bx*Bx+By*By+Bz*Bz) * k2[id];

    }}}
  
    // berechnen globale Größen
    energy_V_loc *= 0.5/double(N*N*N); // 1/N^3: Ortsmittelung, 0.5: zählen nur positive Moden in der Definition in Spektralraum
    energy_V_loc *= 1. /double(N*N*N); // wg DFT Normierung
    MPI_Allreduce(&energy_V_loc, &energy_V, 1, MPI_DOUBLE, MPI_SUM, comm);

    diss_V_loc *= 0.5/double(N*N*N);
    diss_V_loc *= 1. /double(N*N*N);
    MPI_Allreduce(&diss_V_loc, &diss_V, 1, MPI_DOUBLE, MPI_SUM, comm);
    diss_V *= 2.*nu;
    
    energy_B_loc *= 0.5/double(N*N*N);
    energy_B_loc *= 1. /double(N*N*N);
    MPI_Allreduce(&energy_B_loc, &energy_B, 1, MPI_DOUBLE, MPI_SUM, comm);

    diss_B_loc *= 0.5/double(N*N*N);
    diss_B_loc *= 1. /double(N*N*N);
    MPI_Allreduce(&diss_B_loc, &diss_B, 1, MPI_DOUBLE, MPI_SUM, comm);
    diss_B *= 2.*eta;
    
}

void CSpecDyn::set_dt()
{
  bFFT(Vx_F, Vy_F, Vz_F, Vx_R, Vy_R, Vz_R);
  
  double vmax_loc = 0.;
  for(int id = 0; id < size_R_tot; id++)
  {
    double Vx = Vx_R[id];
    double Vy = Vy_R[id];
    double Vz = Vz_R[id];
    
    vmax_loc = std::max( {vmax_loc, Vx*Vx+Vy*Vy+Vz*Vz} );
  }
  
  double vmax;
  MPI_Allreduce(&vmax_loc, &vmax, 1, MPI_DOUBLE, MPI_MAX, comm);
  vmax = sqrt(vmax);
  
  double dt_adv = CFL_ADV * dx / vmax;
  double dt_dif = CFL_DIF * dx * dx / nu;
  dt = std::min(dt_adv, dt_dif);
  
  fFFT(Vx_R, Vy_R, Vz_R, Vx_F, Vy_F, Vz_F);
}

void CSpecDyn::Alvelius()
{
  
  //double P = M_PI*M_PI;
  //~ double P = 7.64;
  //~ int kf = 1;
  //~ double Nf = 1.;
  
  //~ for(int ix = 0; ix<size_F[0]; ix++){
  //~ for(int iy = 0; iy<size_F[1]; iy++){
  //~ for(int iz = 0; iz<size_F[2]; iz++){
    
    //~ int id = ix * size_F[1]*size_F[2] + iy * size_F[2] + iz;
    
    //~ double k = sqrt(k2[id]);
    //~ int k_int  = int(round(k));
    
    //~ if( fabs(kf-k) < 0.001)
    //~ {
    
      //~ double kxx = kx[ix];
      //~ double kyy = ky[iy];
      //~ double kzz = kz[iz];
      //~ double k   = sqrt(k2[id]);
      
      //~ double phi   = atan2(kxx,kzz);
      //~ double theta = atan2(hypot(kxx,kzz), kyy);
      
      //~ double e1[3] = {+           sin(phi), -           cos(phi), 0.         };
      //~ double e2[3] = {-cos(theta)*cos(phi), -cos(theta)*sin(phi), +sin(theta)};
      
      //~ CX xi_1 = Vx_F[id]*e1[0] + Vy_F[id]*e1[1] + Vz_F[id]*e1[2];
      //~ CX xi_2 = Vx_F[id]*e2[0] + Vy_F[id]*e2[1] + Vz_F[id]*e2[2];
      
      //~ double alp = angle(angle_eng);
      //~ double psi = angle(angle_eng);
      //~ double gA = sin(2.*alp);
      //~ double gB = cos(2.*alp);
      
      //~ double theta_1 = atan2( gA * xi_1.real() + gB * ( sin(psi) * xi_2.imag() + cos(psi) * xi_2.real() ) ,
                             //~ -gA * xi_1.imag() + gB * ( sin(psi) * xi_2.real() - cos(psi) * xi_2.imag() ) );
                           
      //~ double theta_2 = theta_1 + psi;
      
      //~ double F =  P/(dt); // delta Forcing! (3 Moden im Band kf=1 oder kf=2)
      
      //~ CX A = sqrt( F ) * exp(IM*theta_1) * gA;
      //~ CX B = sqrt( F ) * exp(IM*theta_2) * gB;
      
      //~ // 2.7 normiert auf eps=1 und sqrt(eps_0) setzt den gewünschten eps-Wert
      //~ double factor = N*N*N * 2.7;   // NS
      //double factor = N*N*N*2.7* sqrt(8.); // MHD
      
      //~ Force_X[id] = factor*(A * e1[0]  + B * e2[0]); 
      //~ Force_Y[id] = factor*(A * e1[1]  + B * e2[1]); 
      //~ Force_Z[id] = factor*(A * e1[2]  + B * e2[2]); 
    
    //~ }
    //~ else
    //~ {
      //~ Force_X[id] = 0.;
      //~ Force_Y[id] = 0.;
      //~ Force_Z[id] = 0.;
    //~ }
  //~ }}}
  
  //~ // only keep real part
  //~ bFFT(Force_X, Force_Y, Force_Z, Jx_R, Jy_R, Jz_R);
  
  //~ for(int id = 0; id<size_F_tot; id++)
  //~ {
    //~ Jx_R[id] = CX( Jx_R[id].real() );
    //~ Jy_R[id] = CX( Jy_R[id].real() );
    //~ Jz_R[id] = CX( Jz_R[id].real() );
  //~ }
  
  //~ fFFT(Jx_R, Jy_R, Jz_R, Force_X, Force_Y, Force_Z);
  
}

void CSpecDyn::calc_RHS(CX* RHSV_X, CX* RHSV_Y, CX* RHSV_Z, CX* V_X, CX* V_Y, CX* V_Z,
                        CX* RHSB_X, CX* RHSB_Y, CX* RHSB_Z, CX* B_X, CX* B_Y, CX* B_Z,
                        double del_t)
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
    
    RHS_Vx_R[id] = Vy_R[id]*Wz_R[id]-Vz_R[id]*Wy_R[id] + (Jy_R[id]*(Bz_R[id]+B0z[id])-Jz_R[id]*(By_R[id]+B0y[id]) );
    RHS_Vy_R[id] = Vz_R[id]*Wx_R[id]-Vx_R[id]*Wz_R[id] + (Jz_R[id]*(Bx_R[id]+B0x[id])-Jx_R[id]*(Bz_R[id]+B0z[id]) );
    RHS_Vz_R[id] = Vx_R[id]*Wy_R[id]-Vy_R[id]*Wx_R[id] + (Jx_R[id]*(By_R[id]+B0y[id])-Jy_R[id]*(Bx_R[id]+B0x[id]) );
    
  }
  
  // RHS_B = VxB
  for(int id = 0; id < size_R_tot; id++){
    
    RHS_Bx_R[id] = Vy_R[id]*(Bz_R[id]+B0z[id])-Vz_R[id]*(By_R[id]+B0y[id]);
    RHS_By_R[id] = Vz_R[id]*(Bx_R[id]+B0x[id])-Vx_R[id]*(Bz_R[id]+B0z[id]);
    RHS_Bz_R[id] = Vx_R[id]*(By_R[id]+B0y[id])-Vy_R[id]*(Bx_R[id]+B0x[id]);
    
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
  
  /************ ALVELIUS - 1999 **********/
  if(setup==2)
  {
    for(int id = 0; id < size_F_tot; id++){
     
      RHSV_X[id] += Force_X[id];
      RHSV_Y[id] += Force_Y[id];
      RHSV_Z[id] += Force_Z[id];
      
    }
  }
  /***************************************/
  
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
  
  std::string file_name  = out_dir + "/vti/step_" + std::to_string(print_count) + ".vti";
  std::ofstream os;
  
  long N_l = N;
  long offset = 0;
	long N_tot = N_l*N_l*N_l;
	long N_bytes_scalar  =   N_tot * sizeof(float);
	long N_bytes_vector  = 3*N_tot * sizeof(float);
  long bin_size_scalar = N_bytes_scalar + sizeof(uint64_t);// 2nd term is the size of the the leading integer announcing the numbers n the data chunk
  long bin_size_vector = N_bytes_vector + sizeof(uint64_t);
  
  // header
  if(myRank==0)
  {
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
         << "\" Direction=\"0 0 1 0 1 0 1 0 0\">" << std::endl; // FORTRAN -> C order
    
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
  
  fFFT(Vx_R, Vy_R, Vz_R, Vx_F, Vy_F, Vz_F);
  fFFT(Bx_R, By_R, Bz_R, Bx_F, By_F, Bz_F);
}

void CSpecDyn::print_mpi_vector(double* field_X, double* field_Y, double* field_Z, long& N_bytes_vector, const char* file_name)
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

void CSpecDyn::print_mpi_scalar(double* field, long& N_bytes_scalar, const char* file_name)
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
  // Hou&Li
  double kmax = N/2*dk;
  
 for(int id = 0; id < size_F_tot; id++){
   
    double k = sqrt(k2[id]);
    double filter = exp(-36. * pow(k/kmax, 36) );
   
    fieldX[id] *= filter;
    fieldY[id] *= filter;
    fieldZ[id] *= filter;
    
  }
    
}

void CSpecDyn::print_EnergySpectrum()
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
      
      energySpectrum_V_loc[id_k] += Vx*Vx+Vy*Vy+Vz*Vz;
      energySpectrum_B_loc[id_k] += Bx*Bx+By*By+Bz*Bz;
      bin_counter_V_loc   [id_k] += 1;
      bin_counter_B_loc   [id_k] += 1;
    }
  }
  
  MPI_Reduce(energySpectrum_V_loc, energySpectrum_V, N_bin, MPI_DOUBLE, MPI_SUM, 0, comm);
  MPI_Reduce(bin_counter_V_loc   , bin_counter_V   , N_bin, MPI_INT   , MPI_SUM, 0, comm);
  MPI_Reduce(energySpectrum_B_loc, energySpectrum_B, N_bin, MPI_DOUBLE, MPI_SUM, 0, comm);
  MPI_Reduce(bin_counter_B_loc   , bin_counter_B   , N_bin, MPI_INT   , MPI_SUM, 0, comm);
  
  if(myRank == 0)
  {
    std::string file_name  = out_dir + "/spectra/spectrum_" + std::to_string(print_count) + ".csv";
    std::ofstream os;
    
    os.open(file_name.c_str(), std::ios::out);
    if(!os){
      std::cout << "Cannot write header to file '" << file_name << "'!\n";
    }
    
    for(int ik = 1; ik < N_bin; ik++)
    {
      energySpectrum_V[ik] /= double(bin_counter_V[ik]); // divide by number of elements in bin to get mean values
      energySpectrum_V[ik] *= 0.5 * 4./3.*M_PI*del_k*del_k* ( (ik+1)*(ik+1)*(ik+1) - ik*ik*ik ); // get discrete Energy density
      energySpectrum_B[ik] /= double(bin_counter_B[ik]);
      energySpectrum_B[ik] *= 0.5 * 4./3.*M_PI*del_k*del_k* ( (ik+1)*(ik+1)*(ik+1) - ik*ik*ik );
      
      os << (ik+0.5)*del_k << ", " << energySpectrum_V[ik] << ", " << energySpectrum_B[ik] << std::endl;
    }
    
    os.close();
  }MPI_Barrier(comm);
  
}

//~ void CSpecDyn::print_Energy()
//~ {
  //~ // Compute mean Energy and Dissipation in Fourier Space
  //~ /** Energie und Dissipation direkt **/
  //~ double energy_V_loc = 0.;
  //~ double energy_B_loc = 0.;
  //~ double energy_V;
  //~ double energy_B;
  
  //~ double diss_V_loc = 0.;
  //~ double diss_B_loc = 0.;
  //~ double diss_V;
  //~ double diss_B;
  
  //~ double ens_V_loc = 0.;
  //~ double ens_B_loc = 0.;
  //~ double ens_V;
  //~ double ens_B;
  
  //~ CX Hc_loc = 0.;
  //~ double Hcr_loc = 0.;
  //~ double Hci_loc = 0.;
  //~ double Hcr;
  //~ double Hci;
  
  //~ double Vx, Vy, Vz;
  //~ double Bx, By, Bz;
  //~ double Wx, Wy, Wz;
  //~ double Jx, Jy, Jz;
  
  //~ for(int ix = 0; ix < size_F[0]; ix++){
  //~ for(int iy = 0; iy < size_F[1]; iy++){
  //~ for(int iz = 0; iz < size_F[2]; iz++){
      
    //~ int id = ix * size_F[1]*size_F[2] + iy * size_F[2] + iz;
    
    //~ Vx = abs(Vx_F[id]);
    //~ Vy = abs(Vy_F[id]);
    //~ Vz = abs(Vz_F[id]);
    //~ Bx = abs(Bx_F[id]);
    //~ By = abs(By_F[id]);
    //~ Bz = abs(Bz_F[id]);
    //~ Wx = abs(ky[iy] * Vz_F[id] - kz[iz] * Vy_F[id]);
    //~ Wy = abs(kz[iz] * Vx_F[id] - kx[ix] * Vz_F[id]);
    //~ Wz = abs(kx[ix] * Vy_F[id] - ky[iy] * Vx_F[id]);
    //~ Jx = abs(ky[iy] * Bz_F[id] - kz[iz] * By_F[id]);
    //~ Jy = abs(kz[iz] * Bx_F[id] - kx[ix] * Bz_F[id]);
    //~ Jz = abs(kx[ix] * By_F[id] - ky[iy] * Bx_F[id]);
    
    //~ energy_V_loc += (Vx*Vx+Vy*Vy+Vz*Vz);
    //~ energy_B_loc += (Bx*Bx+By*By+Bz*Bz);
    
    //~ diss_V_loc   += (Vx*Vx+Vy*Vy+Vz*Vz) * k2[id];
    //~ diss_B_loc   += (Bx*Bx+By*By+Bz*Bz) * k2[id];
    
    //~ ens_V_loc    += (Wx*Wx+Wy*Wy+Wz*Wz);
    //~ ens_B_loc    += (Jx*Jx+Jy*Jy+Jz*Jz);
    
    //~ Hc_loc += Vx_F[id]*std::conj(Bx_F[id])+Vy_F[id]*std::conj(By_F[id])+Vz_F[id]*std::conj(Bz_F[id]);

  //~ }}}
  
  //~ double N3_inv = 1./double(N*N*N);
  
  //~ energy_V_loc *= N3_inv; // Ortsmittelung
  //~ energy_V_loc *= N3_inv; // DFT
  //~ MPI_Reduce(&energy_V_loc, &energy_V, 1, MPI_DOUBLE, MPI_SUM, 0, comm);
  //~ energy_V *= 0.5;
  
  //~ energy_B_loc *= N3_inv;
  //~ energy_B_loc *= N3_inv;
  //~ MPI_Reduce(&energy_B_loc, &energy_B, 1, MPI_DOUBLE, MPI_SUM, 0, comm);
  //~ energy_B *= 0.5;
  
  //~ diss_V_loc *= N3_inv;
  //~ diss_V_loc *= N3_inv;
  //~ MPI_Reduce(&diss_V_loc, &diss_V, 1, MPI_DOUBLE, MPI_SUM, 0, comm);
  //~ diss_V *= nu;
  
  //~ diss_B_loc *= N3_inv;
  //~ diss_B_loc *= N3_inv;
  //~ MPI_Reduce(&diss_B_loc, &diss_B, 1, MPI_DOUBLE, MPI_SUM, 0, comm);
  //~ diss_B *= eta;
  
  //~ ens_V_loc *= N3_inv;
  //~ ens_V_loc *= N3_inv;
  //~ MPI_Reduce(&ens_V_loc, &ens_V, 1, MPI_DOUBLE, MPI_SUM, 0, comm);
  //~ ens_V *= 1.;
  
  //~ ens_B_loc *= N3_inv;
  //~ ens_B_loc *= N3_inv;
  //~ MPI_Reduce(&ens_B_loc, &ens_B, 1, MPI_DOUBLE, MPI_SUM, 0, comm);
  //~ ens_B *= 1.;
  
  //~ Hcr_loc = Hc_loc.real();
  //~ Hcr_loc *= N3_inv;
  //~ Hcr_loc *= N3_inv;
  //~ MPI_Reduce(&Hcr_loc, &Hcr, 1, MPI_DOUBLE, MPI_SUM, 0, comm);
  //~ Hcr *= 1.;
  
  //~ Hci_loc = Hc_loc.imag();
  //~ Hci_loc *= N3_inv;
  //~ Hci_loc *= N3_inv;
  //~ MPI_Reduce(&Hci_loc, &Hci, 1, MPI_DOUBLE, MPI_SUM, 0, comm);
  //~ Hci *= 1.;
  
  //~ // print
  //~ if(myRank == 0)
  //~ {
    //~ std::string file_name  = out_dir + "/energy.csv";
    //~ std::ofstream os;
    
    //~ os.open(file_name.c_str(), std::ios::out | std::ios::app);
    //~ if(!os){
      //~ std::cout << "Cannot write header to file '" << file_name << "'!\n";
    //~ }
      
    //~ os << time << ", " << energy_V << ", " << energy_B  << ", " << diss_V << ", " <<  diss_B 
               //~ << ", " << ens_V    << ", " <<  ens_B    << ", " <<  Hcr    << ", " <<  Hci    << std::endl;
    
    //~ os.close();
  //~ }MPI_Barrier(comm);
  
//~ }

void CSpecDyn::print_Energy()
{
  // Compute mean Energy and Dissipation in Fourier Space
  /** Energie und Dissipation direkt **/
  double energy_V_loc = 0.;
  double energy_B_loc = 0.;
  double energy_V;
  double energy_B;
  
  double diss_V_loc = 0.;
  double diss_B_loc = 0.;
  double diss_V;
  double diss_B;
  
  double ens_V_loc = 0.;
  double ens_B_loc = 0.;
  double ens_V;
  double ens_B;
  
  CX Hc_loc = 0.;
  double Hcr_loc = 0.;
  double Hci_loc = 0.;
  double Hcr;
  double Hci;
  
  double Vx, Vy, Vz;
  double Bx, By, Bz;
  double Wx, Wy, Wz;
  double Jx, Jy, Jz;
  
  double hs; // factor because of hermitian symmetry in z
  
  for(int ix = 0; ix < size_F[0]; ix++){
  for(int iy = 0; iy < size_F[1]; iy++){
  for(int iz = 0; iz < size_F[2]; iz++){
      
    int id = ix * size_F[1]*size_F[2] + iy * size_F[2] + iz;
    
    int kz_id = int(kz[iz]/dk);
    if( 0 < kz_id && kz_id < N/2 )
    {
      hs = 2.;
    }
    else
    {
      hs = 1.;
    }
    
    Vx = abs(Vx_F[id]);
    Vy = abs(Vy_F[id]);
    Vz = abs(Vz_F[id]);
    Bx = abs(Bx_F[id]);
    By = abs(By_F[id]);
    Bz = abs(Bz_F[id]);
    Wx = abs(ky[iy] * Vz_F[id] - kz[iz] * Vy_F[id]);
    Wy = abs(kz[iz] * Vx_F[id] - kx[ix] * Vz_F[id]);
    Wz = abs(kx[ix] * Vy_F[id] - ky[iy] * Vx_F[id]);
    Jx = abs(ky[iy] * Bz_F[id] - kz[iz] * By_F[id]);
    Jy = abs(kz[iz] * Bx_F[id] - kx[ix] * Bz_F[id]);
    Jz = abs(kx[ix] * By_F[id] - ky[iy] * Bx_F[id]);
    
    energy_V_loc += hs*(Vx*Vx+Vy*Vy+Vz*Vz);
    energy_B_loc += hs*(Bx*Bx+By*By+Bz*Bz);
    
    diss_V_loc   += hs*(Vx*Vx+Vy*Vy+Vz*Vz) * k2[id];
    diss_B_loc   += hs*(Bx*Bx+By*By+Bz*Bz) * k2[id];
    
    ens_V_loc    += hs*(Wx*Wx+Wy*Wy+Wz*Wz);
    ens_B_loc    += hs*(Jx*Jx+Jy*Jy+Jz*Jz);
    
    Hc_loc += hs*(Vx_F[id]*std::conj(Bx_F[id])+Vy_F[id]*std::conj(By_F[id])+Vz_F[id]*std::conj(Bz_F[id]));

  }}}
  
  double N3_inv = 1./double(N*N*N);
  
  energy_V_loc *= 0.5*N3_inv; // Ortsmittelung und 0.5 aus Definition der Energie/Definition Energy Spectrum?
  energy_V_loc *= N3_inv; // wg Fourier Space
  MPI_Reduce(&energy_V_loc, &energy_V, 1, MPI_DOUBLE, MPI_SUM, 0, comm);
  
  energy_B_loc *= 0.5*N3_inv;
  energy_B_loc *= N3_inv;
  MPI_Reduce(&energy_B_loc, &energy_B, 1, MPI_DOUBLE, MPI_SUM, 0, comm);
  
  diss_V_loc *= 0.5 *N3_inv;
  diss_V_loc *= N3_inv;
  MPI_Reduce(&diss_V_loc, &diss_V, 1, MPI_DOUBLE, MPI_SUM, 0, comm);
  diss_V *= 2. * nu;
  
  diss_B_loc *= eta*N3_inv;
  diss_B_loc *= N3_inv;
  MPI_Reduce(&diss_B_loc, &diss_B, 1, MPI_DOUBLE, MPI_SUM, 0, comm);
  
  ens_V_loc *= 2. *N3_inv;
  ens_V_loc *= N3_inv;
  MPI_Reduce(&ens_V_loc, &ens_V, 1, MPI_DOUBLE, MPI_SUM, 0, comm);
  
  ens_B_loc *= 2. *N3_inv;
  ens_B_loc *= N3_inv;
  MPI_Reduce(&ens_B_loc, &ens_B, 1, MPI_DOUBLE, MPI_SUM, 0, comm);
  
  Hcr_loc = Hc_loc.real();
  Hcr_loc *= N3_inv;
  Hcr_loc *= N3_inv;
  MPI_Reduce(&Hcr_loc, &Hcr, 1, MPI_DOUBLE, MPI_SUM, 0, comm);
  Hcr *= 1.;
  
  Hci_loc = Hc_loc.imag();
  Hci_loc *= N3_inv;
  Hci_loc *= N3_inv;
  MPI_Reduce(&Hci_loc, &Hci, 1, MPI_DOUBLE, MPI_SUM, 0, comm);
  Hci *= 1.;
  
  // print
  if(myRank == 0)
  {
    std::string file_name  = out_dir + "/energy.csv";
    std::ofstream os;
    
    os.open(file_name.c_str(), std::ios::out | std::ios::app);
    if(!os){
      std::cout << "Cannot write header to file '" << file_name << "'!\n";
    }
      
    os << time << ", " << energy_V << ", " << energy_B  << ", " << diss_V << ", " <<  diss_B 
               << ", " << ens_V    << ", " <<  ens_B    << ", " <<  Hcr    << ", " <<  Hci    << std::endl;
    
    os.close();
  }MPI_Barrier(comm);
  
}

void CSpecDyn::print_scales()
{
  // energy
  double energy_V_loc = 0.;
  double energy_B_loc = 0.;
  double energy_V;
  double energy_B;
  
  // dissipation
  double diss_V_loc = 0.;
  double diss_B_loc = 0.;
  double diss_V;
  double diss_B;
  
  // integral scale
  double L_V_loc = 0.;
  double L_B_loc = 0.;
  double L_V;
  double L_B;
  
  // further scales
  double v0, b0;             // RMS
  double lambda_V, lambda_B; // Taylor micro scale
  double Re_V, Re_B;         // Taylor scale Re
  double T_V, T_B;           // LE turnover time
  double tau_V, tau_B;       // Kolmogorov time scale
  double eta_V, eta_B;       // Kolmogorov length scale
  
  double Vx, Vy, Vz;
  double Bx, By, Bz;
  
  for(int ix = 0; ix < size_F[0]; ix++){
  for(int iy = 0; iy < size_F[1]; iy++){
  for(int iz = 0; iz < size_F[2]; iz++){
      
    int id = ix * size_F[1]*size_F[2] + iy * size_F[2] + iz;
    
    Vx = abs(Vx_F[id]);
    Vy = abs(Vy_F[id]);
    Vz = abs(Vz_F[id]);
    Bx = abs(Bx_F[id]);
    By = abs(By_F[id]);
    Bz = abs(Bz_F[id]);
    
    energy_V_loc += (Vx*Vx+Vy*Vy+Vz*Vz);
    energy_B_loc += (Bx*Bx+By*By+Bz*Bz);
    
    diss_V_loc += (Vx*Vx+Vy*Vy+Vz*Vz) * k2[id];
    diss_B_loc += (Bx*Bx+By*By+Bz*Bz) * k2[id];
    
    if(id!=0) // avoid dividing by zero
    {
      L_V_loc += (Vx*Vx+Vy*Vy+Vz*Vz) / sqrt(k2[id]);
      L_B_loc += (Bx*Bx+By*By+Bz*Bz) / sqrt(k2[id]);
    }

  }}}
  
  energy_V_loc *= 0.5/double(N*N*N); // Ortsmittelung und 0.5 aus Definition der Energie/Definition Energy Spectrum?
  energy_V_loc *= 1. /double(N*N*N); // wg Fourier Space
  MPI_Reduce(&energy_V_loc, &energy_V, 1, MPI_DOUBLE, MPI_SUM, 0, comm);
  
  energy_B_loc *= 0.5/double(N*N*N);
  energy_B_loc *= 1. /double(N*N*N);
  MPI_Reduce(&energy_B_loc, &energy_B, 1, MPI_DOUBLE, MPI_SUM, 0, comm);
  
  diss_V_loc *= 0.5/double(N*N*N);
  diss_V_loc *= 1. /double(N*N*N);
  MPI_Reduce(&diss_V_loc, &diss_V, 1, MPI_DOUBLE, MPI_SUM, 0, comm);
  diss_V *= 2.*nu;
  
  diss_B_loc *= 0.5/double(N*N*N);
  diss_B_loc *= 1. /double(N*N*N);
  MPI_Reduce(&diss_B_loc, &diss_B, 1, MPI_DOUBLE, MPI_SUM, 0, comm);
  diss_B *= 2.*eta;
  
  v0 = sqrt(2./3.*energy_V);
  b0 = sqrt(2./3.*energy_B);
  
  L_V_loc *= 0.5 /double(N*N*N);
  L_V_loc *= 1. /double(N*N*N);
  MPI_Reduce(&L_V_loc, &L_V, 1, MPI_DOUBLE, MPI_SUM, 0, comm);
  //~ //L_V *= M_PI/(2.*v0*v0);
  L_V /= energy_V;
  
  L_B_loc *= 0.5 /double(N*N*N);
  L_B_loc *= 1. /double(N*N*N);
  MPI_Reduce(&L_B_loc, &L_B, 1, MPI_DOUBLE, MPI_SUM, 0, comm);
  //~ //L_B *= M_PI/(2.*b0*b0);
  L_V /= energy_V;
  
  lambda_V = sqrt(15.*nu /diss_V)*v0;
  lambda_B = sqrt(15.*eta/diss_B)*b0;
  
  Re_V = v0*lambda_V/nu ;
  Re_B = b0*lambda_B/eta;
  
  T_V = L_V/v0;
  T_B = L_B/b0;
  
  tau_V = sqrt(nu /diss_V);
  tau_B = sqrt(eta/diss_B);
  
  eta_V = pow(nu *nu *nu / diss_V, 0.25);
  eta_B = pow(eta*eta*eta/ diss_B, 0.25);
  
  // print
  if(myRank==0)
  {
    std::string file_name  = out_dir + "/scales/scales_" + std::to_string(print_count) + ".csv";
    std::ofstream os;
    
    os.open(file_name.c_str(), std::ios::out);
    if(!os){
      std::cout << "Cannot write header to file '" << file_name << "'!\n";
    }
    
    os << std::scientific;
    os << "Scale                  , Velocity     , magnetic Field , Simulation Sizes , , resolved?" << std::endl;
    os << "Energy                 , " << energy_V << " , " << energy_B << std::endl;
    os << "Dissipation            , " << diss_V   << " , " << diss_B   << std::endl;
    os << "RMS                    , " << v0       << " , " << b0       << std::endl;
    os << "micro scale Re         , " << Re_V     << " , " << Re_B     << std::endl;
    os << std::endl;
    os << "Integral Scale         , " << L_V      << " , " << L_B      << "  , L    , " << L            << " , " << (L>L_V && L>L_B)                                   << std::endl;
    os << "LE turnover-time       , " << T_V      << " , " << T_B      << "  , t_out, " << out_interval << " , " << ((0.5*T_V)>out_interval && (0.5*T_B)>out_interval) << std::endl;
    os << "Taylor micro scale     , " << lambda_V << " , " << lambda_B << "  , dx   , " << dx           << " , " << (lambda_V>dx && lambda_B>dx)                       << std::endl;
    os << "Kolmogorov length scale, " << eta_V    << " , " << eta_B    << "  , dx   , " << dx           << " , " << (   eta_V>dx && eta_B   >dx)                       << std::endl;
    os << "Kolmogorov time scale  , "   << tau_V  << " , " << tau_B    << "  , dt   , " << dt           << " , " << (tau_V>dt && tau_B>dt)                             << std::endl;
    
    os.close();
  }MPI_Barrier(comm);
  
}

void CSpecDyn::print()
{
  print_vti();
  print_scales();
  print_EnergySpectrum();
  
  if(myRank==0)
  {
    printf("Printing # %d!\n", print_count);
  }MPI_Barrier(comm);
  
  print_count++;
}

void CSpecDyn::read_binary()
{
  
  std::string binary_file;
  //std::string field_names[] = {"V0", "V1", "V2", "B0", "B1", "B2"};
  std::string field_names[] = {"W0", "W1", "W2", "Z0", "Z1", "Z2"}; // Elsässer Coords
  double* field_ptrs[] = {Vx_R, Vy_R, Vz_R, Bx_R, By_R, Bz_R};
  double* field = NULL;
  
  int size_total[3] = {N,N,N};
  MPI_Datatype subarray;
  MPI_Type_create_subarray(3, size_total, size_R, start_R, MPI_ORDER_C, MPI_FLOAT, &subarray);
  MPI_Type_commit(&subarray);
  MPI_File mpi_file;
  
  float* buffer = new float[size_R_tot];
  
  /****************************/
  // iterate over all 6 field components
  for(int i = 0; i < 6; i++)
  {
  
    // open MPI file
    binary_file = BINARY_DIR + "/" + field_names[i] + ".dat";
    field = field_ptrs[i];
    
    MPI_File_open(comm, binary_file.c_str(), MPI_MODE_RDONLY, MPI_INFO_NULL, &mpi_file);
    
    // set specific file view for each process
    MPI_Offset file_offset = 0;
    MPI_File_set_view(mpi_file, file_offset, MPI_FLOAT, subarray, "native", MPI_INFO_NULL);
    MPI_File_read_all(mpi_file, buffer, size_R_tot, MPI_FLOAT, MPI_STATUS_IGNORE);
    
    for(int id = 0; id < size_R_tot; id++)
    {
      field[id] = buffer[id];
    }
  }
    
  /****************************/
  
  // Elsässer -> V,B
  for(int id = 0; id < size_R_tot; id++)
  {
    double Z_x = Vx_R[id];
    double Z_y = Vy_R[id];
    double Z_z = Vz_R[id];
    double W_x = Bx_R[id];
    double W_y = By_R[id];
    double W_z = Bz_R[id];
    
    Vx_R[id] = 0.5 * (Z_x + W_x);
    Vy_R[id] = 0.5 * (Z_y + W_y);
    Vz_R[id] = 0.5 * (Z_z + W_z);
    Bx_R[id] = 0.5 * (Z_x - W_x);
    By_R[id] = 0.5 * (Z_y - W_y);
    Bz_R[id] = 0.5 * (Z_z - W_z);
  }
  
  
  fFFT(Vx_R, Vy_R, Vz_R, Vx_F, Vy_F, Vz_F);
  fFFT(Bx_R, By_R, Bz_R, Bx_F, By_F, Bz_F);
  
  dealias(Vx_F, Vy_F, Vz_F);
  dealias(Bx_F, By_F, Bz_F);
  
  projection(Vx_F, Vy_F, Vz_F);
  projection(Bx_F, By_F, Bz_F);
  
  delete[] buffer;
  
  // close MPI file
  MPI_File_close(&mpi_file);
  MPI_Barrier(MPI_COMM_WORLD);
}

void CSpecDyn::restart()
{
  
  #ifdef RESTART_DIR
  std::string restart_dir = RESTART_DIR;
  #else
  std::string restart_dir = out_dir;
  #endif
  
  std::string restart_file = restart_dir + "/vti/step_" + std::to_string(RESTART_STEP) + ".vti";
  
  if(myRank==0)
  {
    printf("Restart!\n");
  }
  
  // CHECK RESOLUTION N
  if(myRank==0)
  {
    int N_in = 0;
    
    std::ifstream reader;
    reader.open(restart_file, std::ios::in);
    if(!reader){
      std::cout << "Cannot read header to file '" << restart_file << "'!\n";
    }
    
    char word[] = "012345";
    
    for(int i = 0; i < 200; i++)
    {
      reader.seekg(i, std::ios::beg);
      reader.read(word, 6);

      if(strcmp(word,"Extent")==0)
      {
        reader.seekg(i+10, std::ios::beg); // set to start of N
        
        char N_word[10];
        char c; 
        reader.read(&c, sizeof(char));
        for(int i = 0; ; i++,reader.read(&c, sizeof(char)))
        {
          
          if(c==' ')
          {
            N_word[i] = '\0';
            break;
          }
          else{
            N_word[i] = c;
          }
        }
        
        N_in = std::stoi(N_word)+1;
        
        break;
      }
    }
    reader.close();
    
    if(N==N_in)
    {
      printf("Restart resolution matches!\n");
    }
    else
    {
      printf("Resolution N of vti file for restart does not match parameters!\n");
      exit(EXIT_FAILURE);
    }
  }
  
  // GET STEP OUTPUT #
  int step = 0;
  
  if(myRank==0)
  {
    char number_str[10];
    int str_counter = 0;
    
    for(int i = restart_file.size()-5; restart_file[i]!='_'; i--)
    {
      number_str[str_counter] = restart_file[i];
      str_counter ++;
    }
    number_str[str_counter] = '\0';
    
    std::string number_rev = std::string(number_str);
    std::string number = std::string(number_rev.length(), 'x');
    
    for(int i=0, n=number.length()-1; i<number.length(); i++, n--)
    {
      number[i] = number_rev[n];
    }
    
    print_count = std::stoi(number)+1;
    
  }MPI_Barrier(comm);
  
  MPI_Bcast(&print_count, 1, MPI_INT, 0, comm);
  
  // CHECK TIME
  if(myRank==0)
  {
    std::ifstream reader;
    reader.open(restart_file, std::ios::in);
    if(!reader){
      std::cout << "Cannot read header to file '" << restart_file << "'!\n";
    }
    
    char word[] = "TimeValue";
    
    for(int i = 0; i < 1000; i++)
    {
      reader.seekg(i, std::ios::beg);
      reader.read(word, 9);

      if(strcmp(word,"TimeValue")==0)
      {
        reader.seekg(i+54, std::ios::beg); // set to start of time
        
        char t_word[20];
        char c; 
        reader.read(&c, sizeof(char));
        for(int i = 0; ; i++,reader.read(&c, sizeof(char)))
        {
          
          if(c==' ')
          {
            t_word[i] = '\0';
            break;
          }
          else{
            t_word[i] = c;
          }
        }
        
        //~ //printf("new time: %s\n", t_word);
        time = std::stod(t_word);
        
        break;
      }
    }
    reader.close();
    
  }MPI_Barrier(comm);
  
  MPI_Bcast(&time, 1, MPI_DOUBLE, 0, comm);
  
  // GET FIELDS
  // find header offset
  char character;
  int headerOffset = 0;
  
  if(myRank==0)
  {
    // find length of header
    std::ifstream reader;
    reader.open(restart_file, std::ios::in);
    if(!reader){
      std::cout << "Cannot read header to file '" << restart_file << "'!\n";
    }
    
    for(int i = 0; i < 1e5; i++)
    {
      reader.read(&character, sizeof(char));
      headerOffset += 1;
      if(character == ' ')
      {
        reader.read(&character, sizeof(char));
        
        if(character == '_')
        {
          headerOffset += 1;
          break;
        }
        else
        {
          reader.seekg(headerOffset, std::ios::beg);
        }
          
      }
    }
    
    reader.close();
    
  }MPI_Barrier(MPI_COMM_WORLD);
  
  MPI_Bcast(&headerOffset, 1, MPI_INT, 0, comm);
  
  // open MPI file
  MPI_File mpi_file;
  MPI_File_open(comm, restart_file.c_str(), MPI_MODE_RDONLY, MPI_INFO_NULL, &mpi_file);
  
  float* buffer = new float[3*size_R_tot];
  
  // set specific file view for each process
  MPI_Offset file_offset = headerOffset + sizeof(uint64_t);
  MPI_File_set_view(mpi_file, file_offset, vti_float3, vti_subarray_vector, "native", MPI_INFO_NULL);
  MPI_File_read_all(mpi_file, buffer, size_R_tot, vti_float3, MPI_STATUS_IGNORE);
  
  for(int id = 0; id < size_R_tot; id++)
  {
    Vx_R[id] = buffer[id*3  ];
    Vy_R[id] = buffer[id*3+1];
    Vz_R[id] = buffer[id*3+2];
  }
  
  long N_l = N;
  file_offset = headerOffset + 2*sizeof(uint64_t) + 3*N_l*N_l*N_l*sizeof(float);
  MPI_File_set_view(mpi_file, 0, MPI_CHAR, MPI_CHAR, "native", MPI_INFO_NULL); // reset file view to specify offset in bytes in the next file view
  MPI_File_set_view(mpi_file, file_offset, vti_float3, vti_subarray_vector, "native", MPI_INFO_NULL);
  MPI_File_read_all(mpi_file, buffer, size_R_tot, vti_float3, MPI_STATUS_IGNORE);
  
  for(int id = 0; id < size_R_tot; id++)
  {
    Bx_R[id] = buffer[id*3  ];
    By_R[id] = buffer[id*3+1];
    Bz_R[id] = buffer[id*3+2];
  }
  
  fFFT(Vx_R, Vy_R, Vz_R, Vx_F, Vy_F, Vz_F);
  fFFT(Bx_R, By_R, Bz_R, Bx_F, By_F, Bz_F);
  
  delete[] buffer;
  
  // close MPI file
  MPI_File_close(&mpi_file);
  MPI_Barrier(MPI_COMM_WORLD);
  
}
