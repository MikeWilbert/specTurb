#include "MikeFFT.h"

plan3D_r2c::plan3D_r2c(fftw_plan plan_c2c_x_, fftw_plan plan_c2c_y_, fftw_plan plan_r2c_z_) :
plan_c2c_x(plan_c2c_x_), plan_c2c_y(plan_c2c_y_), plan_r2c_z(plan_r2c_z_)
{
  // nothing
}

plan3D_c2r::plan3D_c2r(fftw_plan plan_c2c_x_, fftw_plan plan_c2c_y_, fftw_plan plan_c2r_z_) :
plan_c2c_x(plan_c2c_x_), plan_c2c_y(plan_c2c_y_), plan_c2r_z(plan_c2r_z_)
{
  // nothing
}

plan3D_c2c::plan3D_c2c(fftw_plan plan_c2c_x_, fftw_plan plan_c2c_y_, fftw_plan plan_c2c_z_) :
plan_c2c_x(plan_c2c_x_), plan_c2c_y(plan_c2c_y_), plan_c2c_z(plan_c2c_z_)
{
  // nothing
}

MikeFFT::MikeFFT(int N_, const int* pdims_) :
N(N_), pdims(pdims_)
{
  // MPI
  int periods[] = {0,0};
  MPI_Cart_create(MPI_COMM_WORLD, 2, pdims, periods, 1, &comm);
  
  MPI_Comm_size(comm, &nprocs);
	MPI_Comm_rank(comm, &myRank);
  MPI_Cart_coords(comm, myRank, 2, cart_coords);
  
  MPI_Type_contiguous(2, MPI_DOUBLE, &mpiComplex);
  MPI_Type_commit(&mpiComplex);
  
  int remain_dims[2];
  
  remain_dims[0] = 1;
  remain_dims[1] = 0;
  MPI_Cart_sub(comm, remain_dims, &comm_xy);
  
  remain_dims[0] = 0;
  remain_dims[1] = 1;
  MPI_Cart_sub(comm, remain_dims, &comm_yz);
  
  // get sizes
  // position space: z-continuous
  // Fourier  space: x continuous
  split(N    ,              0,        1, size_Fx[0], start_Fx[0]);
  split(N    , cart_coords[0], pdims[0], size_Fx[1], start_Fx[1]);
  split(N/2+1, cart_coords[1], pdims[1], size_Fx[2], start_Fx[2]);
  
  split(N    , cart_coords[0], pdims[0], size_Fy[0], start_Fy[0]);
  split(N    ,              0,        1, size_Fy[1], start_Fy[1]);
  split(N/2+1, cart_coords[1], pdims[1], size_Fy[2], start_Fy[2]);
  
  split(N    , cart_coords[0], pdims[0], size_Fz[0], start_Fz[0]);
  split(N    , cart_coords[1], pdims[1], size_Fz[1], start_Fz[1]);
  split(N/2+1,             0,         1, size_Fz[2], start_Fz[2]);
  
  split(N    , cart_coords[0], pdims[0], size_Rz[0], start_Rz[0]);
  split(N    , cart_coords[1], pdims[1], size_Rz[1], start_Rz[1]);
  split(N    ,             0,         1, size_Rz[2], start_Rz[2]);
  
  split(N    ,              0,        1, size_Cx[0], start_Cx[0]);
  split(N    , cart_coords[0], pdims[0], size_Cx[1], start_Cx[1]);
  split(N    , cart_coords[1], pdims[1], size_Cx[2], start_Cx[2]);
  
  split(N    , cart_coords[0], pdims[0], size_Cy[0], start_Cy[0]);
  split(N    ,              0,        1, size_Cy[1], start_Cy[1]);
  split(N    , cart_coords[1], pdims[1], size_Cy[2], start_Cy[2]);
  
  split(N    , cart_coords[0], pdims[0], size_Cz[0], start_Cz[0]);
  split(N    , cart_coords[1], pdims[1], size_Cz[1], start_Cz[1]);
  split(N    ,             0,         1, size_Cz[2], start_Cz[2]);
  
  size_Fx_tot = size_Fx[0]*size_Fx[1]*size_Fx[2];
  size_Fy_tot = size_Fy[0]*size_Fy[1]*size_Fy[2];
  size_Fz_tot = size_Fz[0]*size_Fz[1]*size_Fz[2];
  size_Rz_tot = size_Rz[0]*size_Rz[1]*size_Rz[2];
  size_Cx_tot = size_Cx[0]*size_Cx[1]*size_Cx[2];
  size_Cy_tot = size_Cy[0]*size_Cy[1]*size_Cy[2];
  size_Cz_tot = size_Cz[0]*size_Cz[1]*size_Cz[2];
  
  // init subarrays
  subarrays_F_X  = new MPI_Datatype[pdims[0]];
  subarrays_F_Yx = new MPI_Datatype[pdims[0]];
  subarrays_F_Yz = new MPI_Datatype[pdims[1]];
  subarrays_F_Z  = new MPI_Datatype[pdims[1]];
  subarrays_C_X  = new MPI_Datatype[pdims[0]];
  subarrays_C_Yx = new MPI_Datatype[pdims[0]];
  subarrays_C_Yz = new MPI_Datatype[pdims[1]];
  subarrays_C_Z  = new MPI_Datatype[pdims[1]];
  
  setup_subarrays();
  
  send_count_xy = new int[pdims[0]];
  send_disp_xy  = new int[pdims[0]];
  send_count_yz = new int[pdims[1]];
  send_disp_yz  = new int[pdims[1]];
  
  for(int i = 0; i < pdims[0]; i++)
  {
    send_count_xy[i] = 1;
    send_disp_xy[i]  = 0;
  }
  
  for(int i = 0; i < pdims[1]; i++)
  {
    send_count_yz[i] = 1;
    send_disp_yz[i]  = 0;
  }
  
  // malloc fields
  field_X  = reinterpret_cast<CX*>(fftw_alloc_complex((size_t)size_Cx_tot));
  field_Yx = reinterpret_cast<CX*>(fftw_alloc_complex((size_t)size_Cy_tot));
  field_Yz = reinterpret_cast<CX*>(fftw_alloc_complex((size_t)size_Cy_tot));
  field_Z  = reinterpret_cast<CX*>(fftw_alloc_complex((size_t)size_Cz_tot));
  
  // plan FFT
  plan_r2c = make_plan3D_r2c();
  plan_c2r = make_plan3D_c2r();
  plan_c2c_fwrd = make_plan3D_c2c_fwrd();
  plan_c2c_bwrd = make_plan3D_c2c_bwrd();
  
}

void MikeFFT::split(int N, int rank, int split_size, int& size, int& start)
{
	int q = N / split_size;
	int r = N % split_size;
	
	size  = q + (rank < r);
	start = q * rank + std::min(r, rank);
}

void MikeFFT::setup_subarrays()
{
  int size_tmp1[3];
  int size_tmp2[3];
  int start_tmp1[3];
  int start_tmp2[3];
  
  // REAL position space
  
  // x<->y
  for(int i = 0; i < pdims[0]; i++)
  {
    size_tmp1[0] = size_Fx[0];
    size_tmp1[1] = size_Fx[1];
    size_tmp1[2] = size_Fx[2];
    
    start_tmp1[0] = 0; // starts in partitioned local array!
    start_tmp1[1] = 0;
    start_tmp1[2] = 0;
    
    size_tmp2[0] = size_Fy[0];
    size_tmp2[1] = size_Fy[1];
    size_tmp2[2] = size_Fy[2];
    
    start_tmp2[0] = 0;
    start_tmp2[1] = 0;
    start_tmp2[2] = 0;
    
    split(size_Fx[0], i, pdims[0], size_tmp1[0], start_tmp1[0]); // split continuous dimension
    split(size_Fy[1], i, pdims[0], size_tmp2[1], start_tmp2[1]);
                                                         
    MPI_Type_create_subarray(3, size_Fx, size_tmp1, start_tmp1, MPI_ORDER_C, mpiComplex, &subarrays_F_X [i]);
    MPI_Type_create_subarray(3, size_Fy, size_tmp2, start_tmp2, MPI_ORDER_C, mpiComplex, &subarrays_F_Yx[i]);
    
    MPI_Type_commit(&subarrays_F_X [i]);
    MPI_Type_commit(&subarrays_F_Yx[i]);
  
  }
  
  // y<->z
  for(int i = 0; i < pdims[1]; i++)
  {
    size_tmp1[0] = size_Fy[0];
    size_tmp1[1] = size_Fy[1];
    size_tmp1[2] = size_Fy[2];
    
    start_tmp1[0] = 0; // starts in partitioned local array!
    start_tmp1[1] = 0;
    start_tmp1[2] = 0;
    
    size_tmp2[0] = size_Fz[0];
    size_tmp2[1] = size_Fz[1];
    size_tmp2[2] = size_Fz[2];
    
    start_tmp2[0] = 0;
    start_tmp2[1] = 0;
    start_tmp2[2] = 0;
    
    split(size_Fy[1], i, pdims[1], size_tmp1[1], start_tmp1[1]);
    split(size_Fz[2], i, pdims[1], size_tmp2[2], start_tmp2[2]);
                                                         
    MPI_Type_create_subarray(3, size_Fy, size_tmp1, start_tmp1, MPI_ORDER_C, mpiComplex, &subarrays_F_Yz[i]);
    MPI_Type_create_subarray(3, size_Fz, size_tmp2, start_tmp2, MPI_ORDER_C, mpiComplex, &subarrays_F_Z [i]);
    
    MPI_Type_commit(&subarrays_F_Yz[i]);
    MPI_Type_commit(&subarrays_F_Z [i]);
  
  }
  
  // COMPLEX position space
  
  // x<->y
  for(int i = 0; i < pdims[0]; i++)
  {
    size_tmp1[0] = size_Cx[0];
    size_tmp1[1] = size_Cx[1];
    size_tmp1[2] = size_Cx[2];
    
    start_tmp1[0] = 0; // starts in partitioned local array!
    start_tmp1[1] = 0;
    start_tmp1[2] = 0;
    
    size_tmp2[0] = size_Cy[0];
    size_tmp2[1] = size_Cy[1];
    size_tmp2[2] = size_Cy[2];
    
    start_tmp2[0] = 0;
    start_tmp2[1] = 0;
    start_tmp2[2] = 0;
    
    split(size_Cx[0], i, pdims[0], size_tmp1[0], start_tmp1[0]); // split continuous dimension
    split(size_Cy[1], i, pdims[0], size_tmp2[1], start_tmp2[1]);
                                                         
    MPI_Type_create_subarray(3, size_Cx, size_tmp1, start_tmp1, MPI_ORDER_C, mpiComplex, &subarrays_C_X [i]);
    MPI_Type_create_subarray(3, size_Cy, size_tmp2, start_tmp2, MPI_ORDER_C, mpiComplex, &subarrays_C_Yx[i]);
    
    MPI_Type_commit(&subarrays_C_X [i]);
    MPI_Type_commit(&subarrays_C_Yx[i]);
  
  }
  
   // y<->z
  for(int i = 0; i < pdims[1]; i++)
  {
    size_tmp1[0] = size_Cy[0];
    size_tmp1[1] = size_Cy[1];
    size_tmp1[2] = size_Cy[2];
    
    start_tmp1[0] = 0; // starts in partitioned local array!
    start_tmp1[1] = 0;
    start_tmp1[2] = 0;
    
    size_tmp2[0] = size_Cz[0];
    size_tmp2[1] = size_Cz[1];
    size_tmp2[2] = size_Cz[2];
    
    start_tmp2[0] = 0;
    start_tmp2[1] = 0;
    start_tmp2[2] = 0;
    
    split(size_Cy[1], i, pdims[1], size_tmp1[1], start_tmp1[1]);
    split(size_Cz[2], i, pdims[1], size_tmp2[2], start_tmp2[2]);
                                                         
    MPI_Type_create_subarray(3, size_Cy, size_tmp1, start_tmp1, MPI_ORDER_C, mpiComplex, &subarrays_C_Yz[i]);
    MPI_Type_create_subarray(3, size_Cz, size_tmp2, start_tmp2, MPI_ORDER_C, mpiComplex, &subarrays_C_Z [i]);
    
    MPI_Type_commit(&subarrays_C_Yz[i]);
    MPI_Type_commit(&subarrays_C_Z [i]);
  
  }
  
}

void MikeFFT::get_sizes_real(int size_R_[3], int start_R_[3], int size_F_[3], int start_F_[3])
{
  size_R_[0] = size_Rz[0];
  size_R_[1] = size_Rz[1];
  size_R_[2] = size_Rz[2];
  
  size_F_[0] = size_Fx[0];
  size_F_[1] = size_Fx[1];
  size_F_[2] = size_Fx[2];
  
  start_R_[0] = start_Rz[0];
  start_R_[1] = start_Rz[1];
  start_R_[2] = start_Rz[2];
  
  start_F_[0] = start_Fx[0];
  start_F_[1] = start_Fx[1];
  start_F_[2] = start_Fx[2];
}

void MikeFFT::get_sizes_cplx(int size_R_[3], int start_R_[3], int size_F_[3], int start_F_[3])
{
  size_R_[0] = size_Cz[0];
  size_R_[1] = size_Cz[1];
  size_R_[2] = size_Cz[2];
  
  size_F_[0] = size_Cx[0];
  size_F_[1] = size_Cx[1];
  size_F_[2] = size_Cx[2];
  
  start_R_[0] = start_Cz[0];
  start_R_[1] = start_Cz[1];
  start_R_[2] = start_Cz[2];
  
  start_F_[0] = start_Cx[0];
  start_F_[1] = start_Cx[1];
  start_F_[2] = start_Cx[2];
}

void MikeFFT::get_comm(MPI_Comm& comm_)
{
  comm_ = comm;
}

double*               MikeFFT::malloc_R()
{
  return fftw_alloc_real((size_t)size_Rz_tot);
}

std::complex<double>* MikeFFT::malloc_F()
{
  return reinterpret_cast<CX*>(fftw_alloc_complex((size_t)size_Fx_tot));
}

std::complex<double>* MikeFFT::malloc_R_cplx()
{
  return reinterpret_cast<CX*>(fftw_alloc_complex((size_t)size_Cz_tot));
}

std::complex<double>* MikeFFT::malloc_F_cplx()
{
  return reinterpret_cast<CX*>(fftw_alloc_complex((size_t)size_Cx_tot));
}

void MikeFFT::free_R(double* A)
{
  fftw_free(reinterpret_cast<void*>(A));
}

void MikeFFT::free_C(std::complex<double>* A)
{
  fftw_free(reinterpret_cast<void*>(A));
}

plan3D_r2c MikeFFT::make_plan3D_r2c()
{
  double* in_R  = fftw_alloc_real((size_t)size_Rz_tot);
  fftw_complex* out_F  = fftw_alloc_complex((size_t)size_Fx_tot);
  
  fftw_plan plan_x;
  fftw_plan plan_y;
  fftw_plan plan_z;
  
  fftw_iodim dims[1];
  fftw_iodim howmany_dims[2];
  int stride[3];
  int stride_R[3];

  // Z
  stride_R[0] = size_Rz[1]*size_Rz[2];
  stride_R[1] = size_Rz[2];
  stride_R[2] = 1;
  stride  [0] = size_Fz[1]*size_Fz[2];
  stride  [1] = size_Fz[2];
  stride  [2] = 1;
  
  dims[0].n  = size_Rz [2]; // dims.n needs to be the size of the real array
  dims[0].is = stride_R[2];
  dims[0].os = stride  [2];
  
  howmany_dims[0].n  = size_Rz [0];
  howmany_dims[0].is = stride_R[0]; 
  howmany_dims[0].os = stride  [0]; 
  
  howmany_dims[1].n  = size_Rz [1];
  howmany_dims[1].is = stride_R[1]; 
  howmany_dims[1].os = stride  [1]; 
  
  plan_z = fftw_plan_guru_dft_r2c(1, dims, 2, howmany_dims, 
           in_R, reinterpret_cast<fftw_complex*>(field_Z), FFTW_MEASURE);
  
  // Y
  stride[0] = size_Fy[1]*size_Fy[2];
  stride[1] = size_Fy[2];
  stride[2] = 1;
  
  dims[0].n  = size_Fy[1];
  dims[0].is = stride [1];
  dims[0].os = stride [1];
  
  howmany_dims[0].n  = size_Fy[0];
  howmany_dims[0].is = stride [0]; 
  howmany_dims[0].os = stride [0]; 
  
  howmany_dims[1].n  = size_Fy[2];
  howmany_dims[1].is = stride [2]; 
  howmany_dims[1].os = stride [2]; 
  
  plan_y = fftw_plan_guru_dft(1, dims, 2, howmany_dims, 
          reinterpret_cast<fftw_complex*>(field_Yz), reinterpret_cast<fftw_complex*>(field_Yx), FFTW_FORWARD, FFTW_MEASURE);
  
  // X
  stride[0] = size_Fx[1]*size_Fx[2];
  stride[1] = size_Fx[2];
  stride[2] = 1;
  
  dims[0].n  = size_Fx[0];
  dims[0].is = stride [0];
  dims[0].os = stride [0];
  
  howmany_dims[0].n  = size_Fx[1];
  howmany_dims[0].is = stride [1]; 
  howmany_dims[0].os = stride [1]; 
  
  howmany_dims[1].n  = size_Fx[2];
  howmany_dims[1].is = stride [2]; 
  howmany_dims[1].os = stride [2]; 
           
  plan_x = fftw_plan_guru_dft(1, dims, 2, howmany_dims, 
           reinterpret_cast<fftw_complex*>(field_X), reinterpret_cast<fftw_complex*>(out_F), FFTW_FORWARD, FFTW_MEASURE);
  
  fftw_free(in_R);
  fftw_free(out_F);
  
  return plan3D_r2c(plan_x, plan_y, plan_z);
}

plan3D_c2c MikeFFT::make_plan3D_c2c_fwrd()
{
  fftw_complex* in_C   = fftw_alloc_complex((size_t)size_Cz_tot);
  fftw_complex* out_C  = fftw_alloc_complex((size_t)size_Cx_tot);
  
  fftw_plan plan_x;
  fftw_plan plan_y;
  fftw_plan plan_z;
  
  fftw_iodim dims[1];
  fftw_iodim howmany_dims[2];
  int stride[3];
  int stride_R[3];

  // Z
  stride[0] = size_Cz[1]*size_Cz[2];
  stride[1] = size_Cz[2];
  stride[2] = 1;
  
  dims[0].n  = size_Cz[2];
  dims[0].is = stride [2];
  dims[0].os = stride [2];
  
  howmany_dims[0].n  = size_Cz[0];// dims.n needs to be the size of the real array [OLD COMMENT CROM R2C VERSION]
  howmany_dims[0].is = stride [0]; 
  howmany_dims[0].os = stride [0]; 
  
  howmany_dims[1].n  = size_Cz[1];
  howmany_dims[1].is = stride [1]; 
  howmany_dims[1].os = stride [1]; 
  
  plan_z = fftw_plan_guru_dft(1, dims, 2, howmany_dims, 
         reinterpret_cast<fftw_complex*>(in_C), reinterpret_cast<fftw_complex*>(field_Z), FFTW_FORWARD, FFTW_MEASURE);
  
  // Y
  stride[0] = size_Cy[1]*size_Cy[2];
  stride[1] = size_Cy[2];
  stride[2] = 1;
  
  dims[0].n  = size_Cy[1];
  dims[0].is = stride [1];
  dims[0].os = stride [1];
  
  howmany_dims[0].n  = size_Cy[0];
  howmany_dims[0].is = stride [0]; 
  howmany_dims[0].os = stride [0]; 
  
  howmany_dims[1].n  = size_Cy[2];
  howmany_dims[1].is = stride [2]; 
  howmany_dims[1].os = stride [2]; 
  
  plan_y = fftw_plan_guru_dft(1, dims, 2, howmany_dims, 
          reinterpret_cast<fftw_complex*>(field_Yz), reinterpret_cast<fftw_complex*>(field_Yx), FFTW_FORWARD, FFTW_MEASURE);
  
  // X
  stride[0] = size_Cx[1]*size_Cx[2];
  stride[1] = size_Cx[2];
  stride[2] = 1;
  
  dims[0].n  = size_Cx[0];
  dims[0].is = stride [0];
  dims[0].os = stride [0];
  
  howmany_dims[0].n  = size_Cx[1];
  howmany_dims[0].is = stride [1]; 
  howmany_dims[0].os = stride [1]; 
  
  howmany_dims[1].n  = size_Cx[2];
  howmany_dims[1].is = stride [2]; 
  howmany_dims[1].os = stride [2]; 
           
  plan_x = fftw_plan_guru_dft(1, dims, 2, howmany_dims, 
           reinterpret_cast<fftw_complex*>(field_X), reinterpret_cast<fftw_complex*>(out_C), FFTW_FORWARD, FFTW_MEASURE);
  
  fftw_free(in_C);
  fftw_free(out_C);
  
  return plan3D_c2c(plan_x, plan_y, plan_z);
}

plan3D_c2r MikeFFT::make_plan3D_c2r()
{
  fftw_complex* in_F  = fftw_alloc_complex((size_t)size_Fx_tot);
  double* out_R  = fftw_alloc_real((size_t)size_Rz_tot);
  
  fftw_plan plan_x;
  fftw_plan plan_y;
  fftw_plan plan_z;
  
  fftw_iodim dims[1];
  fftw_iodim howmany_dims[2];
  int stride[3];
  int stride_R[3];
  
  // X
  stride[0] = size_Fx[1]*size_Fx[2];
  stride[1] = size_Fx[2];
  stride[2] = 1;
  
  dims[0].n  = size_Fx[0];
  dims[0].is = stride [0];
  dims[0].os = stride [0];
  
  howmany_dims[0].n  = size_Fx[1];
  howmany_dims[0].is = stride [1]; 
  howmany_dims[0].os = stride [1]; 
  
  howmany_dims[1].n  = size_Fx[2];
  howmany_dims[1].is = stride [2]; 
  howmany_dims[1].os = stride [2]; 
  
  plan_x = fftw_plan_guru_dft(1, dims, 2, howmany_dims, 
           reinterpret_cast<fftw_complex*>(in_F), reinterpret_cast<fftw_complex*>(field_X), FFTW_BACKWARD, FFTW_MEASURE);
  
  // Y
  stride[0] = size_Fy[1]*size_Fy[2];
  stride[1] = size_Fy[2];
  stride[2] = 1;
  
  dims[0].n  = size_Fy[1];
  dims[0].is = stride [1];
  dims[0].os = stride [1];
  
  howmany_dims[0].n  = size_Fy[0];
  howmany_dims[0].is = stride [0]; 
  howmany_dims[0].os = stride [0]; 
  
  howmany_dims[1].n  = size_Fy[2];
  howmany_dims[1].is = stride [2]; 
  howmany_dims[1].os = stride [2]; 
  
  plan_y = fftw_plan_guru_dft(1, dims, 2, howmany_dims, 
          reinterpret_cast<fftw_complex*>(field_Yx), reinterpret_cast<fftw_complex*>(field_Yz), FFTW_BACKWARD, FFTW_MEASURE);
  
  // Z
  stride_R[0] = size_Rz[1]*size_Rz[2];
  stride_R[1] = size_Rz[2];
  stride_R[2] = 1;
  stride  [0] = size_Fz[1]*size_Fz[2];
  stride  [1] = size_Fz[2];
  stride  [2] = 1;
  
  dims[0].n  = size_Rz [2]; // dims.n needs to be the size of the real array
  dims[0].is = stride  [2];
  dims[0].os = stride_R[2];
  
  howmany_dims[0].n  = size_Rz [0]; 
  howmany_dims[0].is = stride  [0]; 
  howmany_dims[0].os = stride_R[0]; 
  
  howmany_dims[1].n  = size_Rz [1];
  howmany_dims[1].is = stride  [1]; 
  howmany_dims[1].os = stride_R[1]; 
  
  plan_z = fftw_plan_guru_dft_c2r(1, dims, 2, howmany_dims, 
         reinterpret_cast<fftw_complex*>(field_Z), out_R, FFTW_MEASURE);
           
  fftw_free(in_F);
  fftw_free(out_R);
           
  return plan3D_c2r(plan_x, plan_y, plan_z);
  
}

plan3D_c2c MikeFFT::make_plan3D_c2c_bwrd()
{
  fftw_complex* in_C   = fftw_alloc_complex((size_t)size_Cx_tot);
  fftw_complex* out_C  = fftw_alloc_complex((size_t)size_Cz_tot);
  
  fftw_plan plan_x;
  fftw_plan plan_y;
  fftw_plan plan_z;
  
  fftw_iodim dims[1];
  fftw_iodim howmany_dims[2];
  int stride[3];
  int stride_R[3];

  // X
  stride[0] = size_Cx[1]*size_Cx[2];
  stride[1] = size_Cx[2];
  stride[2] = 1;
  
  dims[0].n  = size_Cx[0];
  dims[0].is = stride [0];
  dims[0].os = stride [0];
  
  howmany_dims[0].n  = size_Cx[1];
  howmany_dims[0].is = stride [1]; 
  howmany_dims[0].os = stride [1]; 
  
  howmany_dims[1].n  = size_Cx[2];
  howmany_dims[1].is = stride [2]; 
  howmany_dims[1].os = stride [2]; 
  
  plan_x = fftw_plan_guru_dft(1, dims, 2, howmany_dims, 
           reinterpret_cast<fftw_complex*>(in_C), reinterpret_cast<fftw_complex*>(field_X), FFTW_BACKWARD, FFTW_MEASURE);
  
  // Y
  stride[0] = size_Cy[1]*size_Cy[2];
  stride[1] = size_Cy[2];
  stride[2] = 1;
  
  dims[0].n  = size_Cy[1];
  dims[0].is = stride [1];
  dims[0].os = stride [1];
  
  howmany_dims[0].n  = size_Cy[0];
  howmany_dims[0].is = stride [0]; 
  howmany_dims[0].os = stride [0]; 
  
  howmany_dims[1].n  = size_Cy[2];
  howmany_dims[1].is = stride [2]; 
  howmany_dims[1].os = stride [2]; 
  
  plan_y = fftw_plan_guru_dft(1, dims, 2, howmany_dims, 
          reinterpret_cast<fftw_complex*>(field_Yx), reinterpret_cast<fftw_complex*>(field_Yz), FFTW_BACKWARD, FFTW_MEASURE);
  
  // Z
  stride[0] = size_Cz[1]*size_Cz[2];
  stride[1] = size_Cz[2];
  stride[2] = 1;
  
  dims[0].n  = size_Cz[2];
  dims[0].is = stride [2];
  dims[0].os = stride [2];
  
  howmany_dims[0].n  = size_Cz[0];
  howmany_dims[0].is = stride [0]; 
  howmany_dims[0].os = stride [0]; 
  
  howmany_dims[1].n  = size_Cz[1];
  howmany_dims[1].is = stride [1]; 
  howmany_dims[1].os = stride [1]; 
  
  plan_z = fftw_plan_guru_dft(1, dims, 2, howmany_dims, 
           reinterpret_cast<fftw_complex*>(field_X), reinterpret_cast<fftw_complex*>(out_C), FFTW_BACKWARD, FFTW_MEASURE);
 
  fftw_free(in_C);
  fftw_free(out_C);
  
  return plan3D_c2c(plan_x, plan_y, plan_z);
}

void MikeFFT::r2c(double* IN, std::complex<double>* OUT)
{
  
  // FFT Z
  fftw_execute_dft_r2c(plan_r2c.plan_r2c_z, IN, reinterpret_cast<fftw_complex*>(field_Z));
  
  // Z -> Y
  MPI_Alltoallw( field_Z , send_count_yz, send_disp_yz, subarrays_F_Z,
                 field_Yz, send_count_yz, send_disp_yz, subarrays_F_Yz, comm_yz);
                 
  // FFT Y
  fftw_execute(plan_r2c.plan_c2c_y);
  
  // Y -> X
  MPI_Alltoallw( field_Yx, send_count_xy, send_disp_xy, subarrays_F_Yx,
                 field_X , send_count_xy, send_disp_xy, subarrays_F_X, comm_xy);
  
  // FFT X
  fftw_execute_dft(plan_r2c.plan_c2c_x, reinterpret_cast<fftw_complex*>(field_X), reinterpret_cast<fftw_complex*>(OUT));
  
}

void MikeFFT::c2r(std::complex<double>* IN, double* OUT)
{
  
  // FFT X
  fftw_execute_dft(plan_c2r.plan_c2c_x, reinterpret_cast<fftw_complex*>(IN), reinterpret_cast<fftw_complex*>(field_X));
  
  // X -> Y
  MPI_Alltoallw( field_X , send_count_xy, send_disp_xy, subarrays_F_X ,
                 field_Yx, send_count_xy, send_disp_xy, subarrays_F_Yx, comm_xy);
                 
  // FFT Y
  fftw_execute(plan_c2r.plan_c2c_y);
  
  // Y -> Z
  MPI_Alltoallw( field_Yz, send_count_yz, send_disp_yz, subarrays_F_Yz,
                 field_Z , send_count_yz, send_disp_yz, subarrays_F_Z , comm_yz);
                 
  // FFT Z
  fftw_execute_dft_c2r(plan_c2r.plan_c2r_z, reinterpret_cast<fftw_complex*>(field_Z), OUT);
  
  // normalize
  double N3_inv = 1./(N*N*N);
  for(int id = 0; id < size_Rz_tot;id++)
  {
    OUT[id] *= N3_inv;
  }

}

void MikeFFT::c2c_fwrd(std::complex<double>* IN, std::complex<double>* OUT)
{
  
  // FFT Z
  fftw_execute_dft(plan_c2c_fwrd.plan_c2c_z, reinterpret_cast<fftw_complex*>(IN), reinterpret_cast<fftw_complex*>(field_Z));
  
  // Z -> Y
  MPI_Alltoallw( field_Z , send_count_yz, send_disp_yz, subarrays_C_Z,
                 field_Yz, send_count_yz, send_disp_yz, subarrays_C_Yz, comm_yz);
                 
  // FFT Y
  fftw_execute(plan_c2c_fwrd.plan_c2c_y);
  
  // Y -> X
  MPI_Alltoallw( field_Yx, send_count_xy, send_disp_xy, subarrays_C_Yx,
                 field_X , send_count_xy, send_disp_xy, subarrays_C_X, comm_xy);
  
  // FFT X
  fftw_execute_dft(plan_c2c_fwrd.plan_c2c_x, reinterpret_cast<fftw_complex*>(field_X), reinterpret_cast<fftw_complex*>(OUT));
  
}

void MikeFFT::c2c_bwrd(std::complex<double>* IN, std::complex<double>* OUT)
{
  
  // FFT X
  fftw_execute_dft(plan_c2c_bwrd.plan_c2c_x, reinterpret_cast<fftw_complex*>(IN), reinterpret_cast<fftw_complex*>(field_X));
  
  // X -> Y
  MPI_Alltoallw( field_X , send_count_xy, send_disp_xy, subarrays_C_X ,
                 field_Yx, send_count_xy, send_disp_xy, subarrays_C_Yx, comm_xy);
                 
  // FFT Y
  fftw_execute(plan_c2c_bwrd.plan_c2c_y);
  
  // Y -> Z
  MPI_Alltoallw( field_Yz, send_count_yz, send_disp_yz, subarrays_C_Yz,
                 field_Z , send_count_yz, send_disp_yz, subarrays_C_Z , comm_yz);
                 
  // FFT Z
  fftw_execute_dft(plan_c2c_bwrd.plan_c2c_z, reinterpret_cast<fftw_complex*>(field_Z), reinterpret_cast<fftw_complex*>(OUT));
  
  // normalize
  double N3_inv = 1./(N*N*N);
  for(int id = 0; id < size_Cz_tot;id++)
  {
    OUT[id] *= N3_inv;
  }

}
