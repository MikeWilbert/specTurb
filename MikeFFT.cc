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

plan1D::plan1D(fftw_plan plan_y_r2c_, fftw_plan plan_y_c2r_, fftw_plan plan_z_r2c_, fftw_plan plan_z_c2r_) :
plan_y_r2c( plan_y_r2c_),  plan_y_c2r( plan_y_c2r_),
plan_z_r2c( plan_z_r2c_),  plan_z_c2r( plan_z_c2r_)
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
  
  split(N    , cart_coords[0], pdims[0], size_rot_Ry[0], start_rot_Ry[0]);
  split(N    ,              0,        1, size_rot_Ry[1], start_rot_Ry[1]);
  split(N    , cart_coords[1], pdims[1], size_rot_Ry[2], start_rot_Ry[2]);
  
  split(N    , cart_coords[0], pdims[0], size_rot_Fy[0], start_rot_Fy[0]);
  split(N/2+1,              0,        1, size_rot_Fy[1], start_rot_Fy[1]);
  split(N    , cart_coords[1], pdims[1], size_rot_Fy[2], start_rot_Fy[2]);
  
  size_Fx_tot = size_Fx[0]*size_Fx[1]*size_Fx[2];
  size_Fy_tot = size_Fy[0]*size_Fy[1]*size_Fy[2];
  size_Fz_tot = size_Fz[0]*size_Fz[1]*size_Fz[2];
  size_Rz_tot = size_Rz[0]*size_Rz[1]*size_Rz[2];
  size_rot_Ry_tot = size_rot_Ry[0]*size_rot_Ry[1]*size_rot_Ry[2];
  size_rot_Fy_tot = size_rot_Fy[0]*size_rot_Fy[1]*size_rot_Fy[2];
  
  // init subarrays
  subarrays_X  = new MPI_Datatype[pdims[0]];
  subarrays_Yx = new MPI_Datatype[pdims[0]];
  subarrays_Yz = new MPI_Datatype[pdims[1]];
  subarrays_Z  = new MPI_Datatype[pdims[1]];
  subarrays_Y_rot  = new MPI_Datatype[pdims[1]];
  subarrays_Z_rot  = new MPI_Datatype[pdims[1]];
  
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
  field_X  = reinterpret_cast<CX*>(fftw_alloc_complex((size_t)size_Fx_tot));
  field_Yx = reinterpret_cast<CX*>(fftw_alloc_complex((size_t)size_Fy_tot));
  field_Yz = reinterpret_cast<CX*>(fftw_alloc_complex((size_t)size_Fy_tot));
  field_Z  = reinterpret_cast<CX*>(fftw_alloc_complex((size_t)size_Fz_tot));
  
  Field_Y_R_rot = fftw_alloc_real((size_t)size_rot_Ry_tot);
  Field_Y_F_rot = reinterpret_cast<CX*>(fftw_alloc_complex((size_t)size_rot_Fy_tot));
  
  // plan FFT
  plan_r2c = make_plan3D_r2c();
  plan_c2r = make_plan3D_c2r();
  plan_rot = make_plan1D();
  
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
    
    split(size_Fx[0], i, pdims[0], size_tmp1[0], start_tmp1[0]);
    split(size_Fy[1], i, pdims[0], size_tmp2[1], start_tmp2[1]);
                                                         
    MPI_Type_create_subarray(3, size_Fx, size_tmp1, start_tmp1, MPI_ORDER_C, mpiComplex, &subarrays_X [i]);
    MPI_Type_create_subarray(3, size_Fy, size_tmp2, start_tmp2, MPI_ORDER_C, mpiComplex, &subarrays_Yx[i]);
    
    MPI_Type_commit(&subarrays_X [i]);
    MPI_Type_commit(&subarrays_Yx[i]);
  
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
                                                         
    MPI_Type_create_subarray(3, size_Fy, size_tmp1, start_tmp1, MPI_ORDER_C, mpiComplex, &subarrays_Yz[i]);
    MPI_Type_create_subarray(3, size_Fz, size_tmp2, start_tmp2, MPI_ORDER_C, mpiComplex, &subarrays_Z [i]);
    
    MPI_Type_commit(&subarrays_Yz[i]);
    MPI_Type_commit(&subarrays_Z [i]);
  
  }
  
  // real y<->z
  for(int i = 0; i < pdims[1]; i++)
  {
    size_tmp1[0] = size_rot_Ry[0];
    size_tmp1[1] = size_rot_Ry[1];
    size_tmp1[2] = size_rot_Ry[2];
    
    start_tmp1[0] = 0; // starts in partitioned local array!
    start_tmp1[1] = 0;
    start_tmp1[2] = 0;
    
    size_tmp2[0] = size_Rz[0];
    size_tmp2[1] = size_Rz[1];
    size_tmp2[2] = size_Rz[2];
    
    start_tmp2[0] = 0;
    start_tmp2[1] = 0;
    start_tmp2[2] = 0;
    
    split(size_rot_Ry[1], i, pdims[1], size_tmp1[1], start_tmp1[1]);
    split(size_Rz[2]    , i, pdims[1], size_tmp2[2], start_tmp2[2]);
                                                         
    MPI_Type_create_subarray(3, size_rot_Ry, size_tmp1, start_tmp1, MPI_ORDER_C, MPI_DOUBLE, &subarrays_Y_rot[i]);
    MPI_Type_create_subarray(3, size_Rz    , size_tmp2, start_tmp2, MPI_ORDER_C, MPI_DOUBLE, &subarrays_Z_rot[i]);
    
    MPI_Type_commit(&subarrays_Y_rot[i]);
    MPI_Type_commit(&subarrays_Z_rot[i]);
  
  }
}

void MikeFFT::get_sizes(int size_R_[3], int start_R_[3], int size_F_[3], int start_F_[3])
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

void MikeFFT::get_comm(MPI_Comm& comm_)
{
  comm_ = comm;
}

double*               MikeFFT::malloc_R()
{
  return fftw_alloc_real((size_t)size_Rz_tot);
}

std::complex<double>* MikeFFT::malloc_C()
{
  return reinterpret_cast<CX*>(fftw_alloc_complex((size_t)size_Fx_tot));
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
  double* field_R_init  = fftw_alloc_real((size_t)size_Rz_tot);
  fftw_complex* field_C_init  = fftw_alloc_complex((size_t)size_Fx_tot);
  
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
         field_R_init, reinterpret_cast<fftw_complex*>(field_Z), FFTW_MEASURE);
  
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
           reinterpret_cast<fftw_complex*>(field_X), reinterpret_cast<fftw_complex*>(field_C_init), FFTW_FORWARD, FFTW_MEASURE);
  
  fftw_free(field_R_init);
  fftw_free(field_C_init);
  
  return plan3D_r2c(plan_x, plan_y, plan_z);
}

plan3D_c2r MikeFFT::make_plan3D_c2r()
{
  double* field_R_init  = fftw_alloc_real((size_t)size_Rz_tot);
  fftw_complex* field_C_init  = fftw_alloc_complex((size_t)size_Fx_tot);
  
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
           reinterpret_cast<fftw_complex*>(field_C_init), reinterpret_cast<fftw_complex*>(field_X), FFTW_BACKWARD, FFTW_MEASURE);
  
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
         reinterpret_cast<fftw_complex*>(field_Z), field_R_init, FFTW_MEASURE);
           
  fftw_free(field_R_init);
  fftw_free(field_C_init);
           
  return plan3D_c2r(plan_x, plan_y, plan_z);
  
}

plan1D MikeFFT::make_plan1D()
{
  fftw_plan plan_y_f;
  fftw_plan plan_y_b;
  fftw_plan plan_z_f;
  fftw_plan plan_z_b;
  
  double* field_Rz_init  = fftw_alloc_real((size_t)size_Rz_tot);
  
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
  
  // r2c
  dims[0].n  = size_Rz [2]; // dims.n needs to be the size of the real array
  dims[0].is = stride_R[2];
  dims[0].os = stride  [2];
  
  howmany_dims[0].n  = size_Rz [0];
  howmany_dims[0].is = stride_R[0]; 
  howmany_dims[0].os = stride  [0]; 
  
  howmany_dims[1].n  = size_Rz [1];
  howmany_dims[1].is = stride_R[1]; 
  howmany_dims[1].os = stride  [1]; 
  
  plan_z_f = fftw_plan_guru_dft_r2c(1, dims, 2, howmany_dims, 
         field_Rz_init, reinterpret_cast<fftw_complex*>(field_Z), FFTW_MEASURE);
  
  // c2r
  dims[0].n  = size_Rz [2]; // dims.n needs to be the size of the real array
  dims[0].is = stride  [2];
  dims[0].os = stride_R[2];
  
  howmany_dims[0].n  = size_Rz [0]; 
  howmany_dims[0].is = stride  [0]; 
  howmany_dims[0].os = stride_R[0]; 
  
  howmany_dims[1].n  = size_Rz [1];
  howmany_dims[1].is = stride  [1]; 
  howmany_dims[1].os = stride_R[1]; 
  
  plan_z_b = fftw_plan_guru_dft_c2r(1, dims, 2, howmany_dims, 
         reinterpret_cast<fftw_complex*>(field_Z), field_Rz_init, FFTW_MEASURE);

  // Y
  stride_R[0] = size_rot_Ry[1]*size_rot_Ry[2];
  stride_R[1] = size_rot_Ry[2];
  stride_R[2] = 1;
  stride  [0] = size_rot_Fy[1]*size_rot_Fy[2];
  stride  [1] = size_rot_Fy[2];
  stride  [2] = 1;
  
  // r2c
  dims[0].n  = size_rot_Ry [1]; // dims.n needs to be the size of the real array
  dims[0].is = stride_R    [1];
  dims[0].os = stride      [1];
  
  howmany_dims[0].n  = size_rot_Ry [0];
  howmany_dims[0].is = stride_R    [0]; 
  howmany_dims[0].os = stride      [0]; 
  
  howmany_dims[1].n  = size_rot_Ry [2];
  howmany_dims[1].is = stride_R    [2]; 
  howmany_dims[1].os = stride      [2]; 
  
  plan_y_f = fftw_plan_guru_dft_r2c(1, dims, 2, howmany_dims, 
         Field_Y_R_rot, reinterpret_cast<fftw_complex*>(Field_Y_F_rot), FFTW_MEASURE);
  
  // c2r
  dims[0].n  = size_rot_Ry [1]; // dims.n needs to be the size of the real array
  dims[0].is = stride      [1];
  dims[0].os = stride_R    [1];
  
  howmany_dims[0].n  = size_rot_Ry [0]; 
  howmany_dims[0].is = stride      [0]; 
  howmany_dims[0].os = stride_R    [0]; 
  
  howmany_dims[1].n  = size_rot_Ry [2];
  howmany_dims[1].is = stride      [2]; 
  howmany_dims[1].os = stride_R    [2]; 
  
  plan_y_b = fftw_plan_guru_dft_c2r(1, dims, 2, howmany_dims, 
         reinterpret_cast<fftw_complex*>(Field_Y_F_rot), Field_Y_R_rot, FFTW_MEASURE);
           
  fftw_free(field_Rz_init);
  
  return plan1D(plan_y_f, plan_y_b, plan_z_f, plan_z_b);
}

void MikeFFT::r2c(double* IN, std::complex<double>* OUT)
{
  
  // FFT Z
  fftw_execute_dft_r2c(plan_r2c.plan_r2c_z, IN, reinterpret_cast<fftw_complex*>(field_Z));
  
  // Z -> Y
  MPI_Alltoallw( field_Z , send_count_yz, send_disp_yz, subarrays_Z,
                 field_Yz, send_count_yz, send_disp_yz, subarrays_Yz, comm_yz);
                 
  // FFT Y
  fftw_execute(plan_r2c.plan_c2c_y);
  
  // Y -> X
  MPI_Alltoallw( field_Yx, send_count_xy, send_disp_xy, subarrays_Yx,
                 field_X , send_count_xy, send_disp_xy, subarrays_X, comm_xy);
  
  // FFT X
  fftw_execute_dft(plan_r2c.plan_c2c_x, reinterpret_cast<fftw_complex*>(field_X), reinterpret_cast<fftw_complex*>(OUT));
  
}

void MikeFFT::c2r(std::complex<double>* IN, double* OUT)
{
  
  // FFT X
  fftw_execute_dft(plan_c2r.plan_c2c_x, reinterpret_cast<fftw_complex*>(IN), reinterpret_cast<fftw_complex*>(field_X));
  
  // X -> Y
  MPI_Alltoallw( field_X , send_count_xy, send_disp_xy, subarrays_X ,
                 field_Yx, send_count_xy, send_disp_xy, subarrays_Yx, comm_xy);
                 
  // FFT Y
  fftw_execute(plan_c2r.plan_c2c_y);
  
  // Y -> Z
  MPI_Alltoallw( field_Yz, send_count_yz, send_disp_yz, subarrays_Yz,
                 field_Z , send_count_yz, send_disp_yz, subarrays_Z , comm_yz);
                 
  // FFT Z
  fftw_execute_dft_c2r(plan_c2r.plan_c2r_z, reinterpret_cast<fftw_complex*>(field_Z), OUT);
  
  // normalize
  double N3_inv = 1./(N*N*N);
  for(int id = 0; id < size_Rz_tot;id++)
  {
    OUT[id] *= N3_inv;
  }

}

void MikeFFT::rotate_single(double* IN, double* OUT, double angle, double L, double x_beg)
{
  const double dx = L/N;
  const double dk = 2*M_PI/L;
  const double a =  -tan(0.5*angle);
  const double b =   sin(angle);
  const CX IM = CX(0., 1.);
  
  double angle_max = 0.15 * 2.*M_PI;
  
  /*************** Z ***************/
  // FFT Z ->
  fftw_execute_dft_r2c(plan_rot.plan_z_r2c, IN, reinterpret_cast<fftw_complex*>(field_Z));
  
  // shear Z
  for( int ix = 0; ix < size_Fz[0]; ix++){
  for( int iy = 0; iy < size_Fz[1]; iy++){
  for( int iz = 0; iz < size_Fz[2]; iz++){
    
    int id = ix * size_Fz[1] * size_Fz[2] + iy * size_Fz[2] + iz;
    
    double pos_kz = (iz + start_Fz[2]) * dk;
    double pos_y  = (iy + start_Rz[1]) * dx + x_beg;
    
    field_Z[id] *= exp(IM * pos_kz * pos_y * a);
    
    int ik = iz + start_Fz[2];
    //~ field_Z[id] *= 0.5 * (1. + cos( ik *2./N * M_PI ) );
    
  }}}
  
  // FFT Z <-
  fftw_execute_dft_c2r(plan_rot.plan_z_c2r, reinterpret_cast<fftw_complex*>(field_Z), OUT);
  
  /*************** Y ***************/
  // Z -> Y
  MPI_Alltoallw( OUT          , send_count_yz, send_disp_yz, subarrays_Z_rot,
                 Field_Y_R_rot, send_count_yz, send_disp_yz, subarrays_Y_rot, comm_yz);
   
  // FFT Y ->
  fftw_execute_dft_r2c(plan_rot.plan_y_r2c, Field_Y_R_rot, reinterpret_cast<fftw_complex*>(Field_Y_F_rot));
  
  // shear Y
  for( int ix = 0; ix < size_rot_Fy[0]; ix++){
  for( int iy = 0; iy < size_rot_Fy[1]; iy++){
  for( int iz = 0; iz < size_rot_Fy[2]; iz++){
    
    int id = ix * size_rot_Fy[1] * size_rot_Fy[2] + iy * size_rot_Fy[2] + iz;
    
    double pos_ky = (iy + start_rot_Fy[1]) * dk;
    double pos_z  = (iz + start_rot_Ry[2]) * dx + x_beg;
    
    Field_Y_F_rot[id] *= exp(IM * pos_ky * pos_z * b);
    
    int ik = iy + start_rot_Fy[1];
    //~ Field_Y_F_rot[id] *= 0.5 * (1. + cos( ik *2./N * M_PI ) );
    
  }}}
  
  // FFT Y <-
  fftw_execute_dft_c2r(plan_rot.plan_y_c2r, reinterpret_cast<fftw_complex*>(Field_Y_F_rot), Field_Y_R_rot);
                 
  // Y -> Z
  MPI_Alltoallw( Field_Y_R_rot, send_count_yz, send_disp_yz, subarrays_Y_rot,
                 OUT          , send_count_yz, send_disp_yz, subarrays_Z_rot, comm_yz);
  
  /*************** Z ***************/
  
  // FFT Z ->
  fftw_execute_dft_r2c(plan_rot.plan_z_r2c, OUT, reinterpret_cast<fftw_complex*>(field_Z));
  
  // shear Z
  for( int ix = 0; ix < size_Fz[0]; ix++){
  for( int iy = 0; iy < size_Fz[1]; iy++){
  for( int iz = 0; iz < size_Fz[2]; iz++){
    
    int id = ix * size_Fz[1] * size_Fz[2] + iy * size_Fz[2] + iz;
    
    double pos_kz = (iz + start_Fz[2]) * dk;
    double pos_y  = (iy + start_Rz[1]) * dx + x_beg;
    
    field_Z[id] *= exp(IM * pos_kz * pos_y * a);
    
    int ik = iz + start_Fz[2];
    //~ field_Z[id] *= 0.5 * (1. + cos( ik *2./N * M_PI ) );
    
  }}}
  
  // FFT Z <-
  fftw_execute_dft_c2r(plan_rot.plan_z_c2r, reinterpret_cast<fftw_complex*>(field_Z), OUT);
  
  /****************************/
  
  // normalize
  double N3_inv = 1./(N*N*N);
  for(int id = 0; id < size_Rz_tot;id++)
  {
    OUT[id] *= N3_inv;
  }
  
}

void MikeFFT::rotate(double* IN, double* OUT, double angle, double L, double x_beg)
{
  
  double PI_half = 0.5*M_PI;
  double PI_2    = 2. *M_PI;
  // double angle_max =  0.1 * 2.*M_PI; // works 0.15 rotations is the maximum for L = 2.3
  double angle_max =  0.07 * 2.*M_PI; // playing it save
  
  // copy In -> OUT
  for(int id = 0; id < size_Rz_tot; id++)
  {
    OUT[id] = IN[id];
  }
  
  // reduce angle to [-PI,+PI]
  if(angle >  M_PI)
  {
    do
    {
      angle -= PI_2;
    }while(angle >  M_PI);
  }
  else
  if(angle < -M_PI)
  {
    do
    {
      angle += PI_2;
    }while(angle < -M_PI);
  }
  
  // rotate
  if(angle >  angle_max)
  {
    do
    {
      rotate_single(OUT, OUT,  angle_max, L, x_beg);
      angle -= angle_max;
    }while(angle >  angle_max);
  }
  else
  if(angle < -angle_max)
  {
    do
    {
      rotate_single(OUT, OUT, -angle_max, L, x_beg);
      angle += angle_max;
    }while(angle <  angle_max);
  }

  rotate_single(OUT, OUT, angle, L, x_beg);

}
