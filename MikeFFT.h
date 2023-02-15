#pragma once

#include <stdio.h>
#include <stdlib.h>
#include <complex>
#include "math.h" // dont use cmath since openAcc cant work with the STL yet ...

#include <mpi.h>

#include "fftw3.h"

class plan3D_r2c
{
  private:
    
    fftw_plan plan_c2c_x;
    fftw_plan plan_c2c_y;
    fftw_plan plan_r2c_z;
    
    plan3D_r2c(){};
    plan3D_r2c(fftw_plan plan_c2c_x_, fftw_plan plan_c2c_y_, fftw_plan plan_r2c_z_);
    
  friend class MikeFFT;
  
};

class plan3D_c2r
{
  private:
    
    fftw_plan plan_c2c_x;
    fftw_plan plan_c2c_y;
    fftw_plan plan_c2r_z;
    
  public:
    plan3D_c2r(){};
    plan3D_c2r(fftw_plan plan_c2c_x_, fftw_plan plan_c2c_y_, fftw_plan plan_c2r_z_);
    
  friend class MikeFFT;
  
};

class plan1D
{
  private:
    
    fftw_plan plan_y_r2c;
    fftw_plan plan_y_c2r;
    fftw_plan plan_z_r2c;
    fftw_plan plan_z_c2r;
    
  public:
    
    plan1D(){};
    plan1D(fftw_plan plan_y_r2c_, fftw_plan plan_y_c2r_, fftw_plan plan_z_r2c_, fftw_plan plan_z_c2r_);
    
  friend class MikeFFT;
      
};

class MikeFFT
{
  
  private:
    // parameters
    int N;
    const int* pdims;
    
    // complex
    typedef std::complex<double> CX;
    MPI_Datatype mpiComplex;
  
    // MPI
    MPI_Comm comm;
    int myRank, nprocs;
    int cart_coords[2];
    
    MPI_Comm comm_xy;
    MPI_Comm comm_yz;
    
    // sizes/starts
    int size_Fx [3];
    int size_Fy [3];
    int size_Fz [3];
    int size_Rz [3];
    int size_rot_Ry[3];
    int size_rot_Fy[3];
    
    int start_Fx[3];
    int start_Fy[3];
    int start_Fz[3];
    int start_Rz[3];
    int start_rot_Ry[3];
    int start_rot_Fy[3];
    
    int size_Fx_tot;
    int size_Fy_tot;
    int size_Fz_tot;
    int size_Rz_tot; 
    int size_rot_Ry_tot;
    int size_rot_Fy_tot;

    // transpose
    MPI_Datatype* subarrays_X ;
    MPI_Datatype* subarrays_Yx;
    MPI_Datatype* subarrays_Yz;
    MPI_Datatype* subarrays_Z ;
    MPI_Datatype* subarrays_Y_rot ;
    MPI_Datatype* subarrays_Z_rot ;
    
    int* send_count_xy;
    int* send_disp_xy ;
    int* send_count_yz;
    int* send_disp_yz ;
  
    // helper fields
    CX* field_X ;
    CX* field_Yx;
    CX* field_Yz;
    CX* field_Z ;
    
    double* Field_Y_R_rot;
    CX*     Field_Y_F_rot;
    
    // plans
    plan3D_r2c plan_r2c;
    plan3D_c2r plan_c2r;
    plan1D     plan_rot;
    
    // methods
    void split(int N, int rank, int split_size, int& size, int& start);
    void setup_subarrays();
  
    plan3D_r2c make_plan3D_r2c();
    plan3D_c2r make_plan3D_c2r();
    plan1D     make_plan1D();     
  
    void rotate_single(double* IN, double* OUT, double angle, double L, double x_beg); // rotation axis: +x, restricted to angle in [-PI/2, PI/2] (if zero padding would be used)
  
  public:
    MikeFFT(){};
    //~ MikeFFT(int N_, int pdims_[2]);
    MikeFFT(int N_, const int* pdims_);
    void get_sizes(int size_R_[3], int start_R_[3], int size_F_[3], int start_F_[3]);
    void get_comm(MPI_Comm& comm_);
    
    double*               malloc_R();
    std::complex<double>* malloc_C();
    void free_R(double* A);
    void free_C(std::complex<double>* A);
    
    void r2c(double* IN, std::complex<double>* OUT);
    void c2r(std::complex<double>* IN, double* OUT);
    
    void rotate(double* IN, double* OUT, double angle, double L, double x_beg); // rotation axis: +x
  
};
