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
    
    fftw_plan plan_x;
    fftw_plan plan_y;
    fftw_plan plan_z;
    
    plan3D_r2c(){};
    plan3D_r2c(fftw_plan plan_c2c_x_, fftw_plan plan_c2c_y_, fftw_plan plan_r2c_z_);
    
  friend class MikeFFT;
  
};

class plan3D_c2r
{
  private:
    
    fftw_plan plan_x;
    fftw_plan plan_y;
    fftw_plan plan_z;
    
  public:
    plan3D_c2r(){};
    plan3D_c2r(fftw_plan plan_x_, fftw_plan plan_y_, fftw_plan plan_z_);
    
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

    int start_Fx[3];
    int start_Fy[3];
    int start_Fz[3];
    
    int size_Fx_tot;
    int size_Fy_tot;
    int size_Fz_tot;

    // transpose
    MPI_Datatype* subarrays_X ;
    MPI_Datatype* subarrays_Yx;
    MPI_Datatype* subarrays_Yz;
    MPI_Datatype* subarrays_Z ;
    
    int* send_count_xy;
    int* send_disp_xy ;
    int* send_count_yz;
    int* send_disp_yz ;
  
    // helper fields
    CX* field_X ;
    CX* field_Yx;
    CX* field_Yz;
    CX* field_Z ;
    
    // plans
    plan3D_r2c plan_R2F;
    plan3D_c2r plan_F2R;
    
    // methods
    void split(int N, int rank, int split_size, int& size, int& start);
    void setup_subarrays();
  
    plan3D_r2c make_plan3D_R2F();
    plan3D_c2r make_plan3D_F2R();  
  
  public:
    MikeFFT(){};
    MikeFFT(int N_, const int* pdims_);
    void get_sizes(int size_R_[3], int start_R_[3], int size_F_[3], int start_F_[3]);
    void get_comm(MPI_Comm& comm_);
    
    std::complex<double>* malloc();
    void free_R(double* A);
    void free_C(std::complex<double>* A);
    
    void R2F(std::complex<double>* IN, std::complex<double>* OUT);
    void F2R(std::complex<double>* IN, std::complex<double>* OUT);
  
};
