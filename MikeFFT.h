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

class plan3D_c2c
{
  private:
    
    fftw_plan plan_c2c_x;
    fftw_plan plan_c2c_y;
    fftw_plan plan_c2c_z;
    
    plan3D_c2c(){};
    plan3D_c2c(fftw_plan plan_c2c_x_, fftw_plan plan_c2c_y_, fftw_plan plan_c2c_z_);
    
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
    int size_Cx [3]; // C/R: Complex/Real numbers in position space 
    int size_Cy [3];
    int size_Cz [3];
    int size_Fx [3];
    int size_Fy [3];
    int size_Fz [3];
    int size_Rz [3];
    
    int start_Cx[3];
    int start_Cy[3];
    int start_Cz[3];
    int start_Fx[3];
    int start_Fy[3];
    int start_Fz[3];
    int start_Rz[3];
    
    int size_Cx_tot;
    int size_Cy_tot;
    int size_Cz_tot;
    int size_Fx_tot; 
    int size_Fy_tot; 
    int size_Fz_tot; 
    int size_Rz_tot; 

    // transpose
    MPI_Datatype* subarrays_F_X ;
    MPI_Datatype* subarrays_F_Yx;
    MPI_Datatype* subarrays_F_Yz;
    MPI_Datatype* subarrays_F_Z ;
    MPI_Datatype* subarrays_C_X ;
    MPI_Datatype* subarrays_C_Yx;
    MPI_Datatype* subarrays_C_Yz;
    MPI_Datatype* subarrays_C_Z ;
    
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
    plan3D_r2c plan_r2c;
    plan3D_c2r plan_c2r;
    plan3D_c2c plan_c2c_fwrd;
    plan3D_c2c plan_c2c_bwrd;
    
    // methods
    void split(int N, int rank, int split_size, int& size, int& start);
    void setup_subarrays();
  
    plan3D_r2c make_plan3D_r2c();
    plan3D_c2r make_plan3D_c2r();
    plan3D_c2c make_plan3D_c2c_fwrd();
    plan3D_c2c make_plan3D_c2c_bwrd();

  public:
    MikeFFT(){};
    MikeFFT(int N_, const int* pdims_);
    void get_sizes_real(int size_R_[3], int start_R_[3], int size_F_[3], int start_F_[3]);
    void get_sizes_cplx(int size_R_[3], int start_R_[3], int size_F_[3], int start_F_[3]);
    void get_comm(MPI_Comm& comm_);
    
    double*               malloc_R();
    std::complex<double>* malloc_F();
    std::complex<double>* malloc_R_cplx();
    std::complex<double>* malloc_F_cplx();
    
    void free_R(double* A);
    void free_C(std::complex<double>* A);
    
    void r2c(double* IN, std::complex<double>* OUT);
    void c2r(std::complex<double>* IN, double* OUT);
    void c2c_fwrd(std::complex<double>* IN, std::complex<double>* OUT);
    void c2c_bwrd(std::complex<double>* IN, std::complex<double>* OUT);
  
};
