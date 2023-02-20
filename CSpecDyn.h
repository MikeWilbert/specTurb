#pragma once

#include "define.h"
#include "MikeFFT.h"

class CSpecDyn
{
  
  private:
  
    //parameters
    const int N;
    const int* pdims;
    const double cfl;

    const std::string out_dir;
    const double out_interval;
    const double end_simu;

    const double L;
    const double nu;
    const double eta;
    
    const int setup;
    
    // derived quantities
    double XB; // left spatial position
    double dx; // space      discretization
    double dk; // wavenumber discretization
    double dt; // time       discretization
    double time;
  
    // MPI
    int myRank, nprocs;
    int mpi_coords[2];
    MPI_Comm comm;
    
    MPI_Datatype vti_subarray;
    
    // MikeFFT
    MikeFFT FFT;
  
    // local dimensions
    int start_R[3], size_R[3];
    int start_F[3], size_F[3];
    
    int size_R_tot;  
    int size_F_tot;  
    
    // output
    int vti_count = 0;
    
    // fields
    double* kx;
		double* ky;
		double* kz;
		double* k2;
    
    double* Vx_R;  
    double* Vy_R;  
    double* Vz_R;
    
    double* Bx_R;  
    double* By_R;  
    double* Bz_R;
    
    float* float_array;
    float* float_array_vector;
    
    // private methods
    void setup_k();
    void setup_fields();
    
    void print_vti();
    void print_mpi_scalar(double* field, int& N_bytes_scalar, const char* file_name);
    
  public:
  
    // public methods
    CSpecDyn();
    
    void execute();
    void finalize();
    
};
