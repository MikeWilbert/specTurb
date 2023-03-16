#pragma once

#include "define.h"
#include "MikeFFT.h"

class CSpecDyn
{
  
  private:
  
    //parameters
    const int N;
    const int* pdims;
    const double dt;

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
    double time;
  
    // MPI
    int myRank, nprocs;
    int mpi_coords[2];
    MPI_Comm comm;
    
    MPI_Datatype vti_subarray;
    MPI_Datatype vti_subarray_vector;
    MPI_Datatype vti_float3;
    
    // MikeFFT
    MikeFFT FFT;
  
    // local dimensions
    int start_R[3], size_R[3];
    int start_F[3], size_F[3];
    
    int size_R_tot;  
    int size_F_tot;  
    
    // output
    int print_count = 0;
    
    // random
    std::mt19937 normal_eng;
    std::normal_distribution<double> normal;
    
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
    
    // help field for RK3
    CX* Vx_F;  
    CX* Vy_F;  
    CX* Vz_F;
    CX* Vx_F1;  
    CX* Vy_F1;  
    CX* Vz_F1;
    CX* Vx_F2;  
    CX* Vy_F2;  
    CX* Vz_F2;
    CX* Bx_F;  
    CX* By_F;  
    CX* Bz_F;
    CX* Bx_F1;  
    CX* By_F1;  
    CX* Bz_F1;
    CX* Bx_F2;  
    CX* By_F2;  
    CX* Bz_F2;
    CX* RHS_Vx_F;
    CX* RHS_Vy_F;
    CX* RHS_Vz_F;
    CX* RHS_Vx_F1;
    CX* RHS_Vy_F1;
    CX* RHS_Vz_F1;
    CX* RHS_Vx_F2;
    CX* RHS_Vy_F2;
    CX* RHS_Vz_F2;
    CX* RHS_Bx_F;
    CX* RHS_By_F;
    CX* RHS_Bz_F;
    CX* RHS_Bx_F1;
    CX* RHS_By_F1;
    CX* RHS_Bz_F1;
    CX* RHS_Bx_F2;
    CX* RHS_By_F2;
    CX* RHS_Bz_F2;
    
    double* RHS_Vx_R;
    double* RHS_Vy_R;
    double* RHS_Vz_R;
    double* RHS_Bx_R;
    double* RHS_By_R;
    double* RHS_Bz_R;
    
    double* Wx_R;  
    double* Wy_R;  
    double* Wz_R;
    double* Jx_R;  
    double* Jy_R;  
    double* Jz_R;
    CX* Wx_F;  
    CX* Wy_F;  
    CX* Wz_F;
    CX* Jx_F;  
    CX* Jy_F;  
    CX* Jz_F;
    
    // vti output
    float* float_array;
    float* float_array_vector;
    
    // Energy Spectrum
    int N_bin = 50;
    
    double* energySpectrum_V;
    double* energySpectrum_V_loc;
    int*   bin_counter_V;
    int*   bin_counter_V_loc;
    double* energySpectrum_B;
    double* energySpectrum_B_loc;
    int*   bin_counter_B;
    int*   bin_counter_B_loc;
    
    // private methods
    void setup_k();
    void setup_fields();
    
    void time_step();
    void calc_RHS(CX* RHSV_X, CX* RHSV_Y, CX* RHSV_Z, CX* V_X, CX* V_Y, CX* V_Z,
                  CX* RHSB_X, CX* RHSB_Y, CX* RHSB_Z, CX* B_X, CX* B_Y, CX* B_Z,
                  double del_t);
    void diffusion_correction(CX* Vx, CX* Vy, CX* Vz, CX* Bx, CX* By, CX* Bz, double del_t);
    void projection(CX* fieldX, CX* fieldY, CX* fieldZ);
    void dealias(CX* fieldX, CX* fieldY, CX* fieldZ);
    
    void fFFT(double* IN_x, double* IN_y, double* IN_z, CX* OUT_x, CX* OUT_y, CX* OUT_z);
    void bFFT(CX* IN_x, CX* IN_y, CX* IN_z, double* OUT_x, double* OUT_y, double* OUT_z);
    
    void print();
    void print_Energy();
    void print_EnergySpectrum();
    void print_scales();
    void print_vti();
    void print_mpi_scalar(double* field, int& N_bytes_scalar, const char* file_name);
    void print_mpi_vector(double* field_X, double* field_Y, double* field_Z, int& N_bytes_vector, const char* file_name);
    
    void calc_Energy(double& energy_V, double& diss_V);
    
  public:
  
    // public methods
    CSpecDyn();
    
    void execute();
    void finalize();
    
};
