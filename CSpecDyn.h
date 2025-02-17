#pragma once

#include "define.h"
#include "MikeFFT.h"

class CSpecDyn
{
  
  private:
  
    //parameters
    const int N;        // resolution
    const int* pdims;   // processor grid
    const double k_f;   // forcing wavenumber
    const double dk_f;  // width of forcingwavenumber band
    const double c_ref; // resolution of Kolmogorov scale
    const double T;     // large-eddy turnover-time
    const double Pr_m;  // magnetic Prandtl nuber
    const int    hyp;   // hyperviscosity coefficient (ordinary: hyp = 1)

    std::string out_dir;
    double out_interval;
    double end_simu;

    const int setup;
    
    // derived quantities
    double P;  // energy injection rate
    double nu; // kinematic viscosity
    double eta;// magnetic resistivity 

    double XB; // left spatial position
    double ZB; // left spatial position
    double dx; // space      discretization
    double dz; // space      discretization
    double dk; // wavenumber discretization
    double dkz; // wavenumber discretization
    double k_max;
    double kz_max;
    
    // other quantites
    double time;
    double dt;
    double L;
    double Lz;
  
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
    std::mt19937 angle_eng;
    std::uniform_real_distribution<double> angle;
    std::mt19937 length_eng;
    std::uniform_real_distribution<double> length;
    std::mt19937 normal_eng;
    std::normal_distribution<double> normal;
    
    // fields
    double* kx;
		double* ky;
		double* kz;
		double* k2;
		double* k2h;
    
    CX* Vx_R;  
    CX* Vy_R;  
    CX* Vz_R;
    
    CX* Bx_R;  
    CX* By_R;  
    CX* Bz_R;
    
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
    
    CX* RHS_Vx_R;
    CX* RHS_Vy_R;
    CX* RHS_Vz_R;
    CX* RHS_Bx_R;
    CX* RHS_By_R;
    CX* RHS_Bz_R;
    
    CX* Wx_R;  
    CX* Wy_R;  
    CX* Wz_R;
    CX* Jx_R;  
    CX* Jy_R;  
    CX* Jz_R;
    CX* Wx_F;  
    CX* Wy_F;  
    CX* Wz_F;
    CX* Jx_F;  
    CX* Jy_F;  
    CX* Jz_F;
    
    CX* B0x;
    CX* B0y;
    CX* B0z;
    
    // random Forces
    CX* Force_X;
    CX* Force_Y;
    CX* Force_Z;
    
    // vti output
    float* float_array;
    float* float_array_vector;
    
    // Energy Spectrum
    int N_bin = 100;
    
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
    void set_dt();
    void calc_RHS(CX* RHSV_X, CX* RHSV_Y, CX* RHSV_Z, CX* V_X, CX* V_Y, CX* V_Z,
                  CX* RHSB_X, CX* RHSB_Y, CX* RHSB_Z, CX* B_X, CX* B_Y, CX* B_Z,
                  double del_t);
    void projection(CX* fieldX, CX* fieldY, CX* fieldZ);
    void dealias(CX* fieldX, CX* fieldY, CX* fieldZ);
    
    void fFFT(CX* IN_x, CX* IN_y, CX* IN_z, CX* OUT_x, CX* OUT_y, CX* OUT_z);
    void bFFT(CX* IN_x, CX* IN_y, CX* IN_z, CX* OUT_x, CX* OUT_y, CX* OUT_z);
    
    void print();
    void print_Energy();
    void print_scales();
    void print_vti();
    void print_mpi_scalar(CX* field, long& N_bytes_scalar, const char* file_name);
    void print_mpi_vector(CX* field_X, CX* field_Y, CX* field_Z, long& N_bytes_vector, const char* file_name);
    
    void calc_Energy(double& energy_V, double& diss_V);
    void calc_Energy(double& energy_V, double& diss_V, double& energy_B, double& diss_B);
    
    void Alvelius();
    
    void restart();
    void read_binary();
    
  public:
  
    // public methods
    CSpecDyn();
    
    void execute();
    void finalize();
    
};
