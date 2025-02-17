#pragma once

#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <cstring>
#include <complex>
#include <random>
#include <mpi.h>
#include <sys/stat.h> // mkdir
#include <sys/types.h>   
#include <algorithm> // both for max of several numbers
#include <initializer_list>

/** CONSTANTS **/
const double PI2 = 2.*M_PI;
const double PI  =    M_PI;
typedef std::complex<double> CX;
const CX IM = CX(0., 1.);

/** PARAMETERS **/
// spatial resolution
const int NUM = 128;
// processor grid
int const PDIMS[2] = {8,8};

// forced mode
const double K_F = 1.5;
// width of forcing band
const double DK_F = 0.5;
// resolution of turbulence
const double C_REF = 1.5;
// large eddy turnover time
const double T_LE = 1.;
// magnetic Prandtl number
const double PRM = 1.;
// hyperviscosity
const uint HYP = 2;
// output directory
const std::string OUT_DIR = "/home/fs1/mw/MHD/hyper/MHD_N128_hyper2";

// output interval
const double OUT_INTERVAL = 0.5;
// simulation time
// const double END_SIMU = 20.;
const double END_SIMU = 14.;
// restart
const int RESTART_STEP = 0;
//~ #define RESTART_DIR  "/p/scratch/specturb/Turbulence/Tests/singularity_fix"

// choose initial setup: 
// (0):all zero; 
// (1):Orszag-Tang; 
// (2):Energy spectrum [s=11/3] with linear Forcing
// (3):Read from binary data (.dat)
const int SETUP = 2;
const bool BACKGROUND = false;
const double E0_dE = 0.;
const double dE = 10.;
const double BACKGROUND_ENERGY =  E0_dE * dE;
const std::string BINARY_DIR = "/p/project/specturb/synthetic_fields/Mike/mapping_zw_seed-512/lagrangian_mapping/i0";

// choose Forcing
// (0) None
// (1) Alvelius
const int FORCING = 1;

// domain size
const double LENGTH = PI2;
const double Lz_L = sqrt(1 + E0_dE);

/** DEFINES **/

//~ #define NS
