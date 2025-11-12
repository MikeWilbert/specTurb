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
const int NUM = 1024;
// processor grid
int const PDIMS[2] = {96,96};

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
// const std::string OUT_DIR = "/p/scratch/specturb/Turbulence/hyper2/N1024_h2_B0";
// const std::string OUT_DIR = "/p/scratch/specturb/Turbulence/hyper2/N1024_h2_B1";
const std::string OUT_DIR = "/p/scratch/specturb/Turbulence/hyper2/N1024_h2_B10";
// const double E0_dE =  0.;
// const double E0_dE =  1.;
const double E0_dE = 10.;
// const int RESTART_STEP = 18;
// const int RESTART_STEP = 20;
const int RESTART_STEP = 17;

// output interval
const double OUT_INTERVAL = 0.5;
// simulation time
const double END_SIMU = 15.;
// restart
// const int RESTART_STEP = ;
//~ #define RESTART_DIR  "/p/scratch/specturb/Turbulence/Tests/singularity_fix"

// choose initial setup: 
// (0):all zero; 
// (1):Orszag-Tang; 
// (2):Energy spectrum [s=11/3] with linear Forcing
// (3):Read from binary data (.dat)
const int SETUP = 2;
const bool BACKGROUND = true;
// const double E0_dE = 0.;
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
