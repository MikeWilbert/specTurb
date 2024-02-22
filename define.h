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
typedef std::complex<double> CX;
const CX IM = CX(0., 1.);

/** PARAMETERS **/
// spatial resolution
const int NUM = 64;
// processor grid
int const PDIMS[2] = {4,4};
// CFL number
const double DT = 1.e-8;
const double CFL_ADV = 0.5;
const double CFL_DIF = 0.025;

// output directory
const std::string OUT_DIR = "/home/fs1/mw/Turbulence/clean_up/N_128";
// output interval
const double OUT_INTERVAL = 100;
// simulation time
const double END_SIMU = 5000.;
// restart
const int RESTART_STEP = 0;
//~ #define RESTART_DIR  "/p/scratch/specdyn/Turbulence/Restart_tests/Restart_1024_longint"

// choose initial setup: 
// (0):all zero; 
// (1):Orszag-Tang; 
// (2):Energy spectrum [s=11/3] with linear Forcing
// (3):Read from binary data (.dat)
const int SETUP = 2;
const bool BACKGROUND = false;
const double BACKGROUND_ENERGY = 5.0;
const std::string BINARY_DIR = "/p/project/specturb/synthetic_fields/Mike/mapping_zw_seed-512/lagrangian_mapping/i0";

// choose Forcing
// (0) None
// (1) Alvelius
const int FORCING = 1;

// domain size
const double LENGTH = PI2;
//~ const double LENGTH = 1.;
// kinematic viscosity
//~ const double NU  = 0.0009; // decaying Turbulence N=512
//~ const double NU  = 0.0077;
const double NU  = 0.0001;
// magnetic diffusivity
const double ETA = NU;

/** DEFINES **/

//~ #define NS

/** Resolutions **/
/*
 * MHD
 * N =   64 -> nu = 0.045
 * N =  128 -> nu = 0.018
 * N =  256 -> nu = 0.0072
 * N =  512 -> nu = 0.00285
 * N = 1024 -> nu = 0.0012
 * N = 1024 -> nu = 0.0005
 * 
 * NS
 * N =   64 -> nu = 0.05
 * N =  128 -> nu = 0.019
 * N =  256 -> nu = 0.0077
 * N =  512 -> nu = 0.003
 * N = 1024 -> nu = 0.0012
 */
