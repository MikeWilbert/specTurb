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
const int NUM = 1024;
// processor grid
int const PDIMS[2] = {96,96};
//~ int const PDIMS[2] = {48,48};
// CFL number
//~ const double DT = 0.005*(64./NUM);
//~ const double DT = 0.01*(64./NUM);
const double DT = 1.e-8; 

// output directory
const std::string OUT_DIR = "/p/scratch/specdyn/Turbulence/Production/ForcedTurbulence_N1024";
// output interval
const double OUT_INTERVAL = 0.5;
// simulation time
const double END_SIMU = 80.;
// restart
const int RESTART_STEP = 60;
//~ #define RESTART_DIR  "/p/scratch/specdyn/Turbulence/Restart_tests/Restart_1024_longint"

// choose initial setup: 
// (0):all zero; 
// (1):Orszag-Tang; 
// (2):Energy spectrum [s=11/3] with linear Forcing
const int SETUP = 2;
const int BACKGROUND = 0;

// choose Forcing
// (0) Alvelius
// (2) Titov
const int FORCING = 0;

// domain size
const double LENGTH = PI2;
// kinematic viscosity
const double NU  = 0.0012;
//~ const double NU  = 0.0005;
// magnetic diffusivity
const double ETA = NU*1.;
//~ const double ETA = NU*1.e-1;

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
