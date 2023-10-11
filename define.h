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
int const PDIMS[2] = {8,8};
// CFL number
const double DT = 1.e-8;  // TODO! 

// output directory
const std::string OUT_DIR = "/home/fs1/jel/data/Mike/TestRuns/Run_2";
// output interval
const double OUT_INTERVAL = 2.;
// simulation time
const double END_SIMU = 1.;
// restart
const int RESTART_STEP = 0;
//~ #define RESTART_DIR  "/p/scratch/specdyn/Turbulence/Restart_tests/Restart_1024_longint"

// choose initial setup: 
// (0):all zero; 
// (1):Orszag-Tang; 
// (2):Energy spectrum [s=11/3] with linear Forcing
// (3):Read from binary data (.dat)
const int SETUP = 3;
const int BACKGROUND = 0;
const std::string BINARY_DIR = "/home/fs1/jel/data/old/r-64/0";

// choose Forcing
// (0) None
// (1) Alvelius
// (2) Titov
const int FORCING = 0;

// domain size
//~ const double LENGTH = PI2;
const double LENGTH = 1.;
// kinematic viscosity
const double NU  = 0.0143;
//~ const double NU  = 0.0012;
//~ const double NU  = 0.0005;
//~ const double NU  = 0.045;
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
