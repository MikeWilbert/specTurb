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
//~ const double DT = 0.005*(64./NUM);
//~ const double DT = 0.01*(64./NUM);
const double DT = 1.e-3; 

// output directory
const std::string OUT_DIR = "/home/fs1/mw/Turbulence/Tests/restart";
// output interval
const double OUT_INTERVAL = 0.25;
// simulation time
const double END_SIMU = 3.;

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
const double NU  = 0.05;
//~ const double NU  = 0.0012;
// magnetic diffusivity
const double ETA = NU*1.;
//~ const double ETA = NU*1.e-1;

/** DEFINES **/

//~ #define NS
# define RESTART true
//~ # define RESTART false

/** Resolutions **/
/*
 * MHD
 * N =   64 -> nu = 0.045
 * N =  128 -> nu = 0.018
 * N =  256 -> nu = 0.0072
 * N =  512 -> nu = 0.00285
 * N = 1024 -> nu = 0.0012
 * 
 * NS
 * N =   64 -> nu = 0.05
 * N =  128 -> nu = 0.019
 * N =  256 -> nu = 0.0077
 * N =  512 -> nu = 0.003
 * N = 1024 -> nu = 0.0012
 */
