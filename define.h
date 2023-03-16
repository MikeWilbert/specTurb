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
const double DT = 0.01*(64./NUM);

// output directory
//~ const std::string OUT_DIR = "/p/scratch/specdyn/Turbulence/NoForce_N512_Nu5em4_u0p25_F1e5";
//~ const std::string OUT_DIR = "/p/scratch/specdyn/Turbulence/TG_N512_Nu5em4_u0p25_F1";
const std::string OUT_DIR = "/home/fs1/mw/Turbulence/TestProduction/Forcing";
// output interval
const double OUT_INTERVAL = 1.;
// simulation time
const double END_SIMU = 40.*4.;

// choose initial setup: 
// (0):all zero; 
// (1):Orszag-Tang; 
// (2):Energy spectrum [s=11/3] with linear Forcing
const int SETUP = 2;

// domain size
const double LENGTH = PI2;
// kinematic viscosity
const double NU  = 5.e-4;
// magnetic diffusivity
const double ETA = NU;

/** DEFINES **/

// nothing

/** Resolutions **/
/*
 * N =   64 -> nu = 0.02
 * N =  128 -> nu = 0.007
 * N =  256 -> nu = 0.003
 * N =  512 -> nu = 0.001
 * N = 1024 -> nu = 0.0004
 */
