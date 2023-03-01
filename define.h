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
const double DT = 0.005*(64./NUM);

// output directory
const std::string OUT_DIR = "/home/fs1/mw/Turbulence/Tests/PowerSpectrum";
// output interval
const double OUT_INTERVAL = 0.1;
// simulation time
const double END_SIMU = 3.;

// choose initial setup: 
// (0):all zero; 
// (1):Orszag-Tang; 
// (2):Random with energy spectrum [s=11/3]
const int SETUP = 2;

// domain size
const double LENGTH = PI2;
// kinematic viscosity
const double NU  = 0.002;
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
