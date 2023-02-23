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
const int NUM = 256;
// processor grid
int const PDIMS[2] = {48,48};
// CFL number
const double DT = 0.005*(64./NUM);

// output directory
const std::string OUT_DIR = "/p/scratch/specdyn/Turbulence/OrszagTang_Nu200";
// output interval
const double OUT_INTERVAL = 1.;
// simulation time
const double END_SIMU = 20.;

// choose initial setup: (0):all zero; (1):Orszag-Tang
const int SETUP = 1;

// domain size
const double LENGTH = PI2;
// kinematic viscosity
const double NU  = 1./200.;
// magnetic diffusivity
const double ETA = NU;

/** DEFINES **/

// nothing
