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
int const PDIMS[2] = {2,4};
// CFL number
const double CFL = 0.5;

// output directory
const std::string OUT_DIR = "/home/fs1/mw/Turbulence/Test1";
// output interval
const double OUT_INTERVAL = 1.;
// simulation time
const double END_SIMU = 4.;

// choose initial setup: (0):all zero; (1):testing purposes
const int SETUP = 1;

// domain size
const double LENGTH = 1.;
// kinematic viscosity
const double NU = 1.e-6;
// magnetic diffusivity
const double ETA = 1.e-6;

/** DEFINES **/

// nothing
