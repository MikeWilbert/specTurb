# TSpecDyn

Pseudo-Spectral solver for MHD turbulence.

## Third-party libraries

- **OpenMPI** (BSD-style license, external)
  - Used for MPI parallelization
  - Users must install OpenMPI on their system
  - Tested with OpenMPI 4.x
  - License: https://www.open-mpi.org/software/ompi/license/

- **FFTW3** (GPL, external)
  - Used for FFTs
  - Users must install FFTW3
  - License: http://www.fftw.org/

- **nlohmann/json** (MIT License, included)
  - Header-only JSON parser included in `include/nlohmann/json.hpp`