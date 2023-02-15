#include <stdio.h>
#include <stdlib.h>
#include <string>

#include "CSpecDyn.h"

// main
int main(int argc, char** argv)
{
  
  CSpecDyn simu;
  simu.execute();
  simu.finalize();

  return 0;
  
}
