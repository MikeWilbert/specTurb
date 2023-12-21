#include <stdio.h>
#include <stdlib.h>
#include <string>

#include "CSpecDyn.h"

/* TOOO:
 * 
 * - was ist mit dem deafult in setup_fields() los?
 * - B0=dB Simu (dB aus Simu mit oder ohne Hintergrundfeld?)
 *
 * - typedef für int->long
 * - fft -> rfft & Alvelius explizit im Ortsraum 
 * 
 */

// main
int main(int argc, char** argv)
{
  CSpecDyn simu;
  simu.execute();
  simu.finalize();
  return 0;
}
