#include <stdio.h>
#include <stdlib.h>
#include <string>

#include "CSpecDyn.h"

/* TOOO:
 * 
 * README: Parameter erklären (einzeln) und die Sache mit den Skalen und dem Forcing herleiten:
 * -> Auch für Hyperviskosität! -> hab iwo ein Markdown, wo das drin stehen sollte!
 * Und kurz sagen, wie man das benutzt.
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
