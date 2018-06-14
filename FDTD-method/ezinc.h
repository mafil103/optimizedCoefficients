/* ez-field header
[4]	Understanding the Finite-Difference Time-Domain Method,
John B. Schneider, www.eecs.wsu.edu/~schneidj/ufdtd, 2010.
*/
#ifndef _EZINC_H
#define _EZINC_H

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "fdtd-macro.h"

void ezIncInit(Grid *g);
double ezInc(double x1, double y1, double z1);

#endif
