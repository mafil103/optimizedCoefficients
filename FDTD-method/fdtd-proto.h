/* modified prototype functions from the source
[4]	Understanding the Finite-Difference Time-Domain Method,
John B. Schneider, www.eecs.wsu.edu/~schneidj/ufdtd, 2010.
*/

#ifndef _FDTD_PROTO_H
#define _FDTD_PROTO_H

#include "fdtd-grid.h"

// prototype functions to ensure that the functions can be called
void abcInit(Grid *g);
void abc(Grid *g);

void gridInit(double hx, double hy, double hz, Grid *g);

void snapshot3dInit(Grid *g);
void snapshot3d(double dx, Grid *g);

void updateH(double bxy, double bxz, double byx, double byz, double bzx, double bzy, double d1, double d2, double d3, double dx, double dy, double dz, Grid *g);
void updateE(double hx, double hy, double hz, Grid *g);

#endif
