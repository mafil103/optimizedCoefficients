/*  structure of grid defined by
[4]	Understanding the Finite-Difference Time-Domain Method,
John B. Schneider, www.eecs.wsu.edu/~schneidj/ufdtd, 2010.
*/

#ifndef _FDTD_GRID_H
#define _FDTD_GRID_H

enum GRIDTYPE {oneDGrid, teZGrid, tmZGrid, threeDGrid};

// information contained in Grid
struct Grid {
  double *hx, *chxh, *chxe;
  double *hy, *chyh, *chye;
  double *hz, *chzh, *chze;
  double *ex, *cexe, *cexh;
  double *ey, *ceye, *ceyh;
  double *ez, *ceze, *cezh;
  int sizeX, sizeY, sizeZ;
  int time, maxTime;
  int type;
  double cdtds;
};

typedef struct Grid Grid;

#endif
