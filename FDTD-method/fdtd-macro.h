/* macros used to simplify the fdtd-method by
[4]	Understanding the Finite-Difference Time-Domain Method,
John B. Schneider, www.eecs.wsu.edu/~schneidj/ufdtd, 2010.
*/

#ifndef _FDTD_MACRO_H
#define _FDTD_MACRO_H

#include "fdtd-grid.h"

// 3 dimensional grid macros
#define HxG(G,MM,NN,PP)   G->hx[((MM)*(SizeYG(G)-1) + NN)*(SizeZG(G)-1) + PP]
#define ChxhG(G,MM,NN,PP) G->chxh[((MM)*(SizeYG(G)-1) + NN)*(SizeZG(G)-1) + PP]
#define ChxeG(G,MM,NN,PP) G->chxe[((MM)*(SizeYG(G)-1) + NN)*(SizeZG(G)-1) + PP]

#define HyG(G,MM,NN,PP)   G->hy[((MM)*SizeYG(G) + NN)*(SizeZG(G)-1) + PP]
#define ChyhG(G,MM,NN,PP) G->chyh[((MM)*SizeYG(G) + NN)*(SizeZG(G)-1) + PP]
#define ChyeG(G,MM,NN,PP) G->chye[((MM)*SizeYG(G) + NN)*(SizeZG(G)-1) + PP]

#define HzG(G,MM,NN,PP)   G->hz[((MM)*(SizeYG(G)-1) + NN)*SizeZG(G) + PP]
#define ChzhG(G,MM,NN,PP) G->chzh[((MM)*(SizeYG(G)-1) + NN)*SizeZG(G) + PP]
#define ChzeG(G,MM,NN,PP) G->chze[((MM)*(SizeYG(G)-1) + NN)*SizeZG(G) + PP]

#define ExG(G,MM,NN,PP)   G->ex[((MM)*SizeYG(G) + NN)*SizeZG(G) + PP]
#define CexeG(G,MM,NN,PP) G->cexe[((MM)*SizeYG(G) + NN)*SizeZG(G) + PP]
#define CexhG(G,MM,NN,PP) G->cexh[((MM)*SizeYG(G) + NN)*SizeZG(G) + PP]

#define EyG(G,MM,NN,PP)   G->ey[((MM)*(SizeYG(G)-1) + NN)*SizeZG(G) + PP]
#define CeyeG(G,MM,NN,PP) G->ceye[((MM)*(SizeYG(G)-1) + NN)*SizeZG(G) + PP]
#define CeyhG(G,MM,NN,PP) G->ceyh[((MM)*(SizeYG(G)-1) + NN)*SizeZG(G) + PP]

#define EzG(G,MM,NN,PP)   G->ez[((MM)*SizeYG(G) + NN)*(SizeZG(G)-1) + PP]
#define CezeG(G,MM,NN,PP) G->ceze[((MM)*SizeYG(G) + NN)*(SizeZG(G)-1) + PP]
#define CezhG(G,MM,NN,PP) G->cezh[((MM)*SizeYG(G) + NN)*(SizeZG(G)-1) + PP]

#define SizeXG(G)      G->sizeX
#define SizeYG(G)      G->sizeY
#define SizeZG(G)      G->sizeZ
#define TimeG(G)       G->time
#define MaxTimeG(G)    G->maxTime
#define CdtdsG(G)      G->cdtds
#define TypeG(G)       G->type

// one dimensional grid
#define Hy1(MM)     Hy1G(g,MM)
#define Chyh1(MM)   Chyh1G(g,MM)
#define Chye1(MM)   Chye1G(g,MM)

#define Ez1(MM)     Ez1G(g,MM)
#define Ceze1(MM)   Ceze1G(g,MM)
#define Cezh1(MM)   Cezh1G(g,MM)

// 3 dimensional grid by structure Grid
#define Hx(MM,NN,PP)   HxG(g,MM,NN,PP)
#define Chxh(MM,NN,PP) ChxhG(g,MM,NN,PP)
#define Chxe(MM,NN,PP) ChxeG(g,MM,NN,PP)

#define Hy(MM,NN,PP)   HyG(g,MM,NN,PP)
#define Chyh(MM,NN,PP) ChyhG(g,MM,NN,PP)
#define Chye(MM,NN,PP) ChyeG(g,MM,NN,PP)

#define Hz(MM,NN,PP)   HzG(g,MM,NN,PP)
#define Chzh(MM,NN,PP) ChzhG(g,MM,NN,PP)
#define Chze(MM,NN,PP) ChzeG(g,MM,NN,PP)

#define Ex(MM,NN,PP)   ExG(g,MM,NN,PP)
#define Cexe(MM,NN,PP) CexeG(g,MM,NN,PP)
#define Cexh(MM,NN,PP) CexhG(g,MM,NN,PP)

#define Ey(MM,NN,PP)   EyG(g,MM,NN,PP)
#define Ceye(MM,NN,PP) CeyeG(g,MM,NN,PP)
#define Ceyh(MM,NN,PP) CeyhG(g,MM,NN,PP)

#define Ez(MM,NN,PP)   EzG(g,MM,NN,PP)
#define Ceze(MM,NN,PP) CezeG(g,MM,NN,PP)
#define Cezh(MM,NN,PP) CezhG(g,MM,NN,PP)

#define SizeX       SizeXG(g)
#define SizeY       SizeYG(g)
#define SizeZ       SizeZG(g)
#define Time        TimeG(g)
#define MaxTime     MaxTimeG(g)
#define Cdtds       CdtdsG(g)
#define Type        TypeG(g)

#endif
