/* update-equation defined by
[1] Blinne, Alexander & Schinkel, David & Kuschel, Stephan & Elkina, Nina & G. Rykovanov, Sergey & Zepf, Matt. (2017).
A Systematic Approach to Numerical Dispersion in Maxwell Solvers.
Computer Physics Communications. 10.1016/j.cpc.2017.10.010.
*/

#include "fdtd-macro.h"
#include <stdio.h>
#include <math.h>

static double beta[3][3], delta[3], alpha[3];
static double c = 1; /* Speed of light */

//  magnetic field update equations
void updateH(double bxy, double bxz, double byx, double byz, double bzx, double bzy, double d1, double d2, double d3, double dx, double dy, double dz, Grid *g) {
    int mm, nn, pp;
    beta[0][1] = bxy;
    beta[0][2] = bxz;
    beta[1][0] = byx;
    beta[1][2] = byz;
    beta[2][0] = bzx;
    beta[2][1] = bzy;
    delta[0] = d1;
    delta[1] = d2;
    delta[2] = d3;
    alpha[0] = 1 - (2 * beta[0][1]) - (2 * beta[0][2]) - (3 * delta[0]);
    alpha[1] = 1 - (2 * beta[1][0]) - (2 * beta[1][2]) - (3 * delta[1]);
    alpha[2] = 1 - (2 * beta[2][0]) - (2 * beta[2][1]) - (3 * delta[2]);

    for (mm = 1; mm < SizeX - 1; mm++){
        for (nn = 1; nn < SizeY - 2; nn++){
            for (pp = 1; pp < SizeZ - 2; pp++){
                Hx(mm, nn, pp) = Hx(mm, nn, pp) + Chxe(mm, nn, pp)
                  * ((((alpha[2]/dz)   * (Ey(mm  , nn  , pp+1) - Ey(mm  , nn  , pp  )))
                    + ((delta[2]/dz)   * (Ey(mm  , nn  , pp+2) - Ey(mm  , nn  , pp-1)))
                    + ((beta[2][0]/dz) * (Ey(mm+1, nn  , pp+1) - Ey(mm+1, nn  , pp  )))
                    + ((beta[2][0]/dz) * (Ey(mm-1, nn  , pp+1) - Ey(mm-1, nn  , pp  )))
                    + ((beta[2][1]/dz) * (Ey(mm  , nn+1, pp+1) - Ey(mm  , nn+1, pp  )))
                    + ((beta[2][1]/dz) * (Ey(mm  , nn-1, pp+1) - Ey(mm  , nn-1, pp  ))))
                    -(((alpha[1]/dy)   * (Ez(mm  , nn+1, pp  ) - Ez(mm  , nn  , pp  )))
                    + ((delta[1]/dy)   * (Ez(mm  , nn+2, pp  ) - Ez(mm  , nn-1, pp  )))
                    + ((beta[1][0]/dy) * (Ez(mm+1, nn+1, pp  ) - Ez(mm+1, nn  , pp  )))
                    + ((beta[1][0]/dy) * (Ez(mm-1, nn+1, pp  ) - Ez(mm-1, nn  , pp  )))
                    + ((beta[1][2]/dy) * (Ez(mm  , nn+1, pp+1) - Ez(mm  , nn  , pp+1)))
                    + ((beta[1][2]/dy) * (Ez(mm  , nn+1, pp-1) - Ez(mm  , nn  , pp-1)))));
            }}}
    for (mm = 1; mm < SizeX - 2; mm++){
        for (nn = 1; nn < SizeY - 1; nn++){
            for (pp = 1; pp < SizeZ - 2; pp++){
                Hy(mm, nn, pp) = Hy(mm, nn, pp) + Chye(mm, nn, pp)
                  * ((((alpha[0]/dx)   * (Ez(mm+1, nn  , pp  ) - Ez(mm  , nn  , pp  )))
                    + ((delta[0]/dx)   * (Ez(mm+2, nn  , pp  ) - Ez(mm-1, nn  , pp  )))
                    + ((beta[0][1]/dx) * (Ez(mm+1, nn+1, pp  ) - Ez(mm  , nn+1, pp  )))
                    + ((beta[0][1]/dx) * (Ez(mm+1, nn-1, pp  ) - Ez(mm  , nn-1, pp  )))
                    + ((beta[0][2]/dx) * (Ez(mm+1, nn  , pp+1) - Ez(mm  , nn  , pp+1)))
                    + ((beta[0][2]/dx) * (Ez(mm+1, nn  , pp-1) - Ez(mm  , nn  , pp-1))))
                    -(((alpha[2]/dz)   * (Ex(mm  , nn  , pp+1) - Ex(mm  , nn  , pp  )))
                    + ((delta[2]/dz)   * (Ex(mm  , nn  , pp+2) - Ex(mm  , nn  , pp-1)))
                    + ((beta[2][0]/dz) * (Ex(mm+1, nn  , pp+1) - Ex(mm+1, nn  , pp  )))
                    + ((beta[2][0]/dz) * (Ex(mm-1, nn  , pp+1) - Ex(mm-1, nn  , pp  )))
                    + ((beta[2][1]/dz) * (Ex(mm  , nn+1, pp+1) - Ex(mm  , nn+1, pp  )))
                    + ((beta[2][1]/dz) * (Ex(mm  , nn-1, pp+1) - Ex(mm  , nn-1, pp  )))));
            }}}
    for (mm = 1; mm < SizeX - 2; mm++){
        for (nn = 1; nn < SizeY - 2; nn++){
            for (pp = 1; pp < SizeZ - 1; pp++){
                Hz(mm, nn, pp) = Hz(mm, nn, pp) + Chze(mm, nn, pp)
                  * ((((alpha[1]/dy)   * (Ex(mm  , nn+1, pp  ) - Ex(mm  , nn  , pp  )))
                    + ((delta[1]/dy)   * (Ex(mm  , nn+2, pp  ) - Ex(mm  , nn-1, pp  )))
                    + ((beta[1][0]/dy) * (Ex(mm+1, nn+1, pp  ) - Ex(mm+1, nn  , pp  )))
                    + ((beta[1][0]/dy) * (Ex(mm-1, nn+1, pp  ) - Ex(mm-1, nn  , pp  )))
                    + ((beta[1][2]/dy) * (Ex(mm  , nn+1, pp+1) - Ex(mm  , nn  , pp+1)))
                    + ((beta[1][2]/dy) * (Ex(mm  , nn+1, pp-1) - Ex(mm  , nn  , pp-1))))
                    -(((alpha[0]/dx)   * (Ey(mm+1, nn  , pp  ) - Ey(mm  , nn  , pp  )))
                    + ((delta[0]/dx)   * (Ey(mm+2, nn  , pp  ) - Ey(mm-1, nn  , pp  )))
                    + ((beta[0][1]/dx) * (Ey(mm+1, nn+1, pp  ) - Ey(mm  , nn+1, pp  )))
                    + ((beta[0][1]/dx) * (Ey(mm+1, nn-1, pp  ) - Ey(mm  , nn-1, pp  )))
                    + ((beta[0][2]/dx) * (Ey(mm+1, nn  , pp+1) - Ey(mm  , nn  , pp+1)))
                    + ((beta[0][2]/dx) * (Ey(mm+1, nn  , pp-1) - Ey(mm  , nn  , pp-1)))));
            }}}

    return;
}

//  electric field update equations
void updateE(double dx, double dy, double dz, Grid *g) {
    int mm, nn, pp;
    for (mm = 1; mm < SizeX - 2; mm++)
        for (nn = 2; nn < SizeY - 2; nn++)
            for (pp = 2; pp < SizeZ - 2; pp++)
                Ex(mm, nn, pp) = Ex(mm, nn, pp) + Cexh(mm, nn, pp) * c * c
                                * (((1/dy)*(Hz(mm, nn, pp) - Hz(mm, nn-1, pp)))
                                - ((1/dz)*(Hy(mm, nn, pp) - Hy(mm, nn, pp-1))));

    for (mm = 2; mm < SizeX - 2; mm++)
        for (nn = 1; nn < SizeY - 2; nn++)
            for (pp = 2; pp < SizeZ - 2; pp++)
                Ey(mm, nn, pp) =  Ey(mm, nn, pp) + Ceyh(mm, nn, pp) * c * c
                                * (((1/dz)*(Hx(mm, nn, pp) - Hx(mm,nn,pp-1)))
                                -((1/dx)*(Hz(mm, nn, pp) - Hz(mm-1,nn,pp))));

    for (mm = 2; mm < SizeX - 2; mm++)
        for (nn = 2; nn < SizeY - 2; nn++)
            for (pp = 1; pp < SizeZ - 2; pp++)
                Ez(mm, nn, pp) =  Ez(mm, nn, pp) + Cezh(mm, nn, pp) * c * c
                                * (((1/dx)*(Hy(mm, nn, pp) - Hy(mm-1,nn,pp)))
                                -((1/dy)*(Hx(mm, nn, pp) - Hx(mm,nn-1,pp))));

    return;
}
