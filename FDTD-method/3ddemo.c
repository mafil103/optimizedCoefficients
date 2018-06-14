/* FDTD-Methode used to simulate an pulse propagating on the x-Axis*/
/* using the following source as template to create the fdtd-method for own pulse and optimized coefficients
[4]	Understanding the Finite-Difference Time-Domain Method,
John B. Schneider, www.eecs.wsu.edu/~schneidj/ufdtd, 2010.
*/
#include "fdtd-alloc.h"
#include "fdtd-macro.h"
#include "fdtd-proto.h"
#include "ezinc.h"

//initialize variable for the grid
static double b[3][3], d[3];
static double dx, dy, dz, completion;

int main()
{
    int mm, nn, pp;
    Grid *g;

    printf("Enter spatial steps in order (hx, hy, hz) : ");
    scanf(" %lf", &dx);
    scanf(" %lf", &dy);
    scanf(" %lf", &dz);

    ALLOC_1D(g, 1, Grid);
    gridInit(dx, dy, dz, g);        // initialize 3D grid
    printf("Enter beta values in order (betaXY, betaXZ, betaYX, betaYZ, betaZX, betaXY) : ");
    scanf(" %lf", &b[0][1]);
    scanf(" %lf", &b[0][2]);
    scanf(" %lf", &b[1][0]);
    scanf(" %lf", &b[1][2]);
    scanf(" %lf", &b[2][0]);
    scanf(" %lf", &b[2][1]);
    printf("Enter delta values in order (deltaX, deltaY, deltaZ) : ");
    scanf(" %lf", &d[0]);
    scanf(" %lf", &d[1]);
    scanf(" %lf", &d[2]);

    ezIncInit(g);
    snapshot3dInit(g);  // initialize snapshots

// yee-algorithm
    for (Time = 0; Time < MaxTime; Time++) {
        updateH(b[0][1], b[0][2], b[1][0], b[1][2], b[2][0], b[2][1], d[0], d[1], d[2], dx, dy, dz, g);       // update magnetic fields
        updateE(dx, dy, dz, g);       // update electric fields
        if(Time == 0)
        {
        // Generating initial field with pulse
        for (mm = 5; mm < SizeX-5; mm++){
            for (nn = 5; nn < SizeY-5; nn++){
                for (pp = 5; pp < SizeZ-5; pp++){
                    Hx(mm, nn, pp) = 0;
                    Hz(mm, nn, pp) = 0;
                    Ex(mm, nn, pp) = 0;
                    Ey(mm, nn, pp) = 0;
                    Hy(mm, nn, pp) = (-1.0)*ezInc(mm, nn, pp);
                    Ez(mm, nn, pp) = ezInc(mm, nn, pp);
        }}}}
        snapshot3d(dx, g);    // using snapshot method to generate raw files
        completion = ( (double)(Time+1)/MaxTime )*100.; // displaying the progress of the simulation
        printf( "Computation: %4.3f %% completed \r", completion);
    }

    return 0;
}
