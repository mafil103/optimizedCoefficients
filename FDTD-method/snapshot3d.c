/* modified template from the source
[4]	Understanding the Finite-Difference Time-Domain Method,
John B. Schneider, www.eecs.wsu.edu/~schneidj/ufdtd, 2010.
*/

#include <stdio.h>
#include <stdlib.h>
#include "fdtd-macro.h"

static int temporalStride = -2, frameZ = 0, frameY = 0,startTime;
static char basename[80];

// initializing the snapshot method with the necessary variables
void snapshot3dInit(Grid *g) {

    int choice;

    // Activate snapshot module
    printf("Do you want 2D snapshots of the 3D grid? (1=yes, 0=no) ");
    scanf("%d", &choice);
    if (choice == 0) {
        temporalStride = -1;
        return;
    }
    // enter time of first snapshot happening
    printf("Duration of simulation is %d steps.\n", MaxTime);
    // enter interval for the following snapshots
    printf("Enter start time and temporal stride: ");
    scanf(" %d %d", &startTime, &temporalStride);
    printf("Enter the base name: ");
    scanf(" %s", basename);

    return;
}


void snapshot3d(double dx, Grid *g) {
  int mm, nn, pp;
  float dim1, dim2, temp;
  char filename[100];
  FILE *out;
  // variable to save analytical solution
  int tmp = ((Cdtds/dx)*Time);

  // checking if entered variables are reasonable
  if (temporalStride == -1) {
    return;
  } if (temporalStride < -1) {
    fprintf(stderr,
      "snapshot2d: snapshotInit2d must be called before snapshot.\n"
      "            Temporal stride must be set to positive value.\n");
    exit(-1);
  }

// control if limitations are full filled
  if (/*(Time == 0 ||Time== 499|| Time == 999){*/Time >= startTime &&
      (Time - startTime) % temporalStride == 0) {
    /************ write the constant-z slice ************/
    sprintf(filename, "%s-z.%d", basename, frameZ++);
    out = fopen(filename, "wb");

    // define dimensions of the collected data
    dim2 = SizeY;
    dim1 = SizeX-1;
    fwrite(&dim1, sizeof(float), 1, out);
    fwrite(&dim2, sizeof(float), 1, out);

    // write ez-field into the file
    pp = SizeY / 2;
    for (nn = SizeY - 1; nn >= 0; nn--)
      for (mm = 0; mm < SizeX - 1; mm++) {
	temp = (float)Ez(mm, nn, pp); // save ez-field in a temporal variable
	fwrite(&temp, sizeof(float), 1, out); // save the float in a raw-file
      }

    fclose(out);  // close file

    // write ez-component in a fixed window following the analytical solution (WIP)
    // snapshot on analyitical solution
/*    sprintf(filename, "%s-y.%d", basename, frameY++);
    out = fopen(filename, "wb");

    if(Time == 0)
    {
        tmp=0;
    } else if (Time>400)
    {
        tmp=-65;
    }
    if (Time>900)
    {
        tmp=130;
    }
    dim1 = 76;
    dim2 = 76;
    fwrite(&dim1, sizeof(float), 1, out);
    fwrite(&dim2, sizeof(float), 1, out);

    pp = SizeY / 2;
    for (nn = 113; nn >= 37; nn--)
      for (mm = 162+tmp; mm < 238+tmp; mm++) {
	temp = (float)Ez(mm, nn, pp);
	fwrite(&temp, sizeof(float), 1, out);
      }

    fclose(out);  // close file
*/
  }

  return;
}
