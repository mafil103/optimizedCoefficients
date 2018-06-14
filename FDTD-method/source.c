#include "ezinc.h"

// variables for the pulse definition and position
static double widthX, widthY, widthZ, ppw,cdtds;
static double locationX = 200.0, locationY = 75, locationZ = 75;

// query to enter parameters for the pulse
 void ezIncInit(Grid *g){
    cdtds = Cdtds;
    printf("Enter kx: ");
    scanf(" %lf",&ppw);
    printf("Enter widthX: ");
    scanf(" %lf", &widthX);
    printf("Enter widthY: ");
    scanf(" %lf", &widthY);
    printf("Enter widthZ: ");
    scanf(" %lf", &widthZ);
    return;
 }

// calculate simple gaussian pulse
 double ezInc(double x1, double y1, double z1) {
     //return exp(-pow((time - delay - location / cdtds) / width, 2));
     return cos(ppw*x1)*exp(-pow((x1 - locationX) / widthX, 2)-pow((y1 - locationY) / widthY, 2)-pow((z1 - locationZ) / widthZ, 2));
 }
