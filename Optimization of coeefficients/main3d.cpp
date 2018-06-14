#include <iostream >
#include <math.h>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <vector>
#include <ctime>
#include <algorithm>
#include <iterator>
#include <cmath>
#include "discpp.h"
#include <nlopt.h>
#include <string>

using namespace std;

// Initialize variables for minimization process
typedef struct {
    double a, b, hx, hy, hz, ht, owurzel, tORx, xZ, yZ, zZ;
} my_constraint_data;

// Nature constants
#define PI 3.14159265
#define c 299792458

// Variables of program limits
const int gridSIZE = 40;
const double weight = 1.0;
const double dt_multiplier = 1.0;
int counts = 0;

// Abbrevation of A in 3d: Ax = alpha_X + (2 * beta _xy * cos (k_y*hy) + (2 * beta _xz * cos (k_z*hz) + (delta_x * ( 1 + 2 * cos( k_X*hx))) equation out of [1]
double A_spatial3(double x01, double x02, double y01, double y02, double z01, double z02, double alpha0, double beta01, double beta02, double delta0)
{
    double l,m,n;
    l=cos(x01*x02);
    m=cos(y01*y02);
    n=cos(z01*z02);
    return alpha0 + (2 * beta01 * l) + (2 * beta02 * m) + (delta0 *(1+(2*n)));
}

// Function for the condition to fulfill the limit of d to delta for vanishing grid steps in 3d: alpha_X = 1 - 2 * beta_xy - 2 * beta_xz - 3 * delta_X   equation out of [1]
double alpha_spatial3(double beta1,double beta2, double delta)
{
    return 1-(2*beta1)-(2*beta2)-(3*delta);
}

// Function to calculate the abbreviation of s_x, s_y, s_z and s_omega: s_n = sin( 0.5 * k_n * hn)/hn   equation out of [1]
double sAbb(double x, double y)
{
    return sin(0.5*x*y)/y;
}

// Delta-value according to Lehe  equation out of [2]
double deltaX_0 (double DT)
{
    return 0.25*(1.0-((1)/((DT*DT)))*(sin((PI*DT)/(2))));
}

// Beta_value according to Lehe  equation out of [2]
double betaLEHE (double DX, double DY)
{
    return (DX*DX)/(8*DY*DY);
}

// CFL-condition in 3d  equation out of [2]
double CFL3d(double a1, double a2, double a3, double bxy, double bxz, double byx, double byz, double bzx, double bzy, double d1, double d2, double d3, double dx, double dy, double dz, double argumentWurzel)
{
    double cfl, t[7];
    t[0]= (a3 + 2 * bzx + 2 * bzy - dz)/((dz/dx)*(dz/dx));
    t[1]= (a2 + 2 * byx + 2 * byz - dy)/((dy/dx)*(dy/dx));
    t[2]= (a2 + 2 * byx - 2 * byz - dy)/((dy/dx)*(dy/dx))+(a3 + 2 * bzx - 2 * bzy - dz)/((dz/dx)*(dz/dx));
    t[3]= (a1 + 2 * bxy + 2 * bxz - d1);
    t[4]= (a1 + 2 * bxy - 2 * bxz - d1) + (a3 - 2 * bzx + 2 * bzy - d3)/((dz/dx)*(dz/dx));
    t[5]= (a1 - 2 * bxy + 2 * bxz - d1) + (a2 - 2 * byx + 2 * byz - d2)/((dy/dx)*(dy/dx));
    t[6]= (a1 - 2 * bxy - 2 * bxz - d1) + (a1 - 2 * byx - 2 * byz - d2)/((dy/dx)*(dy/dx)) + (a3 - 2 * bzx - 2 * bzy - d3)/((dz/dx)*(dz/dx));

    cfl = std::min({ t[0], t[1], t[3], t[4], t[5], t[6]});

    if(0 < argumentWurzel)
    {
        return cfl;
    } else
    {
        return argumentWurzel;
    }
}

// Time-step constraint for minimization process  equation out of [1]
double constraintDt(double cflVal, double dx, double dt, double maximum)
{
    double ot;

    ot = dt_multiplier/maximum - dt;

    return ot;
}

// Dispersion relation in 3d: s_omega^2 = s_X^2 * A_X + s_Y^2 * A_Y + s_Z^2 * A_Z  equation out of [1]
double dispersion3d(double alpha1, double alpha2, double alpha3, double beta1, double beta2, double beta3, double beta4, double beta5, double beta6, double delta1, double delta2, double delta3, double k1, double k2, double k3, double step1, double step2, double step3, int x, int y, int z)
{
    double s1[3], A1[3];

    A1[0]=A_spatial3(k2*y, step2, k3*z, step3, k1*x, step1, alpha1, beta1, beta2, delta1);
    A1[1]=A_spatial3(k1*x, step1, k3*z, step3, k2*y, step2, alpha2, beta3, beta4, delta2);
    A1[2]=A_spatial3(k1*x, step1, k2*y, step2, k3*z, step3, alpha3, beta5, beta6, delta3);
    s1[0]=sAbb(k1*x, step1);
    s1[1]=sAbb(k2*y, step2);
    s1[2]=sAbb(k3*z, step3);
    return (s1[0]*s1[0]*A1[0])+(s1[1]*s1[1]*A1[1])+(s1[2]*s1[2]*A1[2]);
}

// Calculation of omega resulting out of the dispersion relation
double omegaCalc(double disp, double dt)
{
    return (2.0/dt)*asin((dt)*sqrt(disp));
}

// Norm-function in 3d for minimization process  equation out of [1]
double NormFunc3d(double weight, double beta1, double beta2, double beta3, double beta4, double beta5, double beta6, double delta1, double delta2, double delta3, double dx, double dy, double dz, double dt)
{
    double f0,f1, omega2, argWurzel2, alpha2[3], k1[3];

    k1[0]=PI/(gridSIZE*dx);
    k1[1]=PI/(gridSIZE*dy);
    k1[2]=PI/(gridSIZE*dz);

    alpha2[0]=alpha_spatial3(beta1, beta2, delta1);
    alpha2[1]=alpha_spatial3(beta3, beta4, delta2);
    alpha2[2]=alpha_spatial3(beta5, beta6, delta3);
    f1=0;

    // norm function integral
    for (int x = 0; x < gridSIZE  ; x++)
    {
        for (int y = 0; y < gridSIZE ; y++)
        {
            for (int z = 0; z < gridSIZE ; z++)
            {
                argWurzel2 = dispersion3d(alpha2[0], alpha2[1], alpha2[2], beta1, beta2, beta3, beta4, beta5, beta6, delta1, delta2, delta3, k1[0], k1[1], k1[2], dx, dy, dz, x, y, z);
                omega2 = omegaCalc(argWurzel2, dt);
                f0=weight*((omega2-sqrt((k1[0] * k1[0] * x * x) + (k1[1] * k1[1] * y * y) + (k1[2] * k1[2] * z * z)))*(omega2-sqrt((k1[0]* k1[0] * x * x) + (k1[1] * k1[1] * y * y) + (k1[2] * k1[2] * z * z))));
                f1=f1+f0;
            }
        }
    }
    f1=f1*(PI/(gridSIZE-1))*(PI/(gridSIZE-1))*(PI/(gridSIZE-1))*(dy/dx)*(dy/dx)*(dz/dx)*(dz/dx);
    return f1;
}

// Modified Norm-function to be used with Nlopt libary containing function and derivation of the function in 3d
 double NormFuncOpti3d(unsigned n, const double *x, double *grad, void *my_func_data)
{
    my_constraint_data *f = (my_constraint_data *) my_func_data;
    double a = f->a, b = f->b, hx = f->hx, hy = f->hy, hz = f->hz;

    double func;
    double d = PI/100000;
    ++counts;

    // derivatives of the norm function
    if (grad) {
        grad[0] = (NormFunc3d(weight, x[1], x[2], x[3], x[4], x[5], x[6], x[7], x[8], x[9], hx, hy, hz, x[0]+d)-NormFunc3d(weight, x[1], x[2], x[3], x[4], x[5], x[6], x[7], x[8], x[9], hx, hy, hz, x[0]))/d;
        grad[1] = (NormFunc3d(weight, x[1]+d, x[2], x[3], x[4], x[5], x[6], x[7], x[8], x[9], hx, hy, hz, x[0])-NormFunc3d(weight, x[1], x[2], x[3], x[4], x[5], x[6], x[7], x[8], x[9], hx, hy, hz, x[0]))/d;
        grad[2] = (NormFunc3d(weight, x[1], x[2]+d, x[3], x[4], x[5], x[6], x[7], x[8], x[9], hx, hy, hz, x[0])-NormFunc3d(weight, x[1], x[2], x[3], x[4], x[5], x[6], x[7], x[8], x[9], hx, hy, hz, x[0]))/d;
        grad[3] = (NormFunc3d(weight, x[1], x[2], x[3]+d, x[4], x[5], x[6], x[7], x[8], x[9], hx, hy, hz, x[0])-NormFunc3d(weight, x[1], x[2], x[3], x[4], x[5], x[6], x[7], x[8], x[9], hx, hy, hz, x[0]))/d;
        grad[4] = (NormFunc3d(weight, x[1], x[2], x[3], x[4]+d, x[5], x[6], x[7], x[8], x[9], hx, hy, hz, x[0])-NormFunc3d(weight, x[1], x[2], x[3], x[4], x[5], x[6], x[7], x[8], x[9], hx, hy, hz, x[0]))/d;
        grad[5] = (NormFunc3d(weight, x[1], x[2], x[3], x[4], x[5]+d, x[6], x[7], x[8], x[9], hx, hy, hz, x[0])-NormFunc3d(weight, x[1], x[2], x[3], x[4], x[5], x[6], x[7], x[8], x[9], hx, hy, hz, x[0]))/d;
        grad[6] = (NormFunc3d(weight, x[1], x[2], x[3], x[4], x[5], x[6]+d, x[7], x[8], x[9], hx, hy, hz, x[0])-NormFunc3d(weight, x[1], x[2], x[3], x[4], x[5], x[6], x[7], x[8], x[9], hx, hy, hz, x[0]))/d;
        grad[7] = (NormFunc3d(weight, x[1], x[2], x[3], x[4], x[5], x[6], x[7]+d, x[8], x[9], hx, hy, hz, x[0])-NormFunc3d(weight, x[1], x[2], x[3], x[4], x[5], x[6], x[7], x[8], x[9], hx, hy, hz, x[0]))/d;
        grad[8] = (NormFunc3d(weight, x[1], x[2], x[3], x[4], x[5], x[6], x[7], x[8]+d, x[9], hx, hy, hz, x[0])-NormFunc3d(weight, x[1], x[2], x[3], x[4], x[5], x[6], x[7], x[8], x[9], hx, hy, hz, x[0]))/d;
        grad[9] = (NormFunc3d(weight, x[1], x[2], x[3], x[4], x[5], x[6], x[7], x[8], x[9]+d, hx, hy, hz, x[0])-NormFunc3d(weight, x[1], x[2], x[3], x[4], x[5], x[6], x[7], x[8], x[9], hx, hy, hz, x[0]))/d;
    }
    func= NormFunc3d(weight, x[1], x[2], x[3], x[4], x[5], x[6], x[7], x[8], x[9], hx, hy, hz, x[0]);
    printf("fNorm = %g ", func);
    return func;
}

// Constraints for the minimization process in 3d with derivations
double myconstraint3d(unsigned n, const double *x, double *grad, void *data)
{
    my_constraint_data *g = (my_constraint_data *) data;
    double a = g->a, b = g->b, hx = g->hx, hy = g->hy, hz = g->hz, ht = g->ht, owurzel = g->owurzel, tORx = g->tORx, xZ = g->xZ, yZ = g->yZ, zZ = g->zZ;
    double alpha3[3];
    double h4 = PI/100000;
    double kx = PI/(gridSIZE*hx);
    double ky = PI/(gridSIZE*hy);
    double kz = PI/(gridSIZE*hz);

    alpha3[0]=alpha_spatial3(x[1], x[2], x[7]);
    alpha3[1]=alpha_spatial3(x[3], x[4], x[8]);
    alpha3[2]=alpha_spatial3(x[5], x[6], x[9]);

    // necessary derivatives of the constraint function
    if (grad) {
        grad[0] = 2*x[0]*dispersion3d(alpha3[0], alpha3[1], alpha3[2], x[1], x[2], x[3], x[4], x[5], x[6], x[7], x[8], x[9], kx, ky, kz, hx, hy, hz, xZ, yZ, zZ);
        grad[1] = ((x[0]*x[0]*dispersion3d(alpha3[0], alpha3[1], alpha3[2], x[1]+h4, x[2], x[3], x[4], x[5], x[6], x[7], x[8], x[9], kx, ky, kz, hx, hy, hz, xZ, yZ, zZ))-(x[0]*x[0]*dispersion3d(alpha3[0], alpha3[1], alpha3[2], x[1], x[2], x[3], x[4], x[5], x[6], x[7], x[8], x[9], kx, ky, kz, hx, hy, hz, xZ, yZ, zZ)))/h4;
        grad[2] = ((x[0]*x[0]*dispersion3d(alpha3[0], alpha3[1], alpha3[2], x[1], x[2]+h4, x[3], x[4], x[5], x[6], x[7], x[8], x[9], kx, ky, kz, hx, hy, hz, xZ, yZ, zZ))-(x[0]*x[0]*dispersion3d(alpha3[0], alpha3[1], alpha3[2], x[1], x[2], x[3], x[4], x[5], x[6], x[7], x[8], x[9], kx, ky, kz, hx, hy, hz, xZ, yZ, zZ)))/h4;
        grad[3] = ((x[0]*x[0]*dispersion3d(alpha3[0], alpha3[1], alpha3[2], x[1], x[2], x[3]+h4, x[4], x[5], x[6], x[7], x[8], x[9], kx, ky, kz, hx, hy, hz, xZ, yZ, zZ))-(x[0]*x[0]*dispersion3d(alpha3[0], alpha3[1], alpha3[2], x[1], x[2], x[3], x[4], x[5], x[6], x[7], x[8], x[9], kx, ky, kz, hx, hy, hz, xZ, yZ, zZ)))/h4;
        grad[4] = ((x[0]*x[0]*dispersion3d(alpha3[0], alpha3[1], alpha3[2], x[1], x[2], x[3], x[4]+h4, x[5], x[6], x[7], x[8], x[9], kx, ky, kz, hx, hy, hz, xZ, yZ, zZ))-(x[0]*x[0]*dispersion3d(alpha3[0], alpha3[1], alpha3[2], x[1], x[2], x[3], x[4], x[5], x[6], x[7], x[8], x[9], kx, ky, kz, hx, hy, hz, xZ, yZ, zZ)))/h4;
        grad[5] = ((x[0]*x[0]*dispersion3d(alpha3[0], alpha3[1], alpha3[2], x[1], x[2], x[3], x[4], x[5]+h4, x[6], x[7], x[8], x[9], kx, ky, kz, hx, hy, hz, xZ, yZ, zZ))-(x[0]*x[0]*dispersion3d(alpha3[0], alpha3[1], alpha3[2], x[1], x[2], x[3], x[4], x[5], x[6], x[7], x[8], x[9], kx, ky, kz, hx, hy, hz, xZ, yZ, zZ)))/h4;
        grad[6] = ((x[0]*x[0]*dispersion3d(alpha3[0], alpha3[1], alpha3[2], x[1], x[2], x[3], x[4], x[5], x[6]+h4, x[7], x[8], x[9], kx, ky, kz, hx, hy, hz, xZ, yZ, zZ))-(x[0]*x[0]*dispersion3d(alpha3[0], alpha3[1], alpha3[2], x[1], x[2], x[3], x[4], x[5], x[6], x[7], x[8], x[9], kx, ky, kz, hx, hy, hz, xZ, yZ, zZ)))/h4;
        grad[7] = ((x[0]*x[0]*dispersion3d(alpha3[0], alpha3[1], alpha3[2], x[1], x[2], x[3], x[4], x[5], x[6], x[7]+h4, x[8], x[9], kx, ky, kz, hx, hy, hz, xZ, yZ, zZ))-(x[0]*x[0]*dispersion3d(alpha3[0], alpha3[1], alpha3[2], x[1], x[2], x[3], x[4], x[5], x[6], x[7], x[8], x[9], kx, ky, kz, hx, hy, hz, xZ, yZ, zZ)))/h4;
        grad[8] = ((x[0]*x[0]*dispersion3d(alpha3[0], alpha3[1], alpha3[2], x[1], x[2], x[3], x[4], x[5], x[6], x[7], x[8]+h4, x[9], kx, ky, kz, hx, hy, hz, xZ, yZ, zZ))-(x[0]*x[0]*dispersion3d(alpha3[0], alpha3[1], alpha3[2], x[1], x[2], x[3], x[4], x[5], x[6], x[7], x[8], x[9], kx, ky, kz, hx, hy, hz, xZ, yZ, zZ)))/h4;
        grad[9] = ((x[0]*x[0]*dispersion3d(alpha3[0], alpha3[1], alpha3[2], x[1], x[2], x[3], x[4], x[5], x[6], x[7], x[8], x[9]+h4, kx, ky, kz, hx, hy, hz, xZ, yZ, zZ))-(x[0]*x[0]*dispersion3d(alpha3[0], alpha3[1], alpha3[2], x[1], x[2], x[3], x[4], x[5], x[6], x[7], x[8], x[9], kx, ky, kz, hx, hy, hz, xZ, yZ, zZ)))/h4;
    }

    return (x[0]*x[0]*dispersion3d(alpha3[0], alpha3[1], alpha3[2], x[1], x[2], x[3], x[4], x[5], x[6], x[7], x[8], x[9], kx, ky, kz, hx, hy, hz, xZ, yZ, zZ)-a);
}

// Minimization in 3d with optimizing of dt using libary from [6]
double* minimization3dFull(double dt0, double sqrto, double xsqrt, double dx, double dy, double dz, double query, double startT, double func(unsigned n, const double *x, double *grad, void *my_func_data), double arr[])
{
    double ratioY = dy/dx;
    double ratioZ = dz/dx;
    // Limitations of coefficients
    double lb[10] = { 0, -1, -1, -1, -1, -1, -1, -1, -1, -1};
    double ub[10] = { 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};

    nlopt_opt opt;

    opt = nlopt_create(NLOPT_LD_SLSQP, 10);
    nlopt_set_lower_bounds(opt, lb);
    nlopt_set_upper_bounds(opt, ub);

    //constraint condition parameters
    my_constraint_data my_func_data[1] = { {dx, dy, dx, dy, dz} };

    nlopt_set_min_objective(opt, func, my_func_data);

    // constraint conditions
    my_constraint_data data[4] = { {0, 0, dx, dy, dz, dt0, xsqrt -1., 0, 0, 0}, {1, 0, dx, dy, dz, dt0, sqrto, -1, gridSIZE, gridSIZE, gridSIZE}};

    nlopt_add_inequality_constraint(opt, myconstraint3d, &data[0], 1e-8);
    nlopt_add_inequality_constraint(opt, myconstraint3d, &data[1], 1e-8);

    // initial starting values
    double x[10] = { 0.0055, 0.00068*ratioY, 0.00068*ratioZ,  0.00068*ratioZ, 0.00068*ratioY*ratioZ, 0.00068*ratioZ, 0.00068*ratioY*ratioZ,  -0.154, -0.154, -0.154};
    double minf;

    // abort conditions
    nlopt_set_ftol_rel(opt, 0.000005);
    nlopt_set_maxeval(opt, 250);
    //nlopt_set_xtol_rel(opt, 1e-4);

    if (nlopt_optimize(opt, x, &minf) < 0) {
        printf("nlopt failed!\n");
    }
    else {
            // output of solution in command window
        printf("found minimum at f(%g,%g,%g,%g,%g,%g,%g,%g,%g,%g) = %0.10g\n", x[0], x[1], x[2], x[3], x[4], x[5], x[6], x[7], x[8], x[9], minf);
    }
    printf("found minimum after %d evaluations\n", counts);
    // write solution in x array to return to main program
    for (int t = 0; t < 10  ; t++)
    {
        arr[t]=x[t];

    }
    nlopt_destroy(opt);
    return arr;
}

// Minimization in 3d with an entered dt using libary from [6]
double* minimization3dSetT(double dt0, double sqrto, double xsqrt, double dx, double dy, double dz, double query, double startT, double func(unsigned n, const double *x, double *grad, void *my_func_data), double arr[])
{
    double ratioY = dy/dx;
    double ratioZ = dz/dx;
    // Limitations of coefficients
    double lb[10] = { dt0, -1, -1, -1, -1, -1, -1, -1, -1, -1};
    double ub[10] = { dt0, 1, 1, 1, 1, 1, 1, 1, 1, 1};

    nlopt_opt opt;

    opt = nlopt_create(NLOPT_LD_SLSQP, 10);
    nlopt_set_lower_bounds(opt, lb);
    nlopt_set_upper_bounds(opt, ub);

    //constraint condition parameters
    my_constraint_data my_func_data[1] = { {dx, dy, dx, dy, dz} };

    nlopt_set_min_objective(opt, func, my_func_data);

    // initial starting values
    my_constraint_data data[2] = { {0, 0, dx, dy, dz, dt0, xsqrt, -1., 0, 0, 0}, {1, 0, dx, dy, dz, dt0, sqrto, -1, gridSIZE, gridSIZE, gridSIZE}};

    nlopt_add_inequality_constraint(opt, myconstraint3d, &data[0], 1e-8);
    nlopt_add_inequality_constraint(opt, myconstraint3d, &data[1], 1e-8);

    double x[10] = { dt0, 0.00068*ratioY, 0.00068*ratioZ,  0.00068*ratioZ, 0.00068*ratioY*ratioZ, 0.00068*ratioZ, 0.00068*ratioY*ratioZ,  -0.154, -0.154, -0.154};
    double minf;

    // abort conditions
    nlopt_set_ftol_rel(opt, 0.000005);
    nlopt_set_maxeval(opt, 250);
    //nlopt_set_xtol_rel(opt, 1e-4);

    if (nlopt_optimize(opt, x, &minf) < 0) {
        printf("nlopt failed!\n");
    }
    else {
            // output of solution in command window
        printf("found minimum at f(%g,%g,%g,%g,%g,%g,%g,%g,%g,%g) = %0.10g\n", x[0], x[1], x[2], x[3], x[4], x[5], x[6], x[7], x[8], x[9], minf);
    }

    printf("found minimum after %d evaluations\n", counts);
    // write solution in x array to return to main program
    for (int t = 0; t < 10  ; t++)
    {
        arr[t]=x[t];

    }
    nlopt_destroy(opt);
    return arr;
}

// Method to draw v_ph/c in three direction using libary from [8]
void draw3d(double xmatrix[][gridSIZE], double ymatrix[][gridSIZE], double zmatrix[][gridSIZE],  int query, double axis1[], double dt)
{
    for(int pic = 0; pic <2; pic++)
    {
        Dislin g;
        double min1, max1, min2, max2, min3, max3, graphStep1, graphStep2, graphStep3;
        min1 = 10;
        max1 = -1;
        min2 = 10;
        max2 = -1;
        min3 = 10;
        max3 = -1;
        for(int x1 = 0; x1 < gridSIZE  ; x1++)
        {
            for(int y1 = 0; y1 < gridSIZE  ; y1++)
            {
                if(min1>zmatrix[x1][y1])
                {
                    min1=zmatrix[x1][y1];
                }
                else if(max1<zmatrix[x1][y1])
                {
                    max1 = zmatrix[x1][y1];
                }
            }
        }
        graphStep1= (max1-min1)/20;

        if(pic ==0)
        {
            g.metafl ("png");
        } else {
        g.metafl ("cons");
        }
        g.disini ();
        g.pagera ();
        g.hwfont ();

        g.titlin ("Numerical Dispersion 3d", -1);
        if(query == 1)
        {
            g.titlin ("YEE", 2);
        } else if (query ==2)
        {
            g.titlin ("LEHE", 2);
        } else if (query ==3)
        {
            g.titlin ("PUKHOV", 2);
        } else if (query ==4)
        {
            g.titlin ("minimized", 2);
        }
        g.titlin ("dt = ",3);
        std::string str = std::to_string(dt);
        const char * s = str.c_str();
        g.titlin (s,4);

        g.name   ("k_x*delta_x in PI", "x");
        g.name   ("k_y*delta_y in PI", "y");
        g.name   ("v_ph in c", "z");

        g.intax  ();
        g.autres (gridSIZE, gridSIZE);
        g.axspos (300, 1750);
        g.ax3len (700, 700, 700);

        g.labdig (2, "x");
        g.labdig (2, "y");
        g.labdig (2, "z");

        g.graf3  (0.0, 1.0, 0.0, 0.25, 0.0, 1.0, 0.0, 0.25, min1, max1, min1, graphStep1);
        g.crvmat ((double *) zmatrix, gridSIZE, gridSIZE, 1, 1);

        g.title  ();
        g.endgrf();

        for(int x2 = 0; x2 < gridSIZE  ; x2++)
        {
            for(int y2 = 0; y2 < gridSIZE  ; y2++)
            {
                if(min2>xmatrix[x2][y2])
                {
                    min2=xmatrix[x2][y2];
                }
                else if(max2<xmatrix[x2][y2])
                {
                    max2 = xmatrix[x2][y2];
                }
            }
        }
        graphStep2= (max2-min2)/20;

        g.name   ("k_y*delta_y in PI", "x");
        g.name   ("k_z*delta_z in PI", "y");
        g.name   ("v_ph in c", "z");

        g.intax  ();
        g.autres (gridSIZE, gridSIZE);
        g.axspos (1750, 750);
        g.ax3len (700, 700, 700);

        g.labdig (2, "x");
        g.labdig (2, "y");
        g.labdig (2, "z");

        g.graf3  (0.0, 1.0, 0.0, 0.25, 0.0, 1.0, 0.0, 0.25, min2, max2, min2, graphStep2);
        g.crvmat ((double *) xmatrix, gridSIZE, gridSIZE, 1, 1);
        g.endgrf();

        for(int x3 = 0; x3 < gridSIZE  ; x3++)
        {
            for(int y3 = 0; y3 < gridSIZE  ; y3++)
            {
                if(min3>ymatrix[x3][y3])
                {
                    min3=ymatrix[x3][y3];
                }
                else if(max3<ymatrix[x3][y3])
                {
                    max3 = ymatrix[x3][y3];
                }
            }
        }
        graphStep3= (max3-min3)/20;

        g.name   ("k_x*delta_x in PI", "x");
        g.name   ("k_z*delta_z in PI", "y");
        g.name   ("v_ph in c", "z");

        g.intax  ();
        g.autres (gridSIZE, gridSIZE);
        g.axspos (1750, 1750);
        g.ax3len (700, 700, 700);

        g.labdig (2, "x");
        g.labdig (2, "y");
        g.labdig (2, "z");

        g.graf3  (0.0, 1.0, 0.0, 0.25, 0.0, 1.0, 0.0, 0.25, min3, max3, min3, graphStep3);
        g.crvmat ((double *) ymatrix, gridSIZE, gridSIZE, 1, 1);

        g.disfin ();
    }
}

int main ()
{
    // File to save results to be exported into a plotter
    ofstream file1("results.txt");
    ofstream file2("coefficients.txt");

    // Initializing our variables
    int grid = gridSIZE;
    double xmat[grid][gridSIZE];
    double ymat[grid][gridSIZE];
    double zmat[grid][gridSIZE];
    int query, query2;
    double alpha[3], beta[3][3], delta[3], step[4], k[3], v_ph;
    double omegaFinal, fomega, startDT, u;
    double maximumDisp;
    double minimumDisp = 10000000;
    double close;

    // Reading in our variables through a console
    std :: cout << "enter hx" << std :: endl;
    std :: cin >> step[0];
    std :: cout << "enter hy" << std :: endl;
    std :: cin >> step[1];
    std :: cout << "enter hz" << std :: endl;
    std :: cin >> step[2];
    std :: cout << "enter 1 for YEE, 2 for LEHE, 3 for PUKHOV or 4 for minimized" << std :: endl;
    std :: cin >> query;
    std :: cout << "should dt fullfill CFL-condition 0 for no, 1 for yes"  << std :: endl;
    std :: cin >> query2;
    if(query2 == 0)
    {
        std :: cout << "enter dt" << std :: endl;
        std :: cin >> step[3];
    } else if (query2 == 1)
    {
        step[3] = 0.5;
    }

    double axis1[grid];

    for(int p=0;p<grid; p++)
    {
        axis1[p]=(p*PI)/(grid*step[0]);
    }

    // calculation of k according to the grid ratio
    k[0]=PI/(grid*step[0]);
    k[1]=PI/(grid*step[1]);
    k[2]=PI/(grid*step[2]);

    // Constants for various Numerical Dispersion processes
    if (query == 1)
    {
        //YEE
        delta[0]=0;
        delta[1]=0;
        delta[2]=0;
        beta[0][1]=0;
        beta[0][2]=0;
        beta[1][0]=0;
        beta[1][2]=0;
        beta[2][0]=0;
        beta[2][1]=0;
    }
    else if (query == 2)
    {
        //LEHE
        delta[0]=deltaX_0(step[3]);
        delta[1]=0;
        delta[2]=0;
        beta[0][1]=0.125;
        beta[0][2]=0.125;
        beta[1][0]=0.125/((step[1]/step[0])*(step[1]/step[0]));
        beta[1][2]=beta[1][0];
        beta[2][0]=0.125/((step[2]/step[0])*(step[2]/step[0]));
        beta[2][1]=beta[2][0];
    }
    else if (query == 3)
    {
        //PUKHOV/COWAN
        delta[0]=0;
        delta[1]=0;
        delta[2]=0;
        beta[0][1]=0.125/((step[1]/step[0])*(step[1]/step[0]));
        beta[0][2]=0.125/((step[2]/step[0])*(step[2]/step[0]));
        beta[1][0]=0.125;
        beta[1][2]=0.125/((step[2]/step[0])*(step[2]/step[0]));
        beta[2][0]=0.125;
        beta[2][1]=0.125/((step[1]/step[0])*(step[1]/step[0]));
    }
    else if (query == 4)
    {
        //minimized initial values for algorithm
        delta[0]=-0.154;
        delta[1]=-0.154;
        delta[2]=-0.154;
        beta[0][1]=0.068;
        beta[1][0]=0.068;
        beta[0][2]=0.068;
        beta[1][2]=0.068;
        beta[2][0]=0.068;
        beta[2][1]=0.068;
    }
    else
    {
        std :: cout << "ERROR" << std :: endl;
    }

    // Labeling our txt file
    file1 << "# This file is called results.txt"  << endl;
    file1 << "# x/PI\t" <<  "y/PI\t" <<  "z/PI\t" << "v_ph/c"<< endl;
    file2 << "# This file is called coeffiecients.txt"  << endl;
    file2 << "# dt\t" <<  "hx\t" <<  "hy\t" << "hz\t"<< "BETAxy\t"<< "BETAxz\t"<< "BETAyx\t"<< "BETAyz\t"<< "BETAzx\t"<< "BETAzy\t"<< "DELTAx\t"<< "DELTAy\t"<< "DELTAz\t"<< endl;

    // calculation of alpha values with entered coefficients
    alpha[0]=alpha_spatial3(beta[0][1], beta[0][2], delta[0]);
    alpha[1]=alpha_spatial3(beta[1][0], beta[1][2], delta[1]);
    alpha[2]=alpha_spatial3(beta[2][0], beta[2][1], delta[2]);

    // Calculation of dispersion for different dimension
    for (int x5 = 0; x5 < grid ; x5++)
    {
        for (int y5 = 0; y5 < grid ; y5++)
        {
            for (int z5 = 0; z5 < grid ; z5++)
            {
                u = dispersion3d(alpha[0], alpha[1], alpha[2], beta[0][1], beta[0][2], beta[1][0], beta[1][2], beta[2][0], beta[2][1], delta[0], delta[1], delta[2], k[0], k[1], k[2], step[0], step[1], step[2], x5, y5, z5);
                omegaFinal = omegaCalc(u, step[3]);
                if(((k[0] * k[0] * x5 * x5 )+ (k[1] * k[1] * y5 * y5) + (k[2] * k[2] * z5 * z5))==0)
                {
                    v_ph = 1.;
                }
                else
                {
                    v_ph = omegaFinal/(sqrt((k[0]*k[0]*x5*x5)+(k[1]*k[1]*y5*y5)+(k[2]*k[2]*z5*z5)));
                }
                if(maximumDisp<(step[0]*sqrt(u)))
                {
                    maximumDisp = (step[0]*sqrt(u));
                }
                if (minimumDisp>u)
                {
                    minimumDisp = u;
                }

                if(x5==0)
                {
                    xmat[y5][z5]=v_ph;
                }
                if(y5==0)
                {
                    ymat[x5][z5]=v_ph;
                }
                if(z5==0)
                {
                    zmat[x5][y5]=v_ph;
                }
                if(query2 == 0 && query != 4)
                {
                    file1 << (k[0]*x5*step[0])/PI <<"\t"<< (k[1]*y5*step[1])/PI <<"\t"<< (k[2]*z5*step[2])/PI <<"\t"<< v_ph << endl;
                }
            }
        }
        axis1[x5]=(k[0]*x5*step[0])/PI;
    }
    fomega =NormFunc3d(1.0, beta[0][1], beta[0][2], beta[1][0], beta[1][2], beta[2][0], beta[2][1], delta[0], delta[1], delta[2], step[0], step[1], step[2], step[3]);
    if(query2 ==0 && query != 4)
    {
        file1 << "NORM = "<< fomega << endl;
        printf("fNorm = %g ", fomega);
    }

    // Optimization if requested
    if(query == 4 || query2 == 1)
    {
        double saveArr[10];
        double save[2];
        save[0] = CFL3d(alpha[0], alpha[1], alpha[2], beta[0][1], beta[1][0],beta[1][0], beta[1][2], beta[2][0], beta[2][1], delta[0], delta[1], delta[2], step[0], step[1], step[2], minimumDisp);
        save[1] = constraintDt(save[0], step[0], step[3], maximumDisp);
        if (query==4)
        {
            startDT = save[1];
            startDT = min(startDT, 0.0);
            startDT = max(startDT, 1.0);
            if(query2 == 1)
            {
                step[3]=1/sqrt((1/(step[0]*step[0]))+(1/(step[1]*step[1]))+(1/(step[2]*step[2])));
                minimization3dFull(step[3], minimumDisp, maximumDisp, step[0], step[1], step[2], query2, startDT, NormFuncOpti3d, saveArr);
            } else
            {
                minimization3dSetT(step[3], minimumDisp, maximumDisp, step[0], step[1], step[2], query2, startDT, NormFuncOpti3d, saveArr);
            }
            step[3]=saveArr[0];
            beta[0][1]=saveArr[1];
            beta[0][2]=saveArr[2];
            beta[1][0]=saveArr[3];
            beta[1][2]=saveArr[4];
            beta[2][0]=saveArr[5];
            beta[2][1]=saveArr[6];
            delta[0]=saveArr[7];
            delta[1]=saveArr[8];
            delta[2]=saveArr[9];
            alpha[0]=alpha_spatial3(beta[0][1], beta[0][2], delta[0]);
            alpha[1]=alpha_spatial3(beta[1][0], beta[1][2], delta[1]);
            alpha[2]=alpha_spatial3(beta[2][0], beta[2][1], delta[2]);
        }
        // Recalculating the dispersion after optimization
        if(query == 4 || query2 == 1)
        {
            if(query2 == 1&& query != 4)
            {
                    step[3]= 1/maximumDisp;
            }
                for (int x4 = 0; x4 < grid  ; x4++)
                {
                    for (int y4 = 0; y4 < grid ; y4++)
                    {
                        for (int z4 = 0; z4 <  grid ; z4++)
                        {
                            u = dispersion3d(alpha[0], alpha[1], alpha[2], beta[0][1], beta[0][2], beta[1][0], beta[1][2], beta[2][0], beta[2][1], delta[0], delta[1], delta[2], k[0], k[1], k[2], step[0], step[1], step[2], x4, y4, z4);
                            omegaFinal = omegaCalc(u, step[3]);
                            if((k[0]*k[0]*x4*x4+k[1]*k[1]*y4*y4+k[2]*k[2]*z4*z4)==0)
                            {
                                v_ph = 1.;
                            }
                            else
                            {
                                v_ph = omegaFinal/(sqrt((k[0]*k[0]*x4*x4)+(k[1]*k[1]*y4*y4)+(k[2]*k[2]*z4*z4)));
                            }
                            if(x4==0)
                            {
                                xmat[y4][z4]=v_ph;
                            }
                            if(y4==0)
                            {
                                ymat[x4][z4]=v_ph;
                            }
                            if(z4==0)
                            {
                                zmat[x4][y4]=v_ph;
                            }
                            file1 << (k[0]*x4*step[0])/PI <<"\t"<< (k[1]*y4*step[1])/PI <<"\t"<< (k[2]*z4*step[2])/PI << "\t"<< v_ph <<endl;
                        }
                    }
                    axis1[x4]=(k[0]*x4*step[0])/PI;
                }
                fomega = 0;
                fomega =NormFunc3d(1.0, beta[0][1], beta[0][2], beta[1][0], beta[1][2], beta[2][0], beta[2][1], delta[0], delta[1], delta[2], step[0], step[1], step[2], step[3]);
                file1 << "NORM = "<< fomega << endl;
                printf("fNorm = %g ", fomega);
            }
        }
        // save coefficitents into a file
    file2 << step[3] << "\t" << step[0] << "\t" << step[1] << "\t" << step[2] << "\t" << beta[0][1] << "\t" << beta[0][2] << "\t" << beta[1][0] << "\t" << beta[1][2] << "\t" << beta[2][0] << "\t" << beta[2][1] << "\t" << delta[0] << "\t" << delta[1] << "\t" << delta[2] <<endl;
    // output of a heatmap
    draw3d(xmat, ymat, zmat, query, axis1, step[3]);
    std :: cin >> close;
}
