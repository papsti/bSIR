/* bSIR_PF.c */
/* Author: Irena Papst */

/* The bSIR gradient in C for use with R */
/* State variables are in units of proportions of individuals
 and "healthy" people are assumed to be all those who are not infected */
/* Compile on the command line with
 R CMD SHLIB bSIR_PF.c */

#include <R.h> // include R library

/* define parameter array globally */

static double parms[6];

/* order of params in array must match the order in which they are passed from R to the DLL */

#define beta parms[0]
#define beta_F parms[1]
#define gamma_F parms[2]
#define r_beta parms[3]
#define gamma parms[4]
#define rbetaF parms[5]

/* define initializer for parms passed from R */

void initmod(void (* odeparms)(int *, double *))
{
    int num=6;
    odeparms(&num, parms);
}


/* define derivatives function */

void derivs (int *neq, // points to the number of equations
             double *t, // points to time
             double *vars, // points to a double precision array of length neq that containing the current value of the state variables
             double *varsdot, // points to an array that contains the calculated derivatives
             double *varsout, // points to a vector whose first nout elements are other output variables (different from the state variables), followed by the double vals passed by parameters rpar when calling the integrator
             int *ip // points to an integer vector with length at least 3
// first element contains the number of output values (should be equal to or larger than nout)
// second element contains the length of yout
// third element contains the length of op
// next are integer values, as passed by the parameter ipar when calling the integrator
)
{
    if (ip[0] <1) error("nout should be at least 1."); // check whether nout is at least 1
    varsdot[0] = - beta*vars[2]*vars[0] - beta_F*vars[0]*vars[2] - rbetaF*beta_F*vars[0]*vars[1] + gamma_F*vars[1]*(1-vars[2]); // equation for S
    varsdot[1] = - r_beta*beta*vars[1]*vars[2] + beta_F*vars[0]*vars[2] + rbetaF*beta_F*vars[0]*vars[1] - gamma_F*vars[1]*(1-vars[2]); // equation for S_F
    varsdot[2] = - gamma*vars[2] + beta*vars[0]*vars[2] + r_beta*beta*vars[1]*vars[2]; // equation for I
    varsdot[3] = beta*vars[0]*vars[2] + r_beta*beta*vars[1]*vars[2]; // equation for Z, the cumulative case count
    
    varsout[0] = 1 - vars[0] - vars[1] - vars[2]; // calculate R(t)
    varsout[1] = varsdot[2]; // export dI/dt (to calculate number of peaks)
}
