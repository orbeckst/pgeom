/* include file for the Numerical Recipes library */

#ifndef NR_H
#define NR_H

#include <float.h>
#include <math.h>

#ifndef HAVE_REAL
#define HAVE_REAL
#ifdef DOUBLE
typedef double real;
#define FTOL   sqrt(DBL_EPSILON)
         /* tolerance for minimisation; cannot do better than
	    sqrt(machine precision) (see Numerical Recipes) This is ca
	    1e-8 for DBL, and 1e-4 for FLOAT */
#else
typedef float real;
#define FTOL   sqrt(FLT_EPSILON)
#endif /* DOUBLE */
#endif /* HAVE_REAL */

/*
  memory allocations & such
*/
#include "nrutil.h"

/* 
   minimisation.
*/
/* Multidimensional minimization of the function funk(x) where
   x[1..ndim] is a vector in ndim dimensions, by the downhill simplex
   method of Nelder and Mead. The matrix p[1..ndim+1] [1..ndim] is
   input. Its ndim+1 rows are ndim-dimensional vectors which are the
   vertices of the starting simplex. Also input is the vector
   y[1..ndim+1], whose components must be pre- initialized to the
   values of funk evaluated at the ndim+1 vertices (rows) of p; and
   ftol the fractional convergence tolerance to be achieved in the
   function value (n.b.!). On output, p and y will have been reset to
   ndim+1 new points all within ftol of a minimum function value, and
   nfunk gives the number of function evaluations taken. 
*/
int amoeba(real **p, real y[], int ndim, real ftol, 
	    real (*funk)(real []), int *nfunk);


/*
  integration in 3D
*/

/*
 Given the lower and upper limits of integration x1 and x2, and given
 n, this routine returns arrays x[1..n] and w[1..n] of length n,
 containing the abscissas and weights of the Gauss- Legendre n-point
 quadrature formula.  
*/
void gauleg(real x1, real x2, real x[], real w[], int n);

/*
  user supplied functions (globals in integrate3d.c)
  =========================================================
*/
real func(real x,real y,real z);
real yy1(real x);
real yy2(real x);
real z1(real x,real y);
real z2(real x,real y);

/* Returns the integral of a user-supplied function func over a
   three-dimensional region specified by the limits x1, x2, and by the
   user-supplied functions yy1, yy2, z1, and z2, as defined in
   (4.6.2). (The functions y1 and y2 are here called yy1 and yy2 to
   avoid conflict with the names of Bessel functions in some C
   libraries). Integration is performed by calling qgaus recursively.  
*/
real quad3d(real (*func)(real, real, real), real x1, real x2);

/* Returns the integral of the function func between a and b, by
   ten-point Gauss-Legendre inte- gration: the function is evaluated
   exactly ten times at interior points in the range of integration.  
*/
real qgaus(real (*func)(real), real a, real b);

#endif /* NR_H */
