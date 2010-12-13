#include <math.h>
#include "nr.h"
#define EPS 3.0e-11          /* EPS is the relative precision. */

/*
 Given the lower and upper limits of integration x1 and x2, and given
 n, this routine returns arrays x[1..n] and w[1..n] of length n,
 containing the abscissas and weights of the Gauss- Legendre n-point
 quadrature formula.  */
void gauleg(real x1, real x2, real x[], real w[], int n)
{    
  int m,j,i;
  double z1,z,xm,xl,pp,p3,p2,p1;    /* High precision is a good idea for this rou-
				       tine. */
  m=(n+1)/2;            /*  The roots are symmetric in the interval, so */
  xm=0.5*(x2+x1);       /*  we only have to find half of them. */
  xl=0.5*(x2-x1);
  for (i=1;i<=m;i++) {  /*   Loop over the desired roots. */
    z=cos(3.141592654*(i-0.25)/(n+0.5));
    /* Starting with the above approximation to the ith root, we enter
       the main loop of refinement by Newton's method. */
    do {p1=1.0;
    p2=0.0;
    for (j=1;j<=n;j++) {    /*  Loop up the recurrence relation to get the */
      p3=p2;                /*  Legendre polynomial evaluated at z. */
      p2=p1;
      p1=((2.0*j-1.0)*z*p2-(j-1.0)*p3)/j;
    }
  /* p1 is now the desired Legendre polynomial. We next compute pp,
  its derivative, by a standard relation involving also p2, the
  polynomial of one lower order. */
    pp=n*(z*p1-p2)/(z*z-1.0);
    z1=z;
    z=z1-p1/pp;             /*  Newton's method. */
    } while (fabs(z-z1) > EPS);
    x[i]=xm-xl*z;           /*  Scale the root to the desired interval, */
    x[n+1-i]=xm+xl*z;       /* and put in its symmetric counterpart. */
    w[i]=2.0*xl/((1.0-z*z)*pp*pp);  /*    Compute the weight */
    w[n+1-i]=w[i];          /*  and its symmetric counterpart. */
  }
}


/*
  =========================================================
*/
real func(real x,real y,real z);
real yy1(real x);
real yy2(real x);
real z1(real x,real y);
real z2(real x,real y);


static real xsav,ysav;
static real (*nrfunc)(real,real,real);

/* Returns the integral of a user-supplied function func over a
   three-dimensional region specified by the limits x1, x2, and by the
   user-supplied functions yy1, yy2, z1, and z2, as defined in
   (4.6.2). (The functions y1 and y2 are here called yy1 and yy2 to
   avoid conflict with the names of Bessel functions in some C
   libraries). Integration is performed by calling qgaus recursively.  */

real quad3d(real (*func)(real, real, real), real x1, real x2){ 
  real qgaus(real (*func)(real), real a, real b);
  real f1(real x);
  nrfunc=func;
  return qgaus(f1,x1,x2);
}

/*    This is H of eq. (4.6.5).*/
real f1(real x) {    
  real qgaus(real (*func)(real), real a, real b);
  real f2(real y);
  real yy1(real),yy2(real);

  xsav=x;
  return qgaus(f2,yy1(x),yy2(x));
}


/* This is G of eq. (4.6.4). */
real f2(real y) {    
  real qgaus(real (*func)(real), real a, real b);
  real f3(real z);
  real z1(real,real),z2(real,real);

  ysav=y;
  return qgaus(f3,z1(xsav,y),z2(xsav,y));
}

/* The integrand f(x, y, z) evaluated at fixed x and y. */
real f3(real z)    
{    
  return (*nrfunc)(xsav,ysav,z);
}

/* Returns the integral of the function func between a and b, by
   ten-point Gauss-Legendre inte- gration: the function is evaluated
   exactly ten times at interior points in the range of integration.  */

real qgaus(real (*func)(real), real a, real b) {    
  int j;
  real xr,xm,dx,s;
  /* The abscissas and weights. First value of each array not used. */
  static real x[]={0.0,0.1488743389,0.4333953941, 
		    0.6794095682,0.8650633666,0.9739065285};    
  static real w[]={0.0,0.2955242247,0.2692667193,    
		    0.2190863625,0.1494513491,0.0666713443};

  xm=0.5*(b+a);
  xr=0.5*(b-a);
  s=0;    /* Will be twice the average value of the function, since the */
  for (j=1;j<=5;j++) {  /* ten weights (five numbers above each used twice) */
    dx=xr*x[j];         /* sum to 2. */
    s += w[j]*((*func)(xm+dx)+(*func)(xm-dx));
  }
  return s *= xr;       /* Scale the answer to the range of integration.*/
}
