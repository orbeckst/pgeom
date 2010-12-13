/* 
   $Id: cylint.h,v 1.5 2002/11/07 14:55:15 oliver Exp $
   $Log: cylint.h,v $
   Revision 1.5  2002/11/07 14:55:15  oliver
   - some volume calculation stuff was not checked in
   - added METHYL (distance 0.14nm) to species
   - rewrote species id using enum

   Revision 1.4  2002/08/27 02:20:49  oliver
   works but contains lots of dirt hacks in cylint.c...

   Revision 1.3  2002/08/25 22:19:51  oliver
   - substitued float with real; compiling with -DDOUBLE uses double precision for real, otherwise it is float
   - Boolean is now in the new xgtypes.h header (as is real)

   Revision 1.2  2002/08/24 15:45:54  oliver
   - Added partition function Zsum()
   - Thermodynamics avearge
   - cylint: use explicit initialisation routines
   - exact test case for the thermodynamice average
   - shifted LJ potential (yet another way to compute the volume)
   - set T from commandline
   BUT: volume is still not correct: Problem: Hard walls at Rmax are not accounted for in the density of states (see Labbook II,p81)

   Revision 1.1  2002/08/23 18:12:52  oliver
   Header file for all the things that come with the determination of the
   volume (mostly integration...)

   Integrate a volume in cylindrical coordinates over a cylinder. This
   is a simplified version of 'Numerical Recipes in C, p164: quad3d' +
   the Gauss-Legendre quadrature (p152)

   I specify my boundaries as numbers not functions because of
   cylindrical symmetry.

   How to use the n-point quadrature:

      init_gaussleg(n1);
      cylint(...)
      cylint(...)
      free_gaussleg()

      init_gaussleg(n2)
      ...

*/

#ifndef CYLINT_H
#define CYLINT_H

#include <math.h>
#include "xgtypes.h"
#include "util.h"

#define GAUSSLEG_EPS 3.0e-11          /* EPS is the relative precision. */
#define kBoltzmann 8.314510e-03 /* kJ mol^-1 K^-1 */

#define RAD  0
#define PHI  1
#define ZZZ  2 

struct potpars {    /* parameters for the potential and the volume calculations */
  real (*u1)(real,real,real);  /* potential function for one interaction centre */
  real Temp;                      /* temperature in Kelvin */
  real Nfreedom;                  /* degrees of freedom of the particle */
  real c6, c12;                   /* LJ parameters in kJ mol^-1 nm^-6/-12 */
  real min;                       /* minimum within inner integration
                                     volume; must be found with  find_min() first ! */
  int   ncenters;                  /* number of LJ centers */
  real **centers;                 /* LJ center ( cylindrical) coordinates */
  real **xyzcenters;               /* LJ center ( cartesian) coordinates */
  int   ngaussleg;                 /* number of points in the Gauss-Legendre quadrature per
                                      dimension */
};

/* accessible volume -- the final solution (Labbook II,p84); with 
   major help by Andrew Horsfield */
real xV (real (*U)(real,real,real), real Rmax, real zmin, real zmax);


/* canonical partition function (Zustandssumme) Tr(exp(-beta H)) 
   (configurational part H == U)
 */
real Zsum (real (*U)(real,real,real), real Rmax, real zmin, real zmax);

/* <A> = Tr[A*exp(-beta H)]/Tr[exp(-beta H)] */
real tdavg (real (*A)(real,real,real,real), real (*U)(real,real,real), 
	     real Rmax, real zmin, real zmax);


/* integration in cylindrical coordinates on cylindrical volumes */
extern real cylint(real (*func)(real, real, real), 
	     real r1, real r2, real phi1, real phi2,
	     real z1, real z2);

/* initialise the global variables needed in vljcyl */
extern void init_vljcyl (struct potpars *);

/* -rho*exp(-beta Sum_i V_LJ(|r(rho,phi,z)-R_i|)) */
extern real vljcyl (real r, real phi, real z); 

/* V_LJ again but takes a rvec (but still in cylindrical coordinates)
   as argument (for the minimisation this is more convenient) */
extern real vljvec (rvec);

/* V_LJ(r) - minvlj */
extern real vljshiftcyl (real r, real phi, real z);

/* classical volume operator Theta(E-V) 
   first three args are dummies, 4th is V, E is a global and is
   initialized with init_vljcyl() */
extern real Qvol (real,real,real,real);

/* Zero potential f(r) = 0 */
extern real Zero (real,real,real);

/* f(r) = 6r (cylindrical symm!) */
extern real LinCheck (real,real,real);

/* configurational part of the partition sum i.e. trace over exp(-beta
   SumV(r)) (same as vljcyl which is a misnomer),*/
extern real Z_V (real r, real phi, real z); 

/* (configurational) trace of Theta(E-V(r)),  betaE=N/2, N:DoF */
extern real Q (real r, real phi, real z); 

/* V_LJ split into Repulsive and Attractive part (WCA) */
extern real vRljcyl (real r, real phi, real z); 
extern real vAljcyl (real r, real phi, real z); 



/*************************************************************
 * find the minimum (uses simplex search (amoeba) from NR)
 * if v0 is NULL, it takes a default starting point
 */
real find_min (real (*f)(rvec), real Rmax, real zmin, real zmax,rvec v0);

/* Gauss-Legendre quadrature on n points; cylint calls qngaus which
   knows through global variables (after init_gaussleg) its order and
   abscissas and weights */
extern void init_gaussleg (int);
extern void free_gaussleg (void);

#endif /* CYLINT_H */
