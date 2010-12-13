/* 
   $Id: cylint.c,v 1.9 2002/11/07 14:55:14 oliver Exp $
   $Log: cylint.c,v $
   Revision 1.9  2002/11/07 14:55:14  oliver
   - some volume calculation stuff was not checked in
   - added METHYL (distance 0.14nm) to species
   - rewrote species id using enum

   Revision 1.8  2002/08/27 02:20:48  oliver
   works but contains lots of dirt hacks in cylint.c...

   Revision 1.7  2002/08/26 23:42:01  oliver
   - simplex search (from NR), uses amoeba
   - xV and gE -- the proper way of calculating volume

   Revision 1.6  2002/08/25 22:19:51  oliver
   - substitued float with real; compiling with -DDOUBLE uses double precision for real, otherwise it is float
   - Boolean is now in the new xgtypes.h header (as is real)

   Revision 1.5  2002/08/24 15:45:54  oliver
   - Added partition function Zsum()
   - Thermodynamics avearge
   - cylint: use explicit initialisation routines
   - exact test case for the thermodynamice average
   - shifted LJ potential (yet another way to compute the volume)
   - set T from commandline
   BUT: volume is still not correct: Problem: Hard walls at Rmax are not accounted for in the density of states (see Labbook II,p81)

   Revision 1.4  2002/08/23 18:11:53  oliver
   - added sigma, epsilon as globals (the alternativ Lennar-Jones
     representation)
   - made c6, c12 global (can be changed in the debugger or in software)
   - Weeks, Chandler, Anderson split of the LJ potential: repulsive part
     vRlj and attractive part vAlj

   Revision 1.3  2002/08/22 23:18:50  oliver
   - Use new integration routine from numerical recipes (Gauss-Legendre quadrature) for specified number of points. Pore profiles converge now...
   - more options (for profile and number of points for the integral)
     NB: segmentation faults etc happen if options are not given properly

   Revision 1.2  2002/08/22 02:49:22  oliver
   - Restructuring of headers
   - removed mgeom -- I'm not using it anyway
   - pgeom calculates a statistical mechanics estimate for the pore volume

   Revision 1.1  2002/08/21 16:14:04  oliver
   Raw pcode from Numerical Recipes (online)


   Integrate a volume in cylindrical coordinates over a cylinder. This
   is a simplified version of 'Numerical Recipes in C, p164: quad3d'

   I specify my boundaries as numbers not functions because of
   cylindrical symmetry.

*/


#include "cylint.h"
#include "util.h"    /* for debugging level, mesg */
#include "nr.h"      /* for amoeba */

/* GLOBAL in cylint 

   The basic problem is that in the end, the intgrator wants a
   function of three floats (nrfunc). To have some flexibility in the
   coice of function, and especially retaining the capability of
   parameters in these function (like temperature, position of LJ
   centers -> functional of V(r)) I use a fair amount of gloabl
   variables. They are initialised in the top routines that are
   visible from outside.  */
static real **cyatm;   /* array of triples (R_i, PHI_i, Z_i) */
static real **xyzatm;  /* array of triples (x,y,z) */
static int   nratm;   /* number of atoms in atm */
static real sigma;     /* repulsive core of LJ */
static real epsilon;   /* well depth of LJ */
static real rmin;      /* position of the minimum */
/* for distances, only squares are needed */
static real sigma2;     /* repulsive core of LJ */
static real rmin2;      /* position of the minimum */
static real invrmin2;   /* 1.0/rmin2 --fast comparison in vRlj */
static real minvlj=0;   /* global minimum of the full potential
                           (within the integration volume) */
/* for the volume operator: energy of a test particle (determines the
   classically allowed configurations). For SPC water E=6/2 kT --> beta E = 3
*/
static real betaE;
/* LJ parameters (in kT/nm^-xx) */
static real c12;     /* Lennard-Jones potential parameters */
static real  c6;


/* for the thermodynamic average: template functions */
static real (*UU)(real,real,real);        /* potential in kT */
static real (*AA)(real,real,real,real);  /* observable      */

/* UU, AA need to be initialised */
void init_potential  (real (*U)(real,real,real));
void init_observable (real (*A)(real,real,real,real));

static real   exp_bH (real,real,real);
static real A_exp_bH (real,real,real);
static real   gE (real,real,real);    /* density of states x volume operator */

/* gE fails if the potential was not properly shifted which happens
   all the time because the minimum finder cannot guarantee to find
   the global min. Thus, if we fail we improve our potential and start
   over again. */
static Boolean gE_success = TRUE;
static real new_v0[4] = {0,0,0,1e10};     /* x,y,z,V */


/******************************************************************
 * 3D quadrature from Numerical Recipes in C 
 */
static real f1(real);
static real f2(real);
static real f3(real);
static real qgaus10(real (*)(real), real, real);
static real qngaus(real (*)(real), real, real);
static void gausslegendre(real, real, real [], real [], int);

static real xsav,ysav;
static real (*nrfunc)(real,real,real);

/* abscissas and weights for qngauss (n point Gauss quadrature)*/
static int qnorder;
static real *x=NULL; 
static real *w=NULL; 

/* GLOBAL boundaries (in the original 3D integrator these were functions) */
static real ylimit1, ylimit2;
static real zlimit1, zlimit2;


/*************************************************************
 * find the minimum (uses simplex search (amoeba) from NR)
 */
void display_simplex (int,real **,real *);


/* 
   thermodynamic average <A> = Tr[A exp(-bH)]/Tr[exp(-bH)] 
   
   Here the Hamiltonian never depends on the momenta explicitly so
   this bit could be done analytically (ideal gas) but as it cancels
   in the average anyway I dont bother and calculate the
   configurational parts only:

   <A> = Tr[A exp(-bU)]/Tr[exp(-bU)] 

*/
void init_potential (real (*U)(real,real,real)) {
  UU = U;
  return;
}

void init_observable (real (*A)(real,real,real,real)) {
  AA = A;
  return;
}

/* canonical partition function (Zustandssumme) */
real Zsum (real (*U)(real,real,real), real Rmax, real zmin, real zmax) {
  init_potential(U);
  return cylint(exp_bH, 0,Rmax, 0,2*PI, zmin,zmax);
};

real tdavg (real (*A)(real,real,real,real), real (*U)(real,real,real), 
	     real Rmax, real zmin, real zmax) {
  /* patch A and U into the templates in A_exp_bH and exp_bH 
     U must return in units of kT !
   
     this dirtiness is called for to keep the final functions adhere
     to (*func)(real,real,real) form 
  */
  init_potential(U);   /* both are needed for A_exp_bH */
  init_observable(A);
  return cylint(A_exp_bH, 0,Rmax, 0,2*PI, zmin,zmax)/cylint(exp_bH, 0,Rmax, 0,2*PI, zmin,zmax);
};


/* Volume operator with the density of states (which is taken for a
   free particle in a cylindrical box, i.e. sqrt(E) 
   see Labbook II,p84
*/
#define TWOSQRTPI 1.1283791670955125   /* 2/sqrt(PI)  */
real gE (real r, real phi, real z) {
  /* GLOBAL: func UU 
             bool gE_success
	     real new_v0[4]
   */
  real V, sqrtV;
  V = UU(r,phi,z);
  sqrtV=sqrt(V);
  if (V<0) {
    /* flag failure for post mortem find min again */
    gE_success = FALSE;
    /* store new min if lower than last */
    if (V<new_v0[3]) {
      new_v0[XX]=r*cos(phi);
      new_v0[YY]=r*sin(phi);
      new_v0[ZZ]=z;
      new_v0[3]=V;
    }
#ifdef DEBUG
    printf("gE: V(r=%-f,phi=%-f,z=%-f)=%f\n",r,phi,z,V);
    printf("gE: V(x=%-f,y  =%-f,z=%-f)=%f\n",r*cos(phi),r*sin(phi),z,V);
#endif
  }
  return r * (TWOSQRTPI*sqrtV*exp(-V) + erfc(sqrtV));
}

#define XV_ITER_MAX 20
static int xV_iterations = 0;

/* accessible volume */
real xV (real (*U)(real,real,real), real Rmax, real zmin, real zmax) {
  /* GLOBAL: gE_success
             new_v0
  */
  real vol,min;
  init_potential(U);
  gE_success = TRUE;
  vol = cylint(gE, 0,Rmax, 0,2*PI, zmin,zmax);
  if (gE_success) {
    xV_iterations = 0;
    return vol;
  }

  /* dirty ... should be handled in main */
  min=find_min(vljvec,Rmax,zmin,zmax,new_v0);
  if (min > new_v0[3]) { /* should not happen but does it all the time :( */
    mesg(INPUT,"xV(): find_min=%f > new_v0=%f, so the last one is used.",min,new_v0[3]);
    min = new_v0[3];
  }
  /* shift original potential even more 
     (dirty to manipulate the global minvlj here) */
  minvlj += min;
  if (xV_iterations++ >= XV_ITER_MAX)  return 0; /* "Last Exit For The
                                                    Lost" (Fields of
                                                    the Nephilim) */
  
  mesg(VERBOSE,"Restarted xV() with new min=%f ... [%d]",minvlj,xV_iterations); 
  return xV(U,Rmax,zmin,zmax);
};


real exp_bH (real r, real phi, real z) {
  /* GLOBAL: UU */
  return r * exp(-UU(r,phi,z));
} 

real A_exp_bH (real r, real phi, real z) {
  /* GLOBAL: UU 
             AA
   */
  real u;
  u = UU(r,phi,z);  
  return r * AA(r,phi,z,u)*exp(-u);
} 

real Qvol (real r, real phi, real z,real u) {
  /* GLOBAL: betaE      energy of the test particle (in kT)   */
#ifdef DEBUG2
  if (u+1 >= betaE) printf ("Qvol, u+1 >= betaE: @(%f %f %f) u=%f betaE=%f\n",r,phi,z,u,betaE);
#endif
  return (betaE >= u) ? 1.0 : 0.0;
}


/* initialise the global variables needed in vljcyl */
void init_vljcyl (struct potpars *pp) {
  nratm=pp->ncenters;
  cyatm=pp->centers; 
  xyzatm=pp->xyzcenters;
  c6=pp->c6/(kBoltzmann*pp->Temp);      /* unit is now  kT nm^-6  */
  c12=pp->c12/(kBoltzmann*pp->Temp);    /* unit is now  kT nm^-12  */
  /* for <Q> */
  betaE=pp->Nfreedom/2.0;               /* energy of a test particle */
  /* for the splitted potentials */
  sigma =pow(c12/c6,1.0/6.0);        /* repulsive core radius (in nm) */
  sigma2=pow(c12/c6,1.0/3.0);        /* for distances, only squares are needed */
  rmin  =pow(2*c12/c6,1.0/6.0);      /* minimum of the LJ potential */
  rmin2 =pow(2*c12/c6,1.0/3.0);     
  invrmin2=1.0/rmin2; 
  epsilon=0.25*c6*c6/c12;
  minvlj=pp->min;                    /* from find_min(), kT */
  return;
}


/* 
   configurational part Z_V of the partition function (same as above,
   but proper name):

   Tr exp(-beta H) = Z_p Z_V
   Z_V = Tr exp(-beta V_LJ) 
*/
real Z_V (real r, real phi, real z) {
  real s;   /* 1/s = (X-R_i)^2 (distance between integration point
		X(r,phi,z) and center of a LJ atom squared) */
  real u=0; /* cumulative value of Sum_i U(|X-R_i|) */
  int i;

  for(i=0;i<nratm;i++) {
    s=1.0/(r*r + cyatm[i][RAD]*cyatm[i][RAD] - 2*r*cyatm[i][RAD]*cos(phi-cyatm[i][PHI])
      + (z-cyatm[i][ZZZ])*(z-cyatm[i][ZZZ]));
    s = s*s*s;                    /* r^-6               */
    u += c12*s*s - c6*s;          /* c12/r^12 - c6/r^6  */
  }
  return r*exp(-u);
}

/* 
   trace of the volume operator Q(r,E) = Theta(E-V_LJ(r))
*/
real Q (real r, real phi, real z) {
  /* GLOBAL: betaE      energy of the test particle (in kT) */
  real s;   /* 1/s = (X-R_i)^2 */
  real u=0; /* cumulative value of Sum_i U(|X-R_i|) */
  int i;

  for(i=0;i<nratm;i++) {
    s=1.0/(r*r + cyatm[i][RAD]*cyatm[i][RAD] - 2*r*cyatm[i][RAD]*cos(phi-cyatm[i][PHI])
      + (z-cyatm[i][ZZZ])*(z-cyatm[i][ZZZ]));
    s = s*s*s;                    /* r^-6                     */
    u += c12*s*s - c6*s;          /* c12/r^12 - c6/r^6  in kT */
  }
  return (betaE >= u) ? r*exp(-u) : 0;
}

real Zero (real r,real phi,real z) {
  return 0.0;
}

real LinCheck (real r, real phi, real z) {
  /* parameters chosen so that f(r=0.5)==3 
     f[r] = 6nm^-1 r
   */
  return 6*r;
}

/*******************************************************
 * 
 * Lennard-Jones 12-6 potential, shifted by epsilon, and WCA split
 *
 * for atomic centres R_i (in cylindrical
 * coordinates). c6 and c12 are in units of kT
 *
 * V(x=(r,phi,z)) = Sum_i c12/(x-R_i))^12 - c6/(x-R_i)^6
 *
 */
real vljcyl (real r, real phi, real z) {
  /* GLOBAL: nratm       number of LJ centers
             cyatm         centers (cylindrical coordinates)
  */
  real s;   /* 1/s = (X-R_i)^2 (distance between integration point
		X(r,phi,z) and center of a LJ atom squared) */
  real u=0; /* cumulative value of Sum_i U(|X-R_i|) */
  int i;

  for(i=0;i<nratm;i++) {
    s=1.0/(r*r + cyatm[i][RAD]*cyatm[i][RAD] - 2*r*cyatm[i][RAD]*cos(phi-cyatm[i][PHI])
      + (z-cyatm[i][ZZZ])*(z-cyatm[i][ZZZ]));
    s = s*s*s;                    /* r^-6               */
    u += c12*s*s - c6*s;          /* c12/r^12 - c6/r^6  */
  }
  return u;
}

/* 
   LJ with a vector as argument (needed for amoeba)
   still cylindrical coordinates

   shifted
*/
real vljvec (rvec x) {
  /* GLOBAL: nratm       number of LJ centers
             xyzatm         centers (CARTESIAN coordinates)
  */
  real s;   /* 1/s = (X-R_i)^2 (distance between integration point */
  real u=0; /* cumulative value of Sum_i U(|X-R_i|) */
  int i;

  for(i=0;i<nratm;i++) {
    s=1.0/((x[XX]-xyzatm[i][XX])*(x[XX]-xyzatm[i][XX])
	   + (x[YY]-xyzatm[i][YY])*(x[YY]-xyzatm[i][YY])
	   + (x[ZZ]-xyzatm[i][ZZ])*(x[ZZ]-xyzatm[i][ZZ]));
    s = s*s*s;                    /* r^-6               */
    u += c12*s*s - c6*s;          /* c12/r^12 - c6/r^6  */
  }
  return u - minvlj;
}


/* shifted LJ potential:
   V(r) -> V(r) + epsilon

   Copied&pasted to save another function call (instead of
   return vljcyl() + epsilon )

*/
real vljshiftcyl (real r, real phi, real z) {
  /* GLOBAL: nratm       number of LJ centers
             cyatm         centers (cylindrical coordinates)
	     minvlj        lower bound of potential, in kT
	                   (needs to be found first and initialised!)
  */
  real s;   /* 1/s = (X-R_i)^2 (distance between integration point
		X(r,phi,z) and center of a LJ atom squared) */
  real u=0; /* cumulative value of Sum_i U(|X-R_i|) */
  int i;

  for(i=0;i<nratm;i++) {
    s=1.0/(r*r + cyatm[i][RAD]*cyatm[i][RAD] - 2*r*cyatm[i][RAD]*cos(phi-cyatm[i][PHI])
      + (z-cyatm[i][ZZZ])*(z-cyatm[i][ZZZ]));
    s = s*s*s;                    /* r^-6               */
    u += c12*s*s - c6*s;          /* c12/r^12 - c6/r^6  */
  }
  return u - minvlj;
}


/***************************
 * Splitting the LJ potential according to Weeks, Chandler & Andersen (1971) 
 *
 * WCA LJ:    repulsive part 
 */
real vRljcyl (real r, real phi, real z) {
  /* GLOBAL: nratm       number of LJ centers
             cyatm         centers (cylindrical coordinates)
  */
  real s;   /* 1/s = (X-R_i)^2 (distance between integration point
		X(r,phi,z) and center of a LJ atom squared) */
  real u=0; /* Sum_i U(|X-R_i|) */
  int i;

  for(i=0;i<nratm;i++) {
    s=1/(r*r + cyatm[i][RAD]*cyatm[i][RAD] - 2*r*cyatm[i][RAD]*cos(phi-cyatm[i][PHI])
	      + (z-cyatm[i][ZZZ])*(z-cyatm[i][ZZZ]));
    if (s > invrmin2) {              /* r < rmin: soft repulsive (faster: s=1/r^2 > 1/rmin^2)*/
      s = s*s*s;                     /* r^-6                     */
      u += c12*s*s - c6*s + epsilon; /* vLJ(r) + epsilon         */
    };                               /* r>=rmin: vR = 0          */
  }
#ifdef DEBUG2
  printf ("vRljcyl(r=%f,phi=%f,z=%f)= %f  exp:%f\n",r,phi,z, 
	  r*exp(-u), -u);
#endif
  return r*exp(-u);
}

/***************************
 * WCA LJ: attractive part
 */

real vAljcyl (real r, real phi, real z) {
  /* GLOBAL: nratm       number of LJ centers
             cyatm         centers (cylindrical coordinates)
  */
  real s;   /* 1/s = (X-R_i)^2 (distance between integration point
		X(r,phi,z) and center of a LJ atom squared) */
  real u=0; /*  Sum_i U(|X-R_i|) */
  int i;

  for(i=0;i<nratm;i++) {
    s=sigma2/(r*r + cyatm[i][RAD]*cyatm[i][RAD] - 2*r*cyatm[i][RAD]*cos(phi-cyatm[i][PHI])
      + (z-cyatm[i][ZZZ])*(z-cyatm[i][ZZZ]));
    if (s > invrmin2) {              /* r < rmin: (faster: s=1/r^2 > 1/rmin^2)*/
      u += -epsilon;                 /* const attractive         */
    } else {                         /* attractive tail          */
      s = s*s*s;                     /* r^-6                     */
      u += c12*s*s - c6*s;           /* vLJ(r)                   */
    };              
  }
#ifdef DEBUG2
  printf ("vAljcyl(r=%f,phi=%f,z=%f)= %f  exp:%f\n",r,phi,z, 
	  r*exp(-u), -u);
#endif
  return r*exp(-u);
}

/* determine the minimum (in cylindrical coordinates) with the
   down-hill simplex method 

   NB: currently this is probably only good to find one minimum in the
   whole pore volume. The potential has multiple minima, and I only
   assume that any one is good enogh (there's also a flat maximum with
   vanishing derivatives at the center which would probably give
   gradient methods some grief. 

   Anyway, I have a rough idea where to start looking.
*/

#define FIND_MIN_ITER_MAX   10
static int find_min_iterations = 0;

real find_min (real (*f)(rvec), real Rmax, real zmin, real zmax, rvec v0) {
  /* GLOBAL: int  nratm  atoms in cartesian coordinates
             real **cyatm  -"-
	     real sigma    repulsive core of a single potential
  */
  real **simplex;     /* array of 4 vertices */
  real *fsplx;        /* value of the function at each vertex */
  real *a1,*a2;       /* atoms */
  real min=0;
  real Rfac;          /* typical length: Rmax * Rfac */
  real phi0,r0;       /* starting angle and radius */
  int nevals;         /* number of evaluation steps in amoeba */
  int i;

  simplex=matrix(1,DIM+1,1,DIM);      /* each row is a vertex of the simplex  */
  fsplx=vector(1,DIM+1);              /* function values on the DIM+1
				      corners of the starting simplex */

  /* NB: amoeba wants a function of a vector -->use vljvec(rvec) */
  /* setup the first simplex:
   * use initial starting point v0 if non-null
     1: center
     2: right by 30degrees and down 3/4 of (z2-z1)
     3: left  by 30degrees and down 3/4 of (z2-z1)
     4: up 3/4 and stay in the minimum-zone
   */
#define R_FUDGE 0.35    /* the minimum is typically somewhere there */
  
  if (v0) {
    Rfac=0.1;             /* keep close to the starting point */
    simplex[1][1]=v0[XX];
    simplex[1][2]=v0[YY];
    simplex[1][3]=v0[ZZ];
  }
  else {
    Rfac=R_FUDGE;
    v0=vector(0,2);
    v0[XX]=simplex[1][1]=0.1*Rmax;
    v0[YY]=simplex[1][2]=0.1*Rmax;
    v0[ZZ]=simplex[1][3]=0.5*(zmax+zmin);
  }
  
  simplex[2][1]=v0[XX] + Rfac*Rmax;
  simplex[2][2]=v0[YY];
  simplex[2][3]=v0[ZZ];

  simplex[3][1]=v0[XX];
  simplex[3][2]=v0[YY] + Rfac*Rmax;
  simplex[3][3]=v0[ZZ];

  simplex[4][1]=v0[XX];
  simplex[4][2]=v0[YY];
  simplex[4][3]=v0[ZZ] + Rfac*Rmax;


  for(i=1;i<=DIM+1;i++) 
    /* uff -- simplex is [1..3] but f() expects x[0..2] 
       note: simplex[1]+1 is a _pointer_ to element 1  */
    fsplx[i]=f(simplex[i]+1);

  if (debuglevel >= SUB1) display_simplex(DIM+1,simplex,fsplx);
 
  if (!amoeba(simplex, fsplx, DIM, 0.01*FTOL, f, &nevals)) {
    mesg(WARN,"Minimisation failed after %d steps--trying again [%d]...",
	 nevals,find_min_iterations);
    if(debuglevel>=SUB1) display_simplex(DIM+1,simplex,fsplx);

    if (find_min_iterations++ < FIND_MIN_ITER_MAX) {
      phi0=2*PI * rand()/(real)RAND_MAX;
      r0=Rmax * rand()/(real)RAND_MAX;
      v0[XX]=r0 * cos(phi0);
      v0[YY]=r0 * sin(phi0);
      v0[ZZ]=rand()/(real)RAND_MAX * (zmin+zmax);
      mesg(SUB1,"Restarting with randomised initial vertex (%f,%f,%f)\n",
	   v0[XX],v0[YY],v0[ZZ]);
      min = find_min(vljvec,Rmax,zmin,zmax,v0);
    }
    else 
      mesg(WARN,"find_min(): Could not find a minimum---continue with what we have.");
    find_min_iterations = 0;  
  } 
  else {
    min = fsplx[1];

    mesg(VERBOSE,"Minimisation returned after %d steps.",nevals);
    if (debuglevel >= SUB1)   display_simplex(DIM+1,simplex,fsplx);
  }

  free_matrix(simplex,1,DIM+1,1,DIM);
  free_vector(fsplx,1,DIM+1);
  return min;
}

void display_simplex (int n,real **s,real *y) {
  int i;
  printf("%6s%8s%8s%8s %8s\n","vertex","X","Y","Z","f(X,Y,Z)");
  for (i=1;i<=n;i++) {
    printf("%6d%8.5f%8.5f%8.5f %8.5f\n",i,s[i][1],s[i][2],s[i][3],y[i]);
  }
}

/*******************************************************
 * 
 * Integration routines
 * source: Numerical Recipes in C
 *
 */

/* initialise the abscissas and weights */
void init_gaussleg (int n) {
  qnorder = n;
  if (!(x=(real *)calloc(n+1,sizeof(real))) ||
      !(w=(real *)calloc(n+1,sizeof(real))))
    fatal_error(0,"Failed to allocate memory in init_gaussleg().");
  gausslegendre(0,1,x,w,n);
  return;
}

void free_gaussleg (void) {
  if(x) { 
    free(x); x=NULL; 
  }
  if(w) {
    free(w); w=NULL;
  }
  return;
}

/* Returns the integral of a user-supplied function func over a
   three-dimensional region specified by the limits x1, x2, and by the
   user-supplied functions ylimit1, ylimit2, z1, and z2, as defined in
   (4.6.2). (The functions y1 and y2 are here called ylimit1 and ylimit2 to
   avoid conflict with the names of Bessel functions in some C
   libraries). Integration is performed by calling qgaus recursively.  

   (which I find difficult to follow.. recursion makes my head going gooey.)
*/

real cylint(real (*func)(real, real, real), 
	     real r1, real r2, real phi1, real phi2,
	     real z1, real z2){ 
  /* Don't want to do init_gaussleg here; otherwise I had to free it
     at the end and that means recomputing weights for every
     call. Better to do it once upstairs. Trust me, I know what I'm doing. 
  */
  if (!w || !x) fatal_error(1,"Uh-oh... one should have run init_gaussleg()!");

  /* initialise private vars */
  ylimit1=phi1; ylimit2=phi2;
  zlimit1=z1;   zlimit2=z2;
  nrfunc=func;
  return qngaus(f1,r1,r2);
}


/*    This is H of eq. (4.6.5).*/
real f1(real x) {    
  xsav=x;
  return qngaus(f2,ylimit1,ylimit2);
}


/* This is G of eq. (4.6.4). */
real f2(real y) {    
  ysav=y;
  return qngaus(f3,zlimit1,zlimit2);
}

/* The integrand f(x, y, z) evaluated at fixed x and y. */
real f3(real z)    
{    
  return (*nrfunc)(xsav,ysav,z);
}


/* Variable Gauss-Legendre quadrature: Initialise the abscissas and
   weights by calling init_gaussleg() ONCE

   The abscissas and weights are presumed to be found in the GLOBAL
   arrays x[] and w[] to avoid cluttering the 3D integration and so is
   the order, qnorder. 

   The arrays are reduced to the interval [0,1[ and are scaled in
   qngauss.

   First value of each array not used.    
 */

real qngaus(real (*func)(real), real a, real b) {    
  /* GLOBAL:
     qnorder            order of the Gauss-Legendre quadrature
                        must be even
     x[1..qnorder]      abscissas in [0;1[
     w[1..qnorder]      weights 
  */
     
  int j;
  real xr,dx,s;
  /* real  xm; */
  
  /* xm=0.5*(b+a); */           /* only needed if I use the symmetry of the weights */
  /* xr=0.5*(b-a); */
  xr=b-a;
  s=0;        
  for (j=1;j<=qnorder;j++) {    /* Integral[f(x),a,b] ~= Sum w_j f(x_j) */
    dx=xr*x[j];               
    s += w[j]*func(a+dx);
  }
  return s *= xr;               /* Scale the answer to the range of integration.*/
}


/*
 Given the lower and upper limits of integration x1 and x2, and given
 n, this routine returns arrays x[1..n] and w[1..n] of length n,
 containing the abscissas and weights of the Gauss- Legendre n-point
 quadrature formula.  */
void gausslegendre(real x1, real x2, real x[], real w[], int n)
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
    } while (fabs(z-z1) > GAUSSLEG_EPS);
    x[i]=xm-xl*z;           /*  Scale the root to the desired interval, */
    x[n+1-i]=xm+xl*z;       /* and put in its symmetric counterpart. */
    w[i]=2.0*xl/((1.0-z*z)*pp*pp);  /*    Compute the weight */
    w[n+1-i]=w[i];          /*  and its symmetric counterpart. */
  }
}


/* 
   Fixed-abscissa Gauss-Legendre quadrature (10 points):

   Returns the integral of the function func between a and b, by
   ten-point Gauss-Legendre inte- gration: the function is evaluated
   exactly ten times at interior points in the range of integration.  

   OB: I just take this as a first attempt; perhaps this is not
   appropriate in cylindrical coordinates ?  
*/

real qgaus10(real (*func)(real), real a, real b) {    
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
