/*
  $Id: xg.h,v 1.11 2008/02/18 01:18:48 oliver Exp $
  $Log: xg.h,v $
  Revision 1.11  2008/02/18 01:18:48  oliver
  o Can read geometry from input file:
    - added datafile_geom()
    - updated usage
  o Set MAXDOMAINS to 6 and MAXRINGS to 4 (but cannot be much more otherwise binary
    crashes with segfault/buserror). This should be dynamically allocated
    anyway...
  o use fatal_error() and mesg() instead of printf/exit constructs
  o started cleaning up Log and whitespace

  Revision 1.10  2003/05/17 21:25:16  oliver
  pgeom segfaulted so I reduced the dimensions of the arrays considerably and recompiled; also removed spurious EMPTY definition of SQR

  Revision 1.9  2002/11/13 15:04:25  oliver
  increased MAXLAYERS, MAXRINGS, MAXSITES (so models can be generated
  with -x; with bonds etc. the arrays are too small but I cannot be
  bothered to rewrite this with dynamic alloc)

  Revision 1.8  2002/11/07 14:55:15  oliver
  - some volume calculation stuff was not checked in
  - added METHYL (distance 0.14nm) to species
  - rewrote species id using enum

  Revision 1.7  2002/08/27 02:20:51  oliver
  works but contains lots of dirt hacks in cylint.c...

  Revision 1.6  2002/08/25 22:19:52  oliver
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

  Revision 1.4  2002/08/23 18:16:39  oliver
  - potential function for the volume can be set from the commandline
    (and is stored in struct profile.pot)
  - cleanup of header file mess

  Revision 1.3  2002/08/22 23:18:52  oliver
  - Use new integration routine from numerical recipes (Gauss-Legendre quadrature) for specified number of points. Pore profiles converge now...
  - more options (for profile and number of points for the integral)
    NB: segmentation faults etc happen if options are not given properly

  Revision 1.2  2002/08/22 02:49:27  oliver
  - Restructuring of headers
  - removed mgeom -- I'm not using it anyway
  - pgeom calculates a statistical mechanics estimate for the pore volume

  Revision 1.1  2000/12/13 16:40:52  oliver
  split functions in main (mgeom, pgeom), common (xg), and general (pg) -- stupid names...
  compile mgeom or pgeom in Makefile, use -DPROG
  all prototypes in separate files because my includes ar a crappy mess


  header common to all Xgeom programms

*/

#ifndef XG_H_READ
#define XG_H_READ

#include <stdio.h>
#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include "xgtypes.h"
#include "util.h"
#include "mol_data.h"


/* default force constants */
#define GMX_BOND_C1   1.0e+5  /* force constant, typical: 1e+5 */
#define GMX_ANGLE_C1   1.0e+3 /* normally ~ 1e+2 */

/* bonds */
#define DLMAX 4.5             /* max bond length in layer */
#define DZMAX 5.0             /* max bond length parallel pore axis */
#define DLMIN 1.0             /* min bond length in layer */
#define DZMIN 1.0             /* min bond length parallel pore axis */


/* potential stuff (from cylint.h) */
#define CSIX_OW_MTH     0.00587289    /* OW - CH4 Lennard Jones (ffgmx with mixing */
#define CTWELVE_OW_MTH  9.51217e-06   /* rules); units: kJ/mol nm^-6|12 */
#define TEMPERATURE     300       /* Kelvin */
#define NFREEDOM_SPC    6         /* beta E = 6/2  for SPC water (3x3 - 3 degrees of freedom) */

#define N_GAUSSLEG   40           /* default number of abscissas in
                                     the Gauss-Legendre quadrature */

#define POT_PLOT_SLICES 5    /* number of slices in potential plots */


#define MAXDOMAINS   6        /* max number of domains */
#define MAXLAYERS    4        /* max number of layers in domain */
#define MAXRINGS     4        /* max number of rings per layer */              
#define MAXSITES    140       /* max number of sites on a ring. GUESSED. 
				 ? = Pi/(ArcSin[1/(Rpore/rA + 2*MAXRINGS)]) */
#define TOTALSITES   MAXDOMAINS * MAXLAYERS * MAXRINGS * MAXSITES

#define MAXBONDSSITE 8        /* 12 nearest neighbours in fcc */
#define MAXBONDS     (int) MAXBONDSSITE * TOTALSITES / 2
                     /* upper bound: TOTALSITES * (TOTALSITES-1) / 2 */ 
#define MAXANGLES    MAXBONDS  
                     /* really: B*(B-1)/2 * TOTALSITES, B=MAXBONDSSITE */ 
#define MAXBONDCONSTRAINTS 4


#define STRLEN    2048

extern int debuglevel;

/*
  struct of same name:

  domain
  layer
  geom

  are defined differently for mgeom and pgeom

  WHO CARES -- I'm not using mgeom anyway. I make xg.h the header for pgeom (and remov mgeom)

*/

enum domaintypes { MOUTH, PORE, LAYER }; 

/* 
   Model = Domain > Layer > Ring > Site

   The channel model is built from layers. Each layer consists of concentric 
   rings. On each ring the 'pseudo atoms' (=sites)  sit like beads on a wire.
*/

struct layer {
  double rhomin;      /* inner radius */
  double rhomax;      /* outer radius */
  double z;           /* z coordinate of all sites in this layer */
  int    maxring;     /* number of rings is maxring+1 */
  double dr;          /* radial distance between two rings */
  struct domain *domain; /* is part of that domain */
};

struct ring {
  int    maxsite;  /* number of sites on this ring */
  double dphi;     /* polar angle increment between two sites */
  double rho;      /* radial coordinate of species on that ring */
  struct layer *layer; /* belongs to that layer */
  Boolean  bExposed; /* TRUE if these atoms are facing the pore cavity;
                      needed for quick & dirty volume calculations */
};


struct domain {
  /* 
     a domain is a part of the model, characterized by the radii at
     top r2 and bottom r1, length, and z coordinate of bottom and top
     layer.
     this struc collects also derived quantities like q, dz
   */
  enum domaintypes type;
  char    *description;
  double  l;    /* length of domain */
  double  (*rho)(double z, double z1, double z2, double r1, double r2); 
                /* rho(z) shape function which describes the form of
                   the inner pore wall between z1 and z2 */
  struct spec *species;
                /* default species for this domain (can be changed on a 
		   per site basis as well */
  double  r1;   /* inner boundary r(z1) */
  double  r2;   /* inner boundary r(z2) */
  double  r_outer; /* outer radius; usually = geom->r_outer */
                 /* redundant, but convenient: all in one place */
  double  rho1; /* rho(z1) */
  double  rho2; /* rho(z2) ; actual cylindrical coordinate of the center */
  double  rho_outer; 
                /* center of autermost species */ 
  double  z1;   /* z of center of atom at bottom of barrel */
  double  z2;   /* z " at top */
  int      q;   /* number of layers */
  double  dz;   /* interlayer distance */
  int     first_site;  /* site numbers in model[] */
  int     last_site;
  struct  domain *prev;
  struct  domain *next;
};


struct geom {
  char    *coordfile;       /* name of output file */
  char    *topofile;        /* itp topology file */
  double  r_outer;          /* outer radius */
  int     ndomains;         /* number of domains */
  int     nsites;           /* number of sites in model[] */ 
  int     nbonds;           /* number of bonds in bonds[] */
  int     nangles;          /* number of angles in angles[] */
  int     nbc;              /* number of constraints in bc[] */
  struct bond_constraint *bc;
  double k_bond, k_angle;   /* force constants */
  Boolean connectonly;      /* write only connectivity in itp */
  Boolean atomsonly;        /* ignore bonds and angles */
  Boolean shiftcbox;        /* shift coordinates to cavitybox or unitcell */
  struct domain *domain[MAXDOMAINS];  
         /* N terminal domain at domain[0] */
  struct pdb_CRYST1 unitcell; 
  struct pdb_CRYST1 cavitybox; 
};

struct bond_constraint {
  int    serial;
  double dmin;          /* min distance between atoms to form a bond 
			   in principle, this should be vdW radius 
			   but in practice I use it to tune bonds
			*/
  double dmax;          /* max distance between atoms to form a bond */
  double gamma_min;     /* bonds orientation with respect to pore axis */   
  double gamma_max;     /* min <= gamma <= max: bond  */
};


struct std_input {
  /* sorted command line arguments */
  double r_pore;
  double l_pore;
  double r_mouth;    
  double l_mouth;
  double r_outer;
  int    specid;     /* id of species */
  char   *coordfile; /* pdb */
  char   *topofile;  /* GROMACS itp */
  int    n_bc;       /* number of constraints in bc[] */
  struct bond_constraint bc[MAXBONDCONSTRAINTS];
  double k_bond, k_angle;   /* force constants */
  Boolean connectonly;  
  Boolean atomsonly;
  Boolean shiftcbox;
  char   *datafile;  /* set if reading geometry specs from datafile */
};


struct pprofile {
  char  fn[STRLEN]; /* filename for output */
  Boolean bSet;     /* from commandline: calculate profile ? */
  Boolean bPlot;    /* plot 2D potential slices ? */
  int   nzplot;     /* how many ? */
  struct potpars *pp; 
  int   nz;         /* number of slices along z */
  real Rmax;       /* integrate out to Rmax */
  real z1, z2;     /* within z1 < z < z2 */
  real **r;        /* the profile: r[0] = z, r[1] = Reff */
};



/* hard coded list of species */

enum { METHANE, CARBON, OXYGEN, METHYL, NSPEC };

extern struct spec species [];


#define NZPROF 50  /* number of slices in a profile */

extern void print_species (void);

extern Cartesian atom_diff (struct pdb_ATOM *, struct pdb_ATOM *);
extern struct pdb_ATOM translate_pdb_atom (struct pdb_ATOM *, Cartesian);
extern struct itp_angle new_angle (struct pdb_ATOM *, struct pdb_ATOM *, 
			    struct pdb_ATOM *, double);
extern void center (struct pdb_ATOM [], struct geom *);
extern void center_model (Cartesian, struct pdb_ATOM [], int);
extern int write_pdb (struct pdb_ATOM *, struct itp_bond *, struct geom *);
extern int pdb_write_header (FILE *, struct geom *); 
                     /* different function but same prototype */
extern int pdb_write_cavitybox (FILE *, struct pdb_CRYST1 *);

extern int write_topology (struct pdb_ATOM *, struct itp_bond *, struct itp_angle *,
		    struct geom *);

extern int write_bonds (FILE *, struct itp_bond *, struct geom *);
extern int write_constraints (FILE *, struct itp_bond *, struct geom *);
extern int write_atoms (FILE *, struct pdb_ATOM *, struct geom *);
extern int write_angles (FILE *, struct itp_angle *, struct geom *);
extern int write_moleculetype (FILE *);
extern double find_dmax (struct bond_constraint *, int);
extern int do_bonds (struct itp_bond *, struct pdb_ATOM *, struct geom *);
extern struct bond_constraint *bond_class (Cartesian, struct bond_constraint [], int);
extern int do_angles (struct itp_angle *, struct itp_bond *, struct geom *);

/* from pgeom.h */
extern double constant (double, double, double, double, double);
extern int default_geom (struct domain []);
extern int do_coordinates (struct pdb_ATOM *, struct geom *);


extern int do_layer (struct domain *, struct pdb_ATOM *);
extern int do_ring (struct layer *, struct pdb_ATOM *);
extern int do_site (struct ring *, struct pdb_ATOM *);
extern int input_geom(struct std_input *, struct geom *);
extern struct layer setup_layer (double, struct domain *);
extern double linear(double, double, double, double, double);

extern struct ring setup_ring (double, struct layer *);
extern void setup_domain (struct geom *); 
extern struct domain ini_domain (enum domaintypes, char *, double,
			  double (*)(), struct spec *, double, double,
			  struct domain *, struct domain *);
extern struct std_input default_input (void);
extern int standard_geom (struct std_input *, struct domain []);
extern void unitcell (struct geom *); 
extern void cavitybox (struct geom *);



extern void print_domain(struct domain *);
extern void print_geom(struct geom *);
extern void print_layer (struct layer *);
extern void print_ring (struct ring *);
extern void print_site (struct pdb_ATOM *);

extern void print_usage (char *, char *);    

extern real volume (struct potpars *, struct pdb_ATOM [], struct domain *, struct pprofile *);
extern void calc_profile (struct pprofile *);
extern void plot_potential (real (*)(real,real,real),real,real,real,int);

#endif /* XG_H_READ */
