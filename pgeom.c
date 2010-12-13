/* 
   $Id: pgeom.c,v 1.37 2008/02/18 01:18:48 oliver Exp $

   Ouline of the program:
   ======================

   Calculate positions in atomistic model of a transmembrane pore.

   INPUT: 
   [all lengths in Angstroem = 10^(-10)m,
    all distances are 'bounding-box' (include all pseudoatoms)]

   Router   outer model radius

   For each domain of the model (cylindrical symmetry):

     Rpore    inner pore radius    
     L        total length of  domain
     spec     species (determines rA)
     rho(z)   shapefunction, describing shape of inner wall
  
  
   OUTPUT:
   list of atomic coordinates in pdb format
   GROMACS topology file (itp)
*/

#include <stdio.h>
#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include "xgtypes.h"
#include "xg.h"
#include "util.h"
#include "mol_data.h"
#include "cylint.h"
#include "nr.h"

/*
   functions

   pore shape functions
   given: z(rho), z
   need:  z^{-1}(rho) === rho(z)
*/

#ifdef DEBUG
#ifndef TESTCASE
#define TESTCASE
#endif
#endif

/* internal test for integration routine cylint() */
#ifdef TESTCASE
real cylconst (real r, real phi, real z) {
  return r;
}
real cylcos (real r, real phi, real z) {
  return r*cos(phi);
}
real cyl3 (real r, real phi, real z) {
  return exp(-r)*cos(phi)*cos(phi)*z*z;
}
#define K1 2.333     /* kT nm^-2 */
real harmonic (real r, real phi, real z) {
  return K1/2.0 * r*r;
}  
#endif /* TESTCASE */

double constant (double z, double z1, double z2, double r1, double r2)
{  
  assert(r1 == r2);
  return r1;
};

double linear (double z, double z1, double z2, double r1, double r2)
{
  double m, t;
  
  assert (z1 != z2);
  assert (r1 != r2);

  m = (z2 - z1)/(r2 - r1);
  t = z1 - m*r1;

  return (z - t)/m;
};

struct domain ini_domain (enum domaintypes t, char * desc, double len,
			  double (*shapefunction)(), struct spec *sp,
			  double r1, double r2,
			  struct domain *prev, struct domain *next)
     /* these parameters have to be given by the user. All other
      parameters are calculated by the programme (eventually...)  
     */
{
  double rA;
  struct domain this;

  rA = sp->radius;
  
  this.type        = t;
  this.description = desc;
  this.l           = len;
  this.rho         = shapefunction;
  this.species     = sp;
  this.r1          = r1;
  this.r2          = r2;
  this.rho1        = r1 + rA;
  this.rho2        = r2 + rA;
  this.prev        = prev;
  this.next        = next;

  /* r_outer is initialized with g->r_outer later on. It is not really
     supossed to be variable on a per domain basis, but conceptually
     it belongs here */
  this.r_outer = 0;
  this.rho_outer = 0;

  /* initialize unknowns with zero values */
  this.z1 = 0.0;
  this.z2 = 0.0;
  this.q  = 0;
  this.dz = 0.0;
  this.first_site = 0;
  this.last_site  = 0;

  return this;
}


int standard_geom (struct std_input *in, struct domain domain_list[MAXDOMAINS])
{
  /* setup standard geometry: 3 domains M-P-M, length, radii */
  #define NDOMAINS_STD 3
  
  assert(NDOMAINS_STD <= MAXDOMAINS + 1);

  domain_list[0] = ini_domain(MOUTH, "intracellular mouth region", in->l_mouth,
			      &linear, &species[in->specid],    
			      in->r_mouth, in->r_pore,
			      NULL, &domain_list[1]);
  domain_list[1] = ini_domain(PORE, "transmembrane pore region", in->l_pore,
			      &constant, &species[in->specid], 
			      in->r_pore, in->r_pore,
			      &domain_list[0], &domain_list[2]);
  domain_list[2] = ini_domain(MOUTH, "extracellular mouth region", in->l_mouth,
			      &linear, &species[in->specid], 
			      in->r_pore, in->r_mouth,
			      &domain_list[1], NULL);
  return NDOMAINS_STD;
};

int datafile_geom (struct std_input *in, struct domain domain_list[MAXDOMAINS])
{
  /* read geometry from file:
     # comment (skipped)
     RADIUS r_outer
     MOUTH r_upper r_lower length [species]
     PORE  r_upper r_lower length [species]
     ...
  */
  FILE *geometry;
  char line[STRLEN], buf[STRLEN];
  int ndomains = 0;
  int nmatch, specid, idom;
  float r_upper, r_lower, length; 
  char *type = buf;
  struct domain *dom;
  enum domaintypes dtype;
  double (*profile_function)();

  geometry=fopen(in->datafile,"r");
  if (! geometry) {
    fatal_error(1, "Failed to open datafile %s.\n",in->datafile);
  };
  while(fgets(line,STRLEN,geometry)) {
    if (line[0] == '#') continue;
    nmatch = sscanf(line,"%s %f %f %f %d",type,&r_upper,&r_lower,&length,&specid);
    if (nmatch==2) {  /* set variables; r_upper always holds the value */
      if (!strcmp("RADIUS",type)) in->r_outer = r_upper;
      mesg(INPUT,"parameter: %s=%f",type,r_upper);
      continue;
    } else if (nmatch==3) specid=in->specid;
    mesg(INPUT,"domain %2d: %s r_u=%f r_l=%f  l=%f specid=%d",
	   ndomains,type,r_upper,r_lower,length,specid);    
    /* record type */
    if (!strcmp("MOUTH",type)) dtype = MOUTH;
    else if (!strcmp("PORE",type)) dtype = PORE;
    else if (!strcmp("LAYER",type)) dtype = LAYER;
    else {
      fatal_error(1,"Unknown domain type %s.\n",type);
    }
    if (ndomains-1 == MAXDOMAINS) {
      fatal_error(1,"Internal limit: only %d domains allowed.\n",MAXDOMAINS);
    }
    profile_function =  (r_lower == r_upper) ? &constant : &linear;
    domain_list[ndomains] = ini_domain(dtype, "custom domain", length,
				       profile_function, &species[specid], 
				       r_upper, r_lower,
				       NULL,NULL);
    ndomains++;
  }
  /* fix domain linkage (already initialized to NULL) */
  if (ndomains > 1) {
    domain_list[0].next = &domain_list[1];
    for(idom=1;idom<ndomains-1;idom++) {
      domain_list[idom].prev = &domain_list[idom-1];
      domain_list[idom].next = &domain_list[idom+1];
    }
    domain_list[ndomains-1].prev = &domain_list[ndomains-2];
    domain_list[ndomains-1].next = NULL;
  }
  return ndomains;
}

struct std_input default_input ()
{
  struct std_input this = {
    3.5,
    8.0,
    7.5,
    4.0,
    15.5,
    METHANE,
    "pore.pdb",
    "pore.itp",
    2,
    {
      {0, DZMIN, DZMAX,  0, 44},
      {1, DLMIN, DLMAX, 45, 90}
    },
    GMX_BOND_C1,
    GMX_ANGLE_C1,
    FALSE,
    FALSE,
    FALSE,
    NULL       /* datafile */
  };
  return this;
};
    


int input_geom(struct std_input *in, struct geom *g) {
  /* input values for the overall geometry of the pore model */

  int i;

  static struct domain all_domains[MAXDOMAINS];

  if (in->datafile) {
    /* read from data file */
    g->ndomains = datafile_geom(in, all_domains);
  } else { /* simple input from commandline */    
    g->ndomains = standard_geom(in, all_domains);
  }
  g->coordfile = in->coordfile;
  g->topofile  = in->topofile;
  g->r_outer   = in->r_outer;    /* can be set from datafile */
  g->nbc       = in->n_bc;
  g->bc        = in->bc;
  g->k_bond    = in->k_bond;
  g->k_angle   = in->k_angle;
  g->connectonly = in->connectonly;
  g->atomsonly   = in->atomsonly;
  g->shiftcbox   = in->shiftcbox;
  for (i = 0; i < g->ndomains && i < MAXDOMAINS; i++)
    {
      g->domain[i] = &all_domains[i];

      /* set r_outer also in each domain. Redundant, but more
         convenient because then all domain information is in one
         place 

	 Kludge...
      */
      g->domain[i]->r_outer   = g->r_outer;
      g->domain[i]->rho_outer = g->r_outer - g->domain[i]->species->radius;
    };
  
  /* initialize to 0 */
  g->nsites = 0;
  g->nbonds = 0;

  /* unitcell struc is filled by unitcell() after setup_domain()*/

  return g->ndomains;
};

void setup_domain (struct geom *g) {
  /* - double pass setup; PORE regios determine global 
       inter-layer spacing dz 
  */
  
  int i;
  double l_pore=0;  /* total length of pore region, center-center */
  int    q_pore;   /* number of gaps between layers */
  int    q_check, nborders, lastborder;
  double dz, dA;
  double z1, z2;
  enum domaintypes lastdtype;
  struct domain *dom;

  dA = 0.0;                /* distance between two atoms ('bond length') */
  q_check = 0;             /* internal topology check */
  nborders = 0;            /* number of interfaces between different domains */
  lastdtype = g->domain[0]->type;

  /* first pass:
     total length of pore region and global dz
     also check correct topology: MOUTH -> PORE -> MOUTH
  */
  for (i = 0; i < g->ndomains; i++) {
    dom = g->domain[i];

    if (dom->type != lastdtype) {
      nborders++;           /* increase nborders for each change M<->P */
      lastborder = i - 1;   /* last PORE region before second MOUTH */
    };
    lastdtype   = dom->type;

    if ( dom->type == PORE) {
      l_pore += dom->l;
      dA      = max(dA, 2 * dom->species->radius);
      mesg (SUB2, "setup_domain(): l_pore = %6.3f  2*rA = %4.2f", l_pore, dA);
    };
  };

  /* if nborders != 2 => error!! */
  if ( nborders != 2 ) {
    fatal_error (2,
    "setup_domain (): Fatal error. Topology of model is incorrect. It should\n"
    "                 have been MOUTH->PORE->MOUTH with 2 borders, but there\n"
    "                 were %d borders in %d domains", nborders, g->ndomains);
  };

  q_pore = (int) l_pore/dA;  
  if ( q_pore < 1 ) {
    fatal_error (2,
		 "setup_domain (): Fatal error. Number of layers in the PORE is not positive\n"
		 "                 (q_pore = %d). (l_pore=%4.2f) < (2*rA=%4.2f)\n",
		 q_pore, l_pore, dA);
  };

  dz = l_pore/q_pore;
  /* second pass:
     number of layers q per domain and z coordinates of lower and upper layer
  */
  z2 = -dz;

  for (i = 0; i < g->ndomains; i++) {
    dom = g->domain[i];
    if ( dom->type == PORE) 
      {
	dom->q   = round (dom->l / dz);
	q_check += dom->q;
	dom->dz  = dz;           /* rather pointless to store dz separately 
				  but I still used the old strucs */
      } 
    else 
      {
	dom->q   = round (dom->l / dA);
	dom->dz  = dz;           /* pointless ... */
      };

    if (dom->q < 0) {
      fatal_error (2,
		   "setup_domain (): Fatal error. Number of layers in domain %d "
		   "is not positive\n"
		   "                 (q = %d). Probably (l=%4.2f) < (2*rA=%4.2f)\n",
		   i,dom->q, dom->l, dA);
    };

    z1 = z2 + dz;
    z2 = z1 + (dom->q - 1) * dz;

    dom->z1 = z1;
    dom->z2 = z2;

  };  

  assert (q_check == q_pore);
  return;
};

int do_coordinates (struct pdb_ATOM *model, struct geom *g)
{
  int j;

  /* loop over domains */
  for (j = 0; j < g->ndomains; j++)
    {
      g->domain[j]->first_site = (j == 0 ? 0 : g->domain[j-1]->last_site + 1);
      g->domain[j]->last_site  = do_layer(g->domain[j], model) - 1; 
    };
  
  /* record number of sites in the model */
  g->nsites = g->domain[j-1]->last_site + 1;

  assert (g->nsites <= TOTALSITES);

  return g->nsites;
};


int do_layer (struct domain *dom, struct pdb_ATOM *model)
{
  struct layer curr_layer;

  int    i, new_site=-1;
  double z;

  /* loop over layers in domain */
  for (i = 0; i < dom->q; i++)
    {
      z = dom->z1 + i*dom->dz;
      curr_layer = setup_layer(z, dom);

      print_layer(&curr_layer);
      
      new_site = do_ring(&curr_layer, model);

    };

  return new_site;
};

struct layer setup_layer (double z, struct domain *d)
{
  struct layer l;

  double rA, z1, z2, rho1, rho2;

  rA = d->species->radius;   /* 1/2 interatomic distance (1/2 'bond' length) */

  l.rhomax = d->rho_outer;

  /* limits for calculating the shape of the wall: continuity
     currently quite messy; it should be possible to eliminate rho1/2
     in favor for one rho 
  */
  if (d->prev) {
    z1   = d->prev->z2;
    rho1 = d->prev->rho2;
  } else {
    z1   = d->z1;
    rho1 = d->rho1;
  };

  if (d->next) {
    z2   = d->next->z1;
    rho2 = d->next->rho1;
  } else {
    z2   = d->z2;
    rho2 = d->rho2;
  };

  /* determine inner radius for this layer at z */
  if (d->q > 1) {
    /* ie z1 != z2 */
    l.rhomin = d->rho(z, z1, z2, rho1, rho2);
  } else {
    /* degenerate case */
    l.rhomin = max(rho1, rho2);
  };
  assert(l.rhomin <= l.rhomax);  

  l.z      = z;
  l.maxring= round (l.rhomax - l.rhomin) / (2 * rA);
  /* rounding is not really correct here, but gives better results */

  if ( l.maxring == 0) {
    /* only enough space for the outermost ring */
    l.rhomin = l.rhomax;
  };

  l.dr     = (l.rhomax - l.rhomin) / (double) l.maxring;
  l.domain = d;

  return l;
}



struct ring setup_ring (double rho, struct layer *l)
{
  struct ring this;

  this.maxsite = floor(PI/(asin( l->domain->species->radius / rho )));
  this.dphi    = 2*PI / this.maxsite;
  this.rho     = rho;
  this.layer   = l;

  return this;
};


int do_ring (struct layer *l, struct pdb_ATOM *model)
{
  struct ring curr_ring;
  double rho;
  int new_site=-1;

  for (rho = l->rhomin; rho <= l->rhomax; rho+=(l->dr))
    {
      curr_ring = setup_ring (rho, l);
      /* hack for flagging cavity exposed atoms */
      curr_ring.bExposed = (rho == l->rhomin);
      print_ring(&curr_ring);

      new_site = do_site (&curr_ring, model);
    };

  return new_site;
}; 


int do_site (struct ring *r, struct pdb_ATOM *model)
{
  static int new_site;   /* initialized to 0 by compiler */
  struct pdb_ATOM *s;
  int    i;

  /* actually it is completely brain damaged to have string pointers
     in the pdb struct; one needs the strings to be somewhere between
     fuction calls... anyway, I am using static chars[] for the time being */
  static char RES[]     = "MTH";
  static char EXPOSED[] = "EXPD";
  static char EMPTY[1]  = "";
  char   *segment;

  Cylindrical u;

  /* flag inner wall atoms */
  segment = (r->bExposed) ? EXPOSED : EMPTY;

  for (i = 0, u.phi = 0; i < r->maxsite; i++,  u.phi += r->dphi)
    {
      s = &model[i + new_site]; 
      s->serial     = i + new_site + 1;  /* serial is NOT the offset in the array*/
      s->species    = r->layer->domain->species;
      s->name       = s->species->symbol;
      s->altLoc     = " ";
      s->resName    = RES;
      s->chainID    = " ";
      s->resSeq     = 1;    /* all pore atomes belong to the same res/chain */
      s->iCode      = " ";

      u.rho = r->rho;
      u.z   = r->layer->z;
      s->cpos    = u;
      s->pos     = cyl2cart(u);

      s->occupancy  = 0; /* dont care, dont know... */
      s->tempFactor = 0;
      s->segID      = segment;
      s->element    = " ";
      s->charge     = s->species->charge;
  
      /* GROMACS atoms section */
      /* buffer overflow possible (max length 10): */
      /* sprintf (s->atom, "%s%d", s->species->gmx_name, s->serial); 
         assert (strlen(s->atom) <= 10);
      */
      strncpy (s->atom, s->species->gmx_name, 9);
      s->cgnr = s->serial;     
                       /* all atoms in their own charge group.  Perhaps
			  better to divide it up into domains, but for
			  the time being... 
		       */
      print_site(s);
    };

  return new_site += i;  /* first site of next ring */
}

void unitcell (struct geom *g) {
  struct pdb_CRYST1 *u;

  u = &(g->unitcell);

  /* unitcell size */
  u->a = u->b = 2 * g->domain[0]->rho_outer;
  u->c = fabs (g->domain[g->ndomains - 1]->z2 - g->domain[0]->z1);

  /* tetragonal unitcell */
  u->alpha = u->beta = u->gamma = 90.00;  
  
  u->eGroup = "P 1";
  u->z = 1;
  
  return;
};

void cavitybox (struct geom *g) {
  /* tetragonal bounding box for the cavity of the pore */
  struct pdb_CRYST1 *c;
  struct domain *top, *bottom;

  c = &(g->cavitybox);
  top    = g->domain[g->ndomains - 1];
  bottom = g->domain[0];

  c->a = c->b = 2 * max( max(top->r1,    top->r2),
			 max(bottom->r1, bottom->r2));
  c->c = fabs (top->z2 - bottom->z1) 
       + top->species->radius + bottom->species->radius;

  c->alpha = c->beta = c->gamma = 90.00;  
  
  c->eGroup = "P 1";
  c->z = 1;
  
  return;
};


/* volume of a domain */
real volume (struct potpars *pp, struct pdb_ATOM model[], struct domain *dom, 
	      struct pprofile *prf) {
  real vol, Rmax, zmin, zmax;
  real min;         /* global (hopefully) minimum in the total potential */
  int i,first,last;

  /* take all atoms of this and neighboring domains: "Take them
     out. All of them!" (Senator Palpartine aka Darth Sidious) */
  first = (dom->prev) ? dom->prev->first_site : dom->first_site;
  last  = (dom->next) ? dom->next->last_site  : dom->last_site;

  pp->ncenters=last-first + 1;

  mesg(SUB1,"volume(): Using %d atoms (%d -> %d), centered on domain '%s'.",
       pp->ncenters,first,last,dom->description);

  Rmax = (prf->bSet) ? prf->Rmax :  dom->rho1/10.0;
  zmin=(dom->z1 - dom->dz/2.0)/10.0;
  zmax=(dom->z2 + dom->dz/2.0)/10.0;

  /* setup the function to be integrated */
  /* (1) find the minimum
     - requires a function in cartesian coordinates
     -> setup the centers in cartesian
     - use cartesian vlj()
  */
  { 
    Cartesian pos;
    pp->xyzcenters=grid2_alloc(pp->ncenters,3);
    for(i=0;i<pp->ncenters;i++) {
      pos = cyl2cart(model[i+first].cpos);
      pp->xyzcenters[i][XX]=pos.x/10.0;
      pp->xyzcenters[i][YY]=pos.y/10.0;
      pp->xyzcenters[i][ZZ]=pos.z/10.0;
    }
  
    init_vljcyl(pp);
    mesg(VERBOSE,"Potential: %s with ffgmx OW-CH4 interaction parameters at T=%g K", 
	 "Lennard-Jones 12-6 V_LJ(x,y,z)", pp->Temp);

    pp->min=find_min(vljvec,Rmax,zmin,zmax,NULL);
    mesg(VERBOSE,"Minimum: %f\n",pp->min);
  }
  
  /* (2) calculate the volume
     - this is better done in cylindrical coordinates
     --> setup again (and use cylindrical LJ)
  */
  /* These cylindrical coordinates are buried in structures; I rather
     have them as simple arrays: AND Length has to be in nm (not
     Angstrom) because this is the length unit in the Lennard-Jones
     parameters (currently CH4-OW hardcoded) */
  pp->centers=grid2_alloc(pp->ncenters,3);
  for(i=0;i<pp->ncenters;i++) {
    pp->centers[i][RAD]=model[i+first].cpos.rho/10.0;
    pp->centers[i][PHI]=model[i+first].cpos.phi;
    pp->centers[i][ZZZ]=model[i+first].cpos.z/10.0;
  }
  init_vljcyl(pp);   /* now with the minimum found and with the
                        LJ-centers in cylindrical coordinates! */
  mesg(VERBOSE,"Potential: Shifted Lennard-Jones 12-6 V_LJ(r,phi,z) - %g kT "
       "with ffgmx OW-CH4 interaction parameters at T=%g K", pp->min, pp->Temp);

  init_gaussleg(pp->ngaussleg);
  mesg(VERBOSE,"Using %d-point Gauss-Legendre quadrature for the volume integrals.",
       pp->ngaussleg);

  if (prf->bPlot)  plot_potential(vljcyl,Rmax,zmin,zmax,prf->nzplot);


  /*
  mesg(VERBOSE,"Potential: %s with ffgmx OW-CH4 interaction parameters", 
       (prf->pot == vRljcyl) ? "WCA repulsive V_R,LJ(r)" : "Lennard-Jones 12-6 V_LJ(r)");
  */

#ifdef TESTCASE
  { 
    real vexact;
    /* test case: this should give the exact volume */
    mesg(SUB1,"-------> Internal test of volume integration <-------");

    vexact=PI*Rmax*Rmax*(zmax-zmin);
    mesg(SUB1,"VOLUME_TEST: R[nm] <%f> R_c[nm] <%f> L[nm] <%f> V=¶·R_c²·L[nm³] <%f>", 
	 dom->r1/10.0, Rmax, zmax-zmin, vexact); 

    vol=cylint(cylconst, 0,Rmax, 0, 2*PI, zmin,zmax); 
    mesg(SUB1,"VOLUME_TEST: f=1: Volume <%f> V_exact <%f>",  vol,vexact); 

    vol=cylint(cylcos, 0,Rmax, 0, 2*PI, zmin,zmax); 
    vexact=0;
    mesg(SUB1,"VOLUME_TEST: f=cos(phi): Volume <%f> V_exact <%f>",  vol,vexact); 

    vol=cylint(cyl3, 0,Rmax, 0, 2*PI, 0,zmax-zmin); 
    vexact=(1-exp(-Rmax))*PI*pow(zmax-zmin,3)/3.0;
    mesg(SUB1,"VOLUME_TEST: f=1/r exp(-r)*cos(phi)*cos(phi)*z: Volume <%f> V_exact <%f>",  
	 vol,vexact); 

    vexact=1.0; /* PI*Rmax*Rmax*(zmax-zmin); */
    vol = tdavg(Qvol,Zero, Rmax,zmin,zmax);
    mesg(SUB1,"VOLUME_TEST: <Qvol>, V(r)=0: Volume <%f> V_exact <%f>",  
	 vol,vexact); 

    /* this analytical soln is alraedy specialised for FHG=6 */
    vexact=(1.0-4.0*exp(-3.0))/(1-exp(-6.0*Rmax)*(1.0+6.0*Rmax));
    vol = tdavg(Qvol,LinCheck, Rmax,zmin,zmax);
    mesg(SUB1,"VOLUME_TEST: <Qvol>, V(r)=6r: Volume <%f> V_exact <%f>",  
	 vol,vexact); 

    /* correct one, f = 6, gives the same as above */
#define FHG 6.0
    vexact= -(exp(FHG*(-0.5 + Rmax))*
	      (-2.0 + 2.0*exp(FHG/2.) - FHG))/ (2.*(1.0 - exp(FHG*Rmax) + FHG*Rmax));
    vol = tdavg(Qvol,LinCheck, Rmax,zmin,zmax);
    mesg(SUB1,"VOLUME_TEST: <Qvol>, V(r)=f/beta r nm^-1: Volume <%f> V_exact <%f>",  
	 vol,vexact); 

    /* g(E) * exp(-bE) * f() 
       f(r) = K1/2 r^2        (integral from Mathematica)
     */
    vexact=2.0*PI*(zmax-zmin)*
           (K1*pow(Rmax,2) - (3.0*sqrt(2.0/PI)*sqrt(K1*pow(Rmax,2)))/
	    exp((K1*pow(Rmax,2))/2.) + 
	    (3.0 - K1*pow(Rmax,2))*
	    erf(sqrt(K1*pow(Rmax,2))/sqrt(2)))/(2.*K1);
    vol = xV(harmonic,Rmax,zmin,zmax);
    mesg(SUB1,"VOLUME_TEST: xV: V(r)= k/2b r² nm^-1: Volume <%f> V_exact <%f>",  
	 vol,vexact); 
      
    mesg(SUB1,"------> End of testcases <-------\n");
  }
#endif /* TESTCASE */ 

  /* total volume */

  /* simple & inefficient (two integrations):
     V = <Q> = Tr Qexp(-beta H) / Tr exp(-beta H)
  */
  mesg(INPUT,"VOLUME_INPUT: R[nm] <%f> R_c[nm] <%f> L[nm] <%f> V=¶·R_c²·L[nm³] <%f>",
       dom->r1/10.0, Rmax, zmax-zmin, PI*Rmax*Rmax*(zmax-zmin));
  mesg(INPUT,"VOLUME_PARAMETERS: Rmax[nm] <%f>  z1 <%f>  z2 <%f>",Rmax,zmin, zmax);

  /* accessible volume, Labbook II, p84 */
  vol = xV(vljshiftcyl,Rmax,zmin,zmax);
  mesg(INPUT,"VOLUME_DATA: xV Volume[nm³] <%f> R*[nm] <%f>\n", 
       vol,sqrt(vol/(PI*(zmax-zmin))));
  
  /* normalised configurational volume, Labbok II, p79 (rubbish!) */
  vol = Zsum(vljshiftcyl,Rmax,zmin,zmax);
  mesg(INPUT,"VOLUME_DATA: v1/Z Volume[nm³] <%f> R*[nm] <%f>\n", 
       vol,sqrt(vol/(PI*(zmax-zmin))));


  /* Thermodynamic average at fixed particle energy, Labbook II, p78 */
/*    vol = tdavg(Qvol,prf->pp->u1, Rmax,zmin,zmax); */
/*    mesg(INPUT,"VOLUME_DATA: <Qvol> Volume[nm³] <%f> R*[nm] <%f>\n",  */
/*         vol,sqrt(vol/(PI*(zmax-zmin)))); */

  /* Pure configurational volume (unnormalised) */
/*    vol = cylint(Z_V, 0,Rmax, 0,2*PI, zmin,zmax); */
/*    mesg(INPUT,"VOLUME_DATA: v1 Volume[nm³] <%f> R*[nm] <%f>\n",  */
/*         vol,sqrt(vol/(PI*(zmax-zmin)))); */

  if (prf->bSet)  calc_profile(prf);
  
  free_gaussleg();
  free(pp->centers);
  free(pp->xyzcenters);
  return vol;
}


/* pore profile 
   calculate the effective radius in thin slices, spaced by dz
*/
void calc_profile (struct pprofile *prf)  {
  int i;
  real deltaZ;
  real Rmax=prf->Rmax;
  real zmin=prf->z1;
  real zmax=prf->z2;

  mesg(VERBOSE,"Calculating pore profile (%d z-slices)",prf->nz);
  mesg(INPUT,"pore profile: Rmax <%f>  zmin <%f>  zmax <%f>  [nm]",
       Rmax,zmin,zmax);
  
  prf->r=grid2_alloc(prf->nz,2);

  deltaZ = (zmax - zmin)/(real)prf->nz;
  for(i=0;i<prf->nz;i++) {
    fprintf(stderr,"\rIntegrating slice %6d (z=%2.3fnm)  [%5.1f%%]   ",
	    i,zmin+(i+0.5)*deltaZ,((real)i+1.0)*100/(real)prf->nz);
    prf->r[i][0]=zmin + (i+0.5)*deltaZ;
    prf->r[i][1]=sqrt(
	  xV(vljshiftcyl, Rmax, zmin+i*deltaZ, zmin+(i+1)*deltaZ) / (PI*deltaZ)  );
  }
  fprintf(stderr,"\n");
  return;
}

/* append int+1 levels (upto value max) to the file */
void xf_append_levels (FILE *fData, int ncols, real min, real max) {
  int i;
  real delta;

  delta = (max-min)/(ncols);
  fprintf(fData,"%d\n",ncols);        /* number of levels */
  for(i=1;i<=ncols;i++)                /* equally spaced */
    fprintf(fData, "%d  %f\n",i,min+i*delta);
  return;
}

void xf_append_axes_annotation (FILE *fData, real x1, real x2, 
			     real y1, real y2) {
  fprintf(fData,"%f  %f\n",x1, x2);
  fprintf(fData,"%f  %f\n",y1, y2);
  return;
}

/* plot the potential (xfarbe format) */
void plot_potential (real (*f)(real,real,real),real Rmax,real zmin, real zmax,
		     int nslices)  { 
    int  i,j,k;
    char buf[STRLEN];
    char *fn = buf;
    FILE *XF;
    real xmin,xmax,ymin,ymax,dx,dy,dz;
    real x,y,z;
    int  nx,ny,nz;
#define XF_RES 0.01   /* nm */
    xmin=ymin=-Rmax;
    xmax=ymax=Rmax;
    zmin -= 0.5;
    zmax += 0.5;
    
    nx = (int)((xmax-xmin)/XF_RES)+1;
    dx=XF_RES;
    ny = (int)((ymax-ymin)/XF_RES)+1;
    dy=XF_RES; 
    nz = nslices;
    dz=(zmax-zmin)/nslices;
    
    for(k=0; k<nslices;k++) {
      z=zmin+k*dz;
      sprintf(fn,"xy_%03d.dat",k);
      XF=fopen(fn,"w");
      mesg(VERBOSE,"Writing potential at z=%f into %s",z,fn);
      fprintf(XF,"LJ-potential at z=%f\n",z);
      fprintf(XF,"%d %d\n",nx,ny);
      for(j=0;j<ny;j++){
	y=ymin+j*dy;
	for(i=0;i<nx;i++) {
	  x=xmin+i*dx;
	  fprintf(XF,"%f ",f(sqrt(x*x+y*y),atan2(y,x),z));
	}
	fprintf(XF,"\n");
      }
      xf_append_levels(XF,32,-2,4);
      xf_append_axes_annotation(XF,xmin,xmax,ymin,ymax);
      fclose(XF);
    }

    /* same thing again for xz */
    nx = (int)((xmax-xmin)/XF_RES)+1;
    dx=XF_RES;
    ny = nslices;
    dy=(ymax-ymin)/nslices; 
    nz = (int)((zmax-zmin)/XF_RES)+1;
    dz=XF_RES;

    for(j=0;j<nslices;j++){
      y=ymin+j*dy;
      sprintf(fn,"xz_%03d.dat",j);
      XF=fopen(fn,"w");
      mesg(VERBOSE,"Writing potential at y=%f into %s",y,fn);
      fprintf(XF,"LJ-potential at y=%f\n",y);
      fprintf(XF,"%d %d\n",nx,nz);
      for(k=0; k<nz;k++) {
	z=zmin+k*dz;
	for(i=0;i<nx;i++) {
	  x=xmin+i*dx;
	  fprintf(XF,"%f ",f(sqrt(x*x+y*y),atan2(y,x),z));
	}
	fprintf(XF,"\n");
      }
      xf_append_levels(XF,32,-2,4);
      xf_append_axes_annotation(XF,xmin,xmax,ymin,ymax);
      fclose(XF);
    }


}




int pdb_write_header (FILE *fp, struct geom *g)
{
  int i, error;

  struct pdb_HEADER h = {
    "Atomistic Model for a Transmembrane Channel",
    " ",
    " "
  };

  struct pdb_REMARK r = {
    6,
    ""
  };

  struct domain *d;

  error = pdb_w_header (fp, &h);

  for (i = 0; i < g->ndomains && i <= MAXDOMAINS; i++)
    {
      /* basically print_domain () */
      d = g->domain[i];
      sprintf (r.remark, "  Type %d: %s", d->type, d->description);
      error = pdb_w_remark (fp, &r);
      mesg (SUB1, "writing... [%s]", r.remark);

      sprintf (r.remark, "  Defaultspecies %s <%s> rA <%4.2f> charge <%-2s>", 
			  d->species->name, d->species->symbol, d->species->radius, 
			  d->species->charge);
      error = pdb_w_remark (fp, &r);

      sprintf (r.remark, "  Length <%6.3f>  r1   <%6.3f>   r2 <%6.3f>   r_outer <%6.3f>",
			  d->l, d->r1, d->r2, d->r_outer);
      error = pdb_w_remark (fp, &r);

      sprintf (r.remark, "                   rho1 <%6.3f> rho2 <%6.3f> rho_outer <%6.3f>",
			  d->rho1, d->rho2, d->rho_outer);
      error = pdb_w_remark (fp, &r);

      
      sprintf (r.remark, "  z1 <%6.3f>  z2 <%6.3f>  Delta_z <%6.3f>   q <%d>",
			  d->z1, d->z2, d->dz, d->q);
      error = pdb_w_remark (fp, &r);

      sprintf (r.remark, "  first_site <%d>  last_site <%d>",
			  d->first_site, d->last_site);
      error = pdb_w_remark (fp, &r);
    };
  
  return error;
};


void print_geom(struct geom *g) 
{
  int i;

  if (debuglevel <= OFF) return;

  printf("Geometry data\nR_outer <%5.2f>  N_domains <%d>  N_sites <%d> N_bonds <%d> \n",
	 g->r_outer, g->ndomains, g->nsites, g->nbonds);
  printf("--------------------------------------------------------------------\n");

  if (debuglevel < INPUT) return; 

  for (i = 0; i < g->ndomains && i <= MAXDOMAINS; i++)
    {
      print_domain(g->domain[i]);
    };
  printf("\n");
  
  return;
};

void print_domain (struct domain *d)
{
  if (debuglevel < INPUT) return; 

  printf("  Type %d: %s\n", d->type, d->description);
  printf("  Defaultspecies %s <%s> rA <%4.2f> charge <%-2s>\n", 
	 d->species->name, d->species->symbol, d->species->radius, 
	 d->species->charge);
  printf("  Length <%6.3f>  r1   <%6.3f>   r2 <%6.3f>   r_outer <%6.3f>\n",
	 d->l, d->r1, d->r2, d->r_outer);
  printf("                   rho1 <%6.3f> rho2 <%6.3f> rho_outer <%6.3f>\n",
	       d->rho1, d->rho2, d->rho_outer);
  
  printf("  z1 <%6.3f>  z2 <%6.3f>  Delta_z <%6.3f>   q <%d>\n",
	 d->z1, d->z2, d->dz, d->q);
  printf("  first_site <%d>  last_site <%d>\n",
	 d->first_site, d->last_site);
  printf("\n");

  return;
};



void print_layer (struct layer *l)
{
  extern int debuglevel;

  if (debuglevel < SUB1) return;

  printf("layer: rhomin <%5.3f> rhomax <%5.3f>  z <%6.3f>  maxring <%d> dr <%5.3f>\n",
	 l->rhomin, l->rhomax, l->z, l->maxring, l->dr);
  return;
};


void print_ring (struct ring *r)
{
  if (debuglevel < SUB2) return;

  printf("ring: rho <%6.3f>  maxsite <%d>  delta_phi <%6.3f> %s\n",
	 r->rho, r->maxsite, r->dphi, (r->bExposed) ? "EXPOSED" : "");
  return;
};


void print_site (struct pdb_ATOM *s)
{
  extern int debuglevel;

  if (debuglevel < COORDINATES) return;
  printf("[%4d] rho <%8.3f>  phi <%8.3f>  z <%8.3f> Spec: <%s>\n",
	 s->serial, s->cpos.rho, s->cpos.phi, s->cpos.z, s->name);
  printf("       x   <%8.3f>  y   <%8.3f>  z <%8.3f> Res:  <%s>\n",
	 s->pos.x, s->pos.y, s->pos.z, s->resName);

  return;
};




void print_usage (char *progname,char *usage)
{
  printf("Usage: %s [OPTIONS]\n%s", progname,usage);
  return;
};


/* change on commandline -debug <val> */
int debuglevel = INPUT;


struct spec species[] = {
  { METHANE, "Methane", "CH4", 1.95, 1.95, "0", 
    "CH4", 16.0430, 0.000, "A", 0, 0 },             /* CH4 c6, c12 ????? */
  { CARBON, "Carbon",  "C",   1.4, 1.4, "0" },
  { OXYGEN, "Oxygen",   "O", 1.5, 1.5, "0" },
  { METHYL, "Methyl", "CH3", 1.4, 1.95, "0",
    "CH3", 15.03500, 0.000, "A", 0.88765e-02,  0.26150e-04}
};


	 
int main (int argc, char *argv[]) {
  static char usage[] = {
      "Calculate positions of pseudo-atoms in a model for a transmembrane channel and\n" 
      "write output to FILE or pore.pdb and pore.itp.\n\n"
      "  -h\t\t show help\n"
      "  -v\t\t be verbose (= -debuglevel 30 )\n"
      "  -debug <NUM>\t set debuglevel (0..100)\n"
      "  -spec <NUM>\t default species id\n"
      "  -showspec\t show hard coded species\n\n"
      "  -o FILE\t pdb coordinate file\n"
      "  -s FILE\t itp topology file\n"
      "\n"
      "  -f <file>\t read pore description from file (see below)\n"
      "  -R <r>\t Outer radius of the model\n"
      "  -P <r> <l>\t Pore region: inner radius and length\n"
      "  -M <r> <l>\t Mouth region: largest inner radius and length\n"
      "  -b dmin dmax g_min g_max   (repeatable)\n"
      "\t\tform bonds with angle gamma when atoms are no further apart\n"
      "\t\tthan dmax Ang and g_min <= gamma <= g_max, g from [0°..90°] )\n"
      "  -c\t\t only write connectivity to output, no bond length or kB, kA\n" 
      "  -x\t\t neither connectivity nor bond length to output (isolated atoms)\n" 
      "  -kB <c>\t Force constant of bonds, in kJ mol^-1 nm^-2\n"
      "  -kA <c>\t Force constant of angles, in kJ mol^-1 rad^-2\n"
      "  -cc\t\t Center ccordinates on cavitybox, not on unitcell\n"
      "\n"
      "Pore volume calculation (setting any of these switches on profile calculation):\n"
      "In order to enable volume calculation, set -volume explicitly!\n"
      "ATTENTION: all these LENGTHs are in NANO METRE not Angstrom !\n"
      "  -profile [<file>]  calculate the profile in addition to the volume\n"
      "  -z1, -z2 <z>       profile between z1 and z2\n"
      "  -Rmax <r>          integrate out to Rmax (also use for the total volume\n"
      "                     integration if -profile is set)\n"
      "  -npoints <N>       number of points per dimension in the integrals\n"
      "  -T temp            Temperature in Kelvin [300]\n"
      "  -wca               If set, only use the repulsive part of the Lennard-Jones\n"
      "                     potential (split after Weeks, Chandler & Andersen [1971])\n"
      "  -plot              xfarbe output of the potential in z slices\n"
      "  -nzplot            number of plot slices \n"
      "\nDescription of the input file:\n"
      "------------------------------\n"
      "Instead of using -R (RADIUS), -M (MOUTH), and -P (PORE) one can describe the system\n"
      "in a more flexible manner with a geometry in put file. It can contain up to "
      "MAXDOMAINS domains (i.e. MOUTH and PORE lines). Allowed lines:\n"
      "# comment (skipped)\n"
      "# RADIUS is the global outer radius (in Angstrom)  of the cylinder\n"
      "RADIUS r_outer\n"
      "# domain type and radius at the upper and lower end of the domain;\n"
      "# r_lower of domain i and r_upper of domain i+1 are typically identical\n"
      "MOUTH r_upper r_lower length [species]\n"
      "PORE  r_upper r_lower length [species]\n"

      
  };
  int i, n_atoms, n_bonds, n_angles, n_bc;
  int error;
  real vol;
  FILE *fp;

  /* 
     segmentation fault if arrays to large 
     --->> 
           now that I have learned to  calloc I should rewrite this 
           on purely aesthetical grounds !
     <<---- 
  */
  struct pdb_ATOM model[TOTALSITES];
  struct itp_bond bonds[MAXBONDS]; 
  struct itp_angle angles[MAXANGLES];
  struct geom geometry;
  struct std_input in;
  struct potpars  potential = {
    NULL,
    TEMPERATURE,
    NFREEDOM_SPC,
    CSIX_OW_MTH, CTWELVE_OW_MTH, 0.0,
    0, NULL, NULL,
    N_GAUSSLEG
  };
  struct pprofile profile = {
    "LJprofile.dat",
    FALSE,
    FALSE,
    POT_PLOT_SLICES,
    NULL,
    NZPROF,
    1.5,
    -2,2,
    NULL
  };

  /* argument processing */
  /* *** no sanity checks *** */

  n_atoms  = 0;  /* number of sites ('atoms') in the model */
  n_bonds  = 0;  /* number of bonds */
  n_angles = 0;  /* number of angles between bonds */
  n_bc     = 0;  /* number of constraint conditions for generating bonds */

  if (argc < 2) {
    printf("Running with default values.\n\nType %s -h for help.\n", argv[0]);
    debuglevel = IMPORTANT;
  };

  /* initialize defaults */
  in = default_input ();
  potential.u1 = vljcyl;     /* use the full Lennard-Jones potential in
                               the configurational volume calculations
                               by default */
  profile.pp = &potential;

  /* rudimentary opt-processing.. yarch */
  for (i = 1; i < argc; i++)
    {
      if (!strcmp(argv[i], "-debug")) { 
	debuglevel = atoi(argv[++i]);
      } else if (!strcmp(argv[i], "-v")) {
	debuglevel = VERBOSE; 
      } else if (!strcmp(argv[i], "-h")) {
	print_usage(argv[0],usage);
	exit(1);
      } else if (!strcmp(argv[i], "-showspec")) {
	print_species ();
	exit(1);
      } else if (!strcmp(argv[i], "-o")) {
	in.coordfile = argv[++i]; 
      } else if (!strcmp(argv[i], "-s")) {
	in.topofile = argv[++i];
      } else if (!strcmp(argv[i], "-f")) {
	in.datafile = argv[++i]; 
      } else if (!strcmp(argv[i], "-R")) {
	in.r_outer = atof(argv[++i]);
      } else if (!strcmp(argv[i], "-M")) {
	in.r_mouth = atof(argv[++i]);
	in.l_mouth = atof(argv[++i]);
      } else if (!strcmp(argv[i], "-P")) {
	in.r_pore = atof(argv[++i]);
	in.l_pore = atof(argv[++i]);
      } else if (!strcmp(argv[i], "-spec")) {
	in.specid = atoi(argv[++i]);
      } else if (!strcmp(argv[i], "-b")) {
	if (n_bc > MAXBONDCONSTRAINTS) {
	  fatal_error (1,
		       "Error: too many bond constraints, maximum is %d.\n",  
		       MAXBONDCONSTRAINTS);
	}
	in.bc[n_bc].serial    = n_bc;
	in.bc[n_bc].dmin      = atof(argv[++i]);
	in.bc[n_bc].dmax      = atof(argv[++i]);
	in.bc[n_bc].gamma_min = atof(argv[++i]);
	in.bc[n_bc].gamma_max = atof(argv[++i]);
	in.n_bc               = ++n_bc;
      } else if (!strcmp(argv[i], "-c")) {
	in.connectonly  = TRUE;
      } else if (!strcmp(argv[i], "-x")) {
	in.atomsonly    = TRUE;
      } else if (!strcmp(argv[i], "-cc")) {
	in.shiftcbox    = TRUE;
      } else if (!strcmp(argv[i], "-kB")) {
	in.k_bond  = atof(argv[++i]);
      } else if (!strcmp(argv[i], "-kA")) {
	in.k_angle = atof(argv[++i]);
      } else if (!strcmp(argv[i], "-volume")) {
	profile.bSet = TRUE;
      } else if (!strcmp(argv[i], "-z1")) {
	profile.bSet = TRUE;
	profile.z1 = atof(argv[++i]);
      } else if (!strcmp(argv[i], "-z2")) {
	profile.bSet = TRUE;
	profile.z2 = atof(argv[++i]);
      } else if (!strcmp(argv[i], "-Rmax")) {
	profile.bSet = TRUE;
	profile.Rmax = atof(argv[++i]);
      } else if (!strcmp(argv[i], "-profile")) {
	profile.bSet = TRUE;
	if (i < argc-1 && argv[i+1][0] != '-') 
	  strncpy(profile.fn,argv[++i],STRLEN);
      } else if (!strcmp(argv[i], "-npoints")) {
	potential.ngaussleg = atof(argv[++i]);
      } else if (!strcmp(argv[i], "-T")) {
	potential.Temp = atof(argv[++i]);
      } else if (!strcmp(argv[i], "-wca")) {
	potential.u1 = vRljcyl;
      } else if (!strcmp(argv[i], "-plot")) {
	profile.bPlot = TRUE;
      } else if (!strcmp(argv[i], "-nzplot")) {
	profile.bPlot = TRUE;
	profile.nzplot = atoi(argv[++i]);
      } else {
	mesg(OFF,"Unknown option: %s",argv[i]);
      };
    };

  mesg (ALL, "TOTALSITES %d\n", TOTALSITES);
  mesg (ALL,
   "MAXBONDS %d\nsizeof struct bond %d, sizeof bonds[] %d\n", 
    MAXBONDS, sizeof (struct itp_bond), sizeof (bonds));
  mesg (ALL, 
   "MAXANGLES %d\nsizeof struct angle %d, sizeof angles[] %d\n", 
    MAXANGLES, sizeof (struct itp_angle), sizeof (angles));



  error = input_geom (&in, &geometry);
  setup_domain (&geometry);
  unitcell  (&geometry);
  cavitybox (&geometry);

  n_atoms  = do_coordinates (model, &geometry);

  if (!geometry.atomsonly) {
    n_bonds  = do_bonds (bonds, model, &geometry);
    n_angles = do_angles (angles, bonds, &geometry);
  } else {
    n_bonds = n_angles = 0;
  }

  print_geom (&geometry);    
  
  mesg (WARN, "\nNumber of sites:  %d", n_atoms);
  if (!geometry.atomsonly) {
    mesg (WARN, 
    "Number of bonds:  %d within max cut-off %5.2f Ang (no double-counting)", 
	n_bonds, find_dmax (geometry.bc, geometry.nbc));
    mesg (WARN, "Number of angles: %d (no double-counting)\n", n_angles);
  }

  center (model, &geometry);
  
  error = write_topology (model, bonds, angles, &geometry);
  error = write_pdb (model, bonds, &geometry);

  /* 
     calculate pore volume, based on J. S. Rowlinson, J Chem Soc,
     Faraday Trans. 2, 82 (1986), 1801, (which didnt really work), and
     discussion with Andrew Horsefield.  */

  if (profile.bSet) {
    vol = volume(&potential,model,geometry.domain[1],&profile);

    /* pore profile (already calculated in volume ) */
    fp=fopen(profile.fn,"w");
    if (fp) {
      for(i=0;i<NZPROF;i++) {
        fprintf(fp,"%f  %f\n",profile.r[i][0],profile.r[i][1]);
      }
      fclose(fp);
    }
    free(profile.r);
  }

  return error < 0 ? 1 : 0;
};





