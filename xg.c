/*
  $Id: xg.c,v 1.5 2002/11/07 14:55:15 oliver Exp $
  $Log: xg.c,v $
  Revision 1.5  2002/11/07 14:55:15  oliver
  - some volume calculation stuff was not checked in
  - added METHYL (distance 0.14nm) to species
  - rewrote species id using enum

  Revision 1.4  2002/08/23 18:16:37  oliver
  - potential function for the volume can be set from the commandline
    (and is stored in struct profile.pot)
  - cleanup of header file mess

  Revision 1.3  2002/08/22 23:23:07  oliver
  Cleanup: remove last traces of these stupid *_proto.h files from the dawn of my carreer as a C programmer

  Revision 1.2  2002/08/22 02:49:25  oliver
  - Restructuring of headers
  - removed mgeom -- I'm not using it anyway
  - pgeom calculates a statistical mechanics estimate for the pore volume

  Revision 1.1  2000/12/13 16:40:52  oliver
  split functions in main (mgeom, pgeom), common (xg), and general (pg) -- stupid names...
  compile mgeom or pgeom in Makefile, use -DPROG
  all prototypes in separate files because my includes ar a crappy mess

  More or less the header file for pgeom.
*/

#include "util.h"
#include "mol_data.h"
#include "xg.h"

Cartesian atom_diff (struct pdb_ATOM *a, struct pdb_ATOM *b) {
  return vec_diff (a->pos, b->pos);
};

struct pdb_ATOM translate_pdb_atom (struct pdb_ATOM *a, Cartesian t) {
  struct pdb_ATOM this;

  this     = *a;
  this.pos = vec_add (a->pos, t);

  return this;
};



double find_dmax (struct bond_constraint *bc, int n_bc) {
  int i;
  double dmax[MAXBONDCONSTRAINTS];

  for (i = 0; i < n_bc; i++) {
    dmax[i] = bc[i].dmax;
  };
  
  piksrt (n_bc, dmax);

  return dmax[n_bc - 1];
};
  


int do_bonds (struct itp_bond *bonds, struct pdb_ATOM *model, struct geom *g) {
  int       n_bonds = 0;
  int       ai, aj, i;
  double    dd;      /* (bond length)^2 = r^2 */
  double     d;      /* bond length */
  Cartesian r;       /* bond vector = ri - rj */
  double dmax;
  struct bond_constraint *bc;
  int bc_class_count[MAXBONDCONSTRAINTS];

  for (i = 0; i < g->nbc; i++) {
    bc_class_count[i] = 0;
  };
  
  dmax = find_dmax (g->bc, g->nbc);   /* maximum cut-off */
  mesg (SUB1, "do_bonds(): initial bond cut-off dmax = %6.4f", dmax);

  for (ai = 0; ai < g->nsites - 1; ai++) {
    for (aj = ai + 1; aj < g->nsites;  aj++) {
      r  = atom_diff ( bonds[n_bonds].a1 = &(model[ai]), 
		       bonds[n_bonds].a2 = &(model[aj])  );
      dd = vec_dot  (r, r);
      if (dd > dmax*dmax) {
	continue;
      };

      /* is a nearest neighbour */
      d = sqrt(dd);
      
      mesg (COORDINATES, "bond: %d <-> %d, length <%6.4f>", 
	    ai, aj, d);
      mesg (INTERMEDIATES, 
	    "bond: (%6.3f, %6.3f, %6.3f) <-> (%6.3f, %6.3f, %6.3f)",
	    model[ai].pos.x, model[ai].pos.y, model[ai].pos.z,
	    model[aj].pos.x, model[aj].pos.y, model[aj].pos.z);
      mesg (INTERMEDIATES, "bond: r = (%6.3f, %6.3f, %6.3f)",
	    r.x, r.y, r.z);
      
      /* right direction ? */
      if ( (bc = bond_class (r, g->bc, g->nbc)) 
	         && d <= bc->dmax && d >= bc->dmin) {
	/* and right bond length for that direction ? */
	bonds[n_bonds].type = GMX_BOND_TYPE; 
	bonds[n_bonds].c0   = d/10;      /* equilib. bond length in nm !*/
	bonds[n_bonds].c1   = g->k_bond; /* force constant */

	n_bonds++;
	assert (n_bonds <= MAXBONDS);

	bc_class_count[bc->serial]++;
      };
    };
  };

  for (i = 0; i < g->nbc; i++) {
    mesg (IMPORTANT, 
	  "Number of bonds in bond class %3d (dmax = %5.3f Ang): %5d",
	  i, g->bc[i].dmax, bc_class_count[i]);
  };

  return g->nbonds = n_bonds;
};

struct bond_constraint *bond_class (Cartesian r, struct bond_constraint bc[],
				    int n_bc) {
  struct bond_constraint *this;
  int i;
  double gamma;
  const Cartesian z_axis = {0, 0, 1};
  
  this = NULL;
  gamma = rad2deg ( acos( vec_cos(r, z_axis) ) );
  gamma = gamma > 90 ? 180 - gamma : gamma;
    
  for (i = 0; i < n_bc; i++) {
    mesg (ALL, "bond_class(): %6.3f <= (g=%6.3f) <= %6.3f: ",
	  bc[i].gamma_min, gamma, bc[i].gamma_max); 
    if ( gamma <= bc[i].gamma_max && gamma >= bc[i].gamma_min ) {
      this = &(bc[i]);
      mesg (ALL, "%15sbond_class = %d (serial %d)", " ",
	    i, bc[i].serial); 
      break;
    };
  };

  return this;
};

int do_angles (struct itp_angle *angles, struct itp_bond *bonds,
	       struct geom *g) {
  int n_angles;
  int i, j;

  n_angles = 0;

  /* this works only for a sorted list of bonds ! */
  for (i = 0; i < g->nbonds; i++) {
    j = i;
    while ( bonds[++j].a1 == bonds[i].a1 ) {
      angles[n_angles] = new_angle ( bonds[i].a2, bonds[i].a1, bonds[j].a2,
				     g->k_angle);

      mesg (COORDINATES, "angle: %3d %3d %3d, theta %6.2f°",
	    angles[n_angles].a1->serial, angles[n_angles].a2->serial, 
	    angles[n_angles].a3->serial, angles[n_angles].theta);

      n_angles++;
      assert (n_angles <= MAXANGLES);
    };
  };

  return g->nangles = n_angles;
};


struct itp_angle new_angle (struct pdb_ATOM *a, struct pdb_ATOM *s, 
			    struct pdb_ATOM *b, double k_angle) {
  /* angle (origin s) between a and b */
  struct itp_angle this;
  Cartesian r1, r2;

  r1 = atom_diff (a, s);
  r2 = atom_diff (b, s);

  this.a1 = a;
  this.a2 = s;
  this.a3 = b;
  this.type = GMX_ANGLE_TYPE;
  this.theta  = rad2deg ( acos (vec_cos(r1, r2)) ); /* GROMACS in deg! */
  this.ctheta = k_angle;  /* force constant */

  return this;
};

void center (struct pdb_ATOM site[], struct geom *g) {
  Cartesian t;

  /* t =  c_box - c_model (c: center);
     model: c1 = c2 = 0; c3 = (z_top + z_bottom)/2 
     box:   cn = an/2
  */

  if (g->shiftcbox) {
    t.x = g->cavitybox.a / 2;
    t.y = g->cavitybox.b / 2;
    /* t.z = -z_bottom + rA_bottom */
    t.z = g->domain[0]->species->radius;
  } else {      /* center in unit cell */
    /* center in x-y plane */
    t.x = g->unitcell.a / 2;
    t.y = g->unitcell.b / 2;
    /* z is ok */
    t.z = 0;
  };
  center_model (t, site, g->nsites);
  return;
};

void center_model (Cartesian t, struct pdb_ATOM site[], int nsites) {
  /* translate coordinates to center of unitcell so that lower left
     corner is (0,0,0) */
  int i;

  mesg (VERBOSE, 
	"center_model(): Shifting all coordinates by (%5.3f, %5.3f, %5.3f).",
	t.x, t.y, t.z);

  for (i = 0; i < nsites; i++) {
    site[i] = translate_pdb_atom ( &(site[i]), t);
  };

  return;
};


int write_pdb (struct pdb_ATOM *model, struct itp_bond *bonds, struct geom *g)
{
  int i, error;
  FILE  *fp;
  
  error = 1;
  mesg (IMPORTANT, "write_pdb (): filename [%s]", g->coordfile);

  if ((fp = fopen(g->coordfile, "w")) == NULL) {
    fprintf(stderr, "write_pdb(): can't open %s",
	    g->coordfile);
    return 1;
  };

  error = pdb_write_header (fp, g);
  error = pdb_write_cavitybox (fp, &(g->cavitybox)); 
  error = pdb_w_cryst1 (fp, &(g->unitcell));
  for (i = 0; i < g->nsites && error > 0; i++) {
    mesg (ALL, "write_pdb (): write site [%d]", i);
    error = pdb_w_atom (fp, &(model[i]));
  };
  
  if (! g->atomsonly) {
    for (i = 0; i < g->nbonds && error > 0; i++) {
      mesg (ALL, "write_pdb (): write bond [%d]", i);
      error = pdb_w_conect (fp, &(bonds[i]));
    };
  };

  error = pdb_w_end (fp);
  fclose(fp);
  return error;
};

int pdb_write_cavitybox (FILE *fp, struct pdb_CRYST1 *cbox) {
  int error;
  struct pdb_REMARK r = {
    6,
    ""
  };
  
  sprintf (r.remark, "CAVITYBOX %9.3f%9.3f%9.3f", cbox->a, cbox->b, cbox->c);
  error = pdb_w_remark (fp, &r);

  return error;
};

int write_topology (struct pdb_ATOM *model, struct itp_bond *bonds, 
		    struct itp_angle *angles, struct geom *geometry){
  /* write GROMACS itp topology file */
  int error;
  FILE *fp;

  error = 1;

  mesg (IMPORTANT, "write_topology (): filename [%s]", geometry->topofile);

  if ((fp = fopen(geometry->topofile, "w")) == NULL) {
    fprintf(stderr, "write_topology(): can't open %s",
	    geometry->topofile);
    return 1;
  };

  error = write_moleculetype (fp);
  error = write_atoms  (fp, model, geometry);

  if (! geometry->atomsonly) {
    error = fprintf (fp, "#ifdef BONDS\n");
    error = write_bonds  (fp, bonds, geometry);
    error = fprintf (fp, "#else\n");
    error = write_constraints (fp, bonds, geometry);
    error = fprintf (fp, "#endif\n\n");
    error = write_angles (fp, angles, geometry);
  } else {
    mesg (IMPORTANT, 
	  "write_topology (): no bonds and angles information written");
  };

  return fclose(fp);
};


int write_moleculetype (FILE *fp){
  /* currently hard coded through GMX_MOLECULE_* */
  int error;

  error = itp_w_header_moleculetype (fp);
  error = itp_w_moleculetype (fp);

  return fprintf (fp, "\n");
};


int write_atoms (FILE *fp, struct pdb_ATOM *model, struct geom *g){
  int i, error;

  error = itp_w_header_atoms (fp);

  for (i = 0; i < g->nsites && error > 0; i++) {
    error = itp_w_atom (fp, &(model[i]));
  };

  return fprintf (fp, "\n");
};

int write_bonds (FILE *fp, struct itp_bond *bonds, struct geom *g){
  int  i, error;

  error = itp_w_header_bonds (fp);

  for (i = 0; i < g->nbonds && error > 0; i++) {
    mesg (ALL, "write_bonds (): write bond [%d]", i);
    error = itp_w_bond (fp, &(bonds[i]), g->connectonly);
  };

  return error;
};

int write_constraints (FILE *fp, struct itp_bond *bonds, struct geom *g){
  int  i, error;

  error = itp_w_header_constraints (fp);

  for (i = 0; i < g->nbonds && error > 0; i++) {
    mesg (ALL, "write_constraints (): write constraint [%d]", i);
    error = itp_w_constraint (fp, &(bonds[i]));
  };

  return error;
};


int write_angles (FILE *fp, struct itp_angle *angles, struct geom *g){
  int  i, error;

  error = itp_w_header_angles (fp);

  for (i = 0; i < g->nangles && error > 0; i++) {
    mesg (ALL, "write_angles (): write bond [%d]", i);
    error = itp_w_angle (fp, &(angles[i]), g->connectonly);
  };

  return error;
};

void print_species (void)
{
  int i;

  
  printf("Species (hardcoded):\n--------------------\n");
  printf("\tSpecID\tName\t\tSymbol\tr\tvdW radius [Ang]\n");
  for (i=0; i < NSPEC; i++){
    printf("\t<%d>\t<%s>\t<%s>\t<%.2f>\t<%.2f>\n",
	   species[i].id, species[i].name, species[i].symbol,
	   species[i].radius,species[i].vdwradius);
  };
  printf("\n");
  return;
};

