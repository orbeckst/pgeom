/*
   $Id: mol_data.h,v 1.5 2002/08/25 22:19:51 oliver Exp $
*/


#ifndef MOL_DATA_H_READ
#define MOL_DATA_H_READ

#include <stdio.h>
#include "xgtypes.h"
#include "util.h"

#define GMX_CONSTRAINT_1 1    /* constraints count as bonds */
#define GMX_CONSTRAINT_2 2    /* constraints do NOT count as bonds */
/* bonds */
#define GMX_BOND_TYPE 1       
/* angles */
#define GMX_ANGLE_TYPE 1      

/* [moleculetype] section of topology file */
#define GMX_MOLECULE_NAME  "CHL"
#define GMX_MOLECULE_NREXCL 3  /* exclude non-bonded interactions
                                  between atoms no further than nrexcl
                                  bonds apart */
struct spec {
  int    id;             /* access groups by a serial number */
  char   *name;          /* full chemical name */
  char   *symbol;        /* chemical formula */ 
  double radius;         /* typical bondlength */
  double vdwradius;      /* vdW radius -- not used */
  char   *charge;        /* eg +2, -1 */
  /* more stuff which is relevant to GROMACS */
  char   *gmx_name;
  double gmx_mass;
  double gmx_charge;    /* can be partial unlike pdb charge */
  char   *gmx_ptype;
  double gmx_c6;
  double gmx_c12;
};

struct pdb_ATOM {
  char *identifier;  /* "ATOM  " */
  int  serial;         /* %5d  start at 1 (because GROMACS likes it) */
  char *name;
  char *altLoc;
  char *resName;
  char *chainID;
  int  resSeq;         /* %4d */
  char *iCode;
  Cartesian pos;       /*   double x, y, z;     8.3f */
  double occupancy;    /* 6.2f */
  double tempFactor;   /* 6.2f */
  char *segID;
  char *element;
  char *charge; 
  /* additional entries for the purpose of pgeom */
  Cylindrical cpos;    /* position in cylindrical coordinates */
  struct spec *species;/* species at this site, element of species[] */
  /* GROMACS specific [ atoms ] (for topology)
     everything else is already in pdb fields */
  char atom[10];
  int  cgnr;
};

struct pdb_HEADER {
  char *classification;  /* 11-50 */
  char *date;            /* 51 - 59 */
  char *IDcode;          /* 63 - 66 */
};

struct pdb_REMARK {
  int  number;           /* 8-10 flushright */
  char remark[80];       /* 12-70 */
};

struct pdb_CRYST1 {
  double a;                /* 9.3 7-15 unitcell dim in Angstroem */
  double b;                /*    16-24 */
  double c;                /*    25-33 */
  double alpha;            /* 7.2 34-40 angles in degrees */
  double beta;             /*     41-47 */
  double gamma;            /*     48-54 */
  char   *eGroup;          /* 7   56-66 space group   */
  int    z;                /* 4   67 70 Z value */
};

struct itp_bond {
  struct pdb_ATOM *a1, *a2;  /* bond between atoms 1 and 2 */
  int    type;            /* GROMACS type of functional*/
  double c0;              /* GROMACS equilibrium bond length in nm */
  double c1;              /* GROMACS harmonic force constant kJ mol^-1 nm^-2 */
};

struct itp_angle {
  struct pdb_ATOM *a1, *a2, *a3; /* angle (atom1, atom2, atom3) */
  int    type;                   /* GROMACS type of functional*/
  double theta;                  /* equilib angle in deg */
  double ctheta;                 /* force const kJ mol^-1 rad^-2 */
};



extern int pdb_w_atom (FILE *, struct pdb_ATOM *);
extern int pdb_w_conect (FILE *, struct itp_bond *);
extern int pdb_w_cryst1 (FILE *, struct pdb_CRYST1 *);
extern int pdb_w_remark (FILE *, struct pdb_REMARK *);
extern int pdb_w_header (FILE *, struct pdb_HEADER *);
extern int pdb_w_end (FILE *);


extern int itp_w_header_moleculetype (FILE *);
extern int itp_w_header_atomtypes (FILE *);
extern int itp_w_header_atoms (FILE *);
extern int itp_w_header_bonds (FILE *);
extern int itp_w_header_constraints (FILE *);
extern int itp_w_header_angles (FILE *);
extern int itp_w_bond (FILE *, struct itp_bond *, Boolean);
extern int itp_w_constraint (FILE *, struct itp_bond *);
extern int itp_w_atom (FILE *, struct pdb_ATOM *);
extern int itp_w_angle (FILE *, struct itp_angle *, Boolean);
extern int itp_w_moleculetype (FILE *);

#endif /* MOL_DATA_H_READ */
