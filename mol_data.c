/*
  $Id: mol_data.c,v 1.3 2000/12/07 13:04:21 oliver Exp $
  $Log: mol_data.c,v $
  Revision 1.3  2000/12/07 13:04:21  oliver
  - writes coorect CRYST1 unitcell to pdb file
  - shifts coordinates to center of unitcell before writing pdb
  - new functions vec_add() and translate_pdb_atom()

  Revision 1.1  2000/12/05 19:22:00  oliver
  Split pgeom seems to compile and to work (as buggy as usual...)


  - important data structures
  - output routines for writing gromacs topologies and pdb files
*/

#include "mol_data.h"

int pdb_w_atom (FILE *fp, struct pdb_ATOM *s)
{
  return fprintf(fp, 
     "%-6s%5d %4s%1s%3s %1s%4d%1s   %8.3f%8.3f%8.3f%6.2f%6.2f      %-4s%-2s%-2s\n",
	  "ATOM", s->serial, s->name, s->altLoc,
	  s->resName, s->chainID, s->resSeq, s->iCode,
	  s->pos.x, s->pos.y, s->pos.z, s->occupancy, s->tempFactor,
	  s->segID, s->element, s->charge);
};

int pdb_w_conect (FILE *fp, struct itp_bond *b) {
  /* simplified: one line per bond. For correct format of CONECT see: 
   http://www.rcsb.org/pdb/docs/format/pdbguide2.2/guide2.2_frame.html
  */
  return fprintf (fp, "%-6s%5d%5d\n",
		  "CONECT", b->a1->serial, b->a2->serial);
}; 

int pdb_w_cryst1 (FILE *fp, struct pdb_CRYST1 *c) {
  return fprintf (fp, "%-6s%9.3f%9.3f%9.3f%7.2f%7.2f%7.2f %-7s%4d\n",
		  "CRYST1", c->a, c->b, c->c, c->alpha, c->beta, c->gamma,
		  c->eGroup, c->z);
};

int pdb_w_header (FILE *fp, struct pdb_HEADER *h)
{
  return fprintf(fp, "%6s    %40s%9s   %4s\n", 
		  "HEADER", h->classification, h->date, h->IDcode);
}

int pdb_w_remark (FILE *fp, struct pdb_REMARK *r)
{
  return fprintf (fp, "%6s %3d %-59s\n", 
		  "REMARK", r->number, r->remark);
}

int pdb_w_end (FILE *fp)
{
  return fprintf (fp, "END\n");
}





int itp_w_header_moleculetype (FILE *fp) {
  /* "%5s %7d" */
  return fprintf (fp, "[ moleculetype ]\n;%4s %7s\n",
		  "name", "nrexcl");
};

int itp_w_header_atomtypes (FILE *fp) {
  /* "%5s  %12.4f %8.3f %6s %12g %12g" */
  return fprintf (fp, "[ atomtypes ]\n;%4s %12.4s %8s %6s %12s %12s\n",
		  "name", "mass", "charge", "ptype", "c6", "c12");
};


int itp_w_header_atoms (FILE *fp) {
  /* %6d %8s %8d %6s %6s %6d %6.3f %12.4f */
  return fprintf (fp, "[ atoms ]\n;%5s %8s %8s %6s %-6s %6s %6s %12s\n",
		  "nr", "type", "resnr", "res", "atom", "cgnr", "charge",
		  "mass");
};


int itp_w_header_bonds (FILE *fp) {
  /* "%5d %5d %5d %12g %12g" */
  return fprintf (fp, "[ bonds ]\n;%4s %5s %5s %12s %12s\n",
		  "ai", "aj", "funct", "d0 (nm)", "c1");
};

int itp_w_header_angles (FILE *fp) {
  /* "%5d %5d %5d %5d %12g %12g" */
  return fprintf (fp, "[ angles ]\n;%4s %5s %5s %5s %12s %12s\n",
		  "ai", "aj", "ak", "funct", "theta", "ctheta");
};

int itp_w_header_constraints (FILE *fp) {
  return fprintf (fp, "[ constraints ]\n;%4s %5s %5s %12s\n",
		  "ai", "aj", "funct", "dist (nm)");
};


int itp_w_bond (FILE *fp, struct itp_bond *b, Boolean connectonly) {
  return connectonly ? 
    fprintf (fp, "%5d %5d %5d\n",
		  b->a1->serial, b->a2->serial, 
		  b->type)
    :
    fprintf (fp, "%5d %5d %5d %12.4f %12e\n",
		  b->a1->serial, b->a2->serial, 
		  b->type, b->c0, b->c1);
};

int itp_w_atom (FILE *fp, struct pdb_ATOM *s) {
  return fprintf (fp, "%6d %8s %8d %6s %-6s %6d %6.3f\n", 
		  s->serial, s->species->gmx_name, s->resSeq,
		  s->resName, s->atom, s->cgnr, s->species->gmx_charge);
}; 

int itp_w_angle (FILE *fp, struct itp_angle *a, Boolean connectonly) {
  return connectonly ? 
    fprintf (fp, "%5d %5d %5d %5d\n",
		  a->a1->serial, a->a2->serial, a->a3->serial,
		  a->type)
    :
    fprintf (fp, "%5d %5d %5d %5d %12.3f %12e\n",
		  a->a1->serial, a->a2->serial, a->a3->serial,
		  a->type, a->theta, a->ctheta);
};

int itp_w_constraint (FILE *fp, struct itp_bond *b) {
  return fprintf (fp, "%5d %5d %5d %12.4f\n",
		  b->a1->serial, b->a2->serial, 
		  GMX_CONSTRAINT_1, b->c0);
};


int itp_w_moleculetype (FILE *fp) {
  /* currently hard coded until I have proper data structures for topology */
  return fprintf (fp, "%5s %7d\n",
		  GMX_MOLECULE_NAME, GMX_MOLECULE_NREXCL);
};  

