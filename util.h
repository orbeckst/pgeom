/* $Id: util.h,v 1.9 2008/02/17 16:55:54 oliver Exp $
   $Log: util.h,v $
   Revision 1.9  2008/02/17 16:55:54  oliver
   - fixed DOUBLE compilation
   - default BIN_DIR to ~/bin

   Revision 1.8  2003/05/17 21:25:16  oliver
   pgeom segfaulted so I reduced the dimensions of the arrays considerably and recompiled; also removed spurious EMPTY definition of SQR

   Revision 1.7  2002/08/27 02:20:49  oliver
   works but contains lots of dirt hacks in cylint.c...

   Revision 1.6  2002/08/25 22:19:52  oliver
   - substitued float with real; compiling with -DDOUBLE uses double precision for real, otherwise it is float
   - Boolean is now in the new xgtypes.h header (as is real)

   Revision 1.5  2002/08/22 02:23:54  oliver
   - added grid_alloc, rewritten for real ans using calloc instead of
     snew()
   - cleanup

   Revision 1.4  2000/12/13 16:40:52  oliver
   split functions in main (mgeom, pgeom), common (xg), and general (pg) -- stupid names...
   compile mgeom or pgeom in Makefile, use -DPROG
   all prototypes in separate files because my includes ar a crappy mess

   Revision 1.3  2000/12/07 13:04:21  oliver
   - writes coorect CRYST1 unitcell to pdb file
   - shifts coordinates to center of unitcell before writing pdb
   - new functions vec_add() and translate_pdb_atom()

   Revision 1.2  2000/12/06 17:15:32  oliver
   - added structure linking (*prev, *next) in domain
   - rewrote setup_domain():
     . global dz
     . removed stupid fixes
     also removed indentation of layers (rho1, rho2) in ini_domain()
   - rewrote setup_layer(): actually now the inner shape of the pore is
     continous but the code is still messy
   - linear(): removed z1==z2 fix

   Revision 1.1  2000/12/05 19:22:00  oliver
   Split pgeom seems to compile and to work (as buggy as usual...)


   utility routines, basic building blocks
*/

#ifndef UTIL_H_READ
#define UTIL_H_READ

#include <stdlib.h>
#include <stdio.h>
#include <float.h>
#include <math.h>
#include <stdarg.h>
#include "xgtypes.h"

#define PI 3.14159265358979323844
#define min(A, B) ((A) < (B) ? (A) : (B))
#define max(A, B) ((A) > (B) ? (A) : (B))

/* debug levels in mesg() (and few other places) */
#define OFF            0
#define WARN           5
#define INPUT         10
#define IMPORTANT     20
#define VERBOSE       40 
#define SUB1          60
#define SUB2          80
#define COORDINATES   90
#define INTERMEDIATES 95
#define ALL          100

/* 
   have a global variable indicating the current debug level to all functions 
*/
extern int debuglevel;

typedef struct {
  double x, y, z;
} Cartesian;

typedef struct {
  double rho, phi, z;
} Cylindrical;


extern void mesg (const int, const char *, ...);
extern void fatal_error (const int,  const char *, ...);
extern void piksrt (int, double []);
extern Cartesian new_vec  (double, double, double);
extern Cartesian vec_scalar (double, Cartesian);
extern Cartesian vec_diff (Cartesian, Cartesian);
extern Cartesian vec_add  (Cartesian, Cartesian);
extern double vec_dot (Cartesian, Cartesian);
extern double vec_cos (Cartesian, Cartesian);
extern double rad2deg (double);
extern double deg2rad (double);
extern Cartesian cyl2cart (Cylindrical);
extern Cylindrical cart2cyl (Cartesian);
/* extern int round (double); */ // is in standard math nowadays (or just on OS X?)
extern int sgn (double);
extern real ***grid3_alloc(int nx, int ny, int nz);
extern real **grid2_alloc(int nx, int ny);

#endif /* UTIL_H_READ */
