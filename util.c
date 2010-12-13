/* 
   $Id: util.c,v 1.9 2008/02/18 01:18:48 oliver Exp $
   $Log: util.c,v $
   Revision 1.9  2008/02/18 01:18:48  oliver
   o Can read geometry from input file:
     - added datafile_geom()
     - updated usage
   o Set MAXDOMAINS to 6 and MAXRINGS to 4 (but cannot be much more otherwise binary
     crashes with segfault/buserror). This should be dynamically allocated
     anyway...
   o use fatal_error() and mesg() instead of printf/exit constructs
   o started cleaning up Log and whitespace

   Revision 1.8  2008/02/17 16:55:54  oliver
   - fixed DOUBLE compilation
   - default BIN_DIR to ~/bin

   Revision 1.7  2002/08/27 02:20:49  oliver
   works but contains lots of dirt hacks in cylint.c...

   Revision 1.6  2002/08/25 22:19:52  oliver
   - substitued float with real; compiling with -DDOUBLE uses double precision for real, otherwise it is float
   - Boolean is now in the new xgtypes.h header (as is real)

   Revision 1.5  2002/08/22 02:23:54  oliver
   - added grid_alloc, rewritten for float ans using calloc instead of
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


   utility functions, basic mathematical stuff
*/

#include "util.h"

Cartesian cyl2cart (Cylindrical u)
{
  Cartesian t;

  t.x = u.rho * cos(u.phi);
  t.y = u.rho * sin(u.phi);
  t.z = u.z;

  return t;
};

Cylindrical cart2cyl (Cartesian t)
{
  Cylindrical u;

  u.rho = sqrt(t.x*t.x + t.y*t.y);
  u.phi = atan2(t.y,t.x);
  u.z   = t.z;

  return u;
};

Cartesian new_vec (double x, double y, double z) {
  Cartesian this;

  this.x = x;
  this.y = y;
  this.z = z;

  return this;
};


Cartesian vec_diff (Cartesian a, Cartesian b) {
  Cartesian this;
  
  this.x = a.x - b.x;
  this.y = a.y - b.y;
  this.z = a.z - b.z;
  
  return this;
};

Cartesian vec_add (Cartesian a, Cartesian b) {
  Cartesian this;
  
  this.x = a.x + b.x;
  this.y = a.y + b.y;
  this.z = a.z + b.z;
  
  return this;
};

Cartesian vec_scalar (double a, Cartesian u) {
  Cartesian this;
  
  this.x = a * u.x;
  this.y = a * u.y;
  this.z = a * u.z;
  
  return this;
};


double vec_dot (Cartesian a, Cartesian b) {
  return a.x*b.x + a.y*b.y + a.z*b.z;
};

double vec_cos (Cartesian a, Cartesian b) {
  return vec_dot (a, b) / sqrt (vec_dot(a,a) * vec_dot(b,b));
};

double deg2rad (double alpha) {
  return alpha * PI/180.0;
}

double rad2deg (double x) {
  return x * 180.0/PI;
}


void piksrt (int n, double arr[]) {
  /* numerical recipes, p330 */
  int    i, j;
  double a;

  for (j = 1; j < n; j++) {
    a = arr[j];
    i = j - 1;
    while (i >= 0 && arr[i] > a) {
      arr[i + 1] = arr[i];
      i--;
    };
    arr[i + 1] = a;
  };
};


void mesg (const int threshold, const char *format, ...)
{
  va_list args;
  va_start(args, format);
  
  if (debuglevel >= threshold) {
    vprintf(format, args);
    printf("\n");
  };
  va_end(args);
  return;
};

void fatal_error (const int err,  const char *format, ...)
{
  char fmt[1024];
  va_list args;
  va_start(args, format);

  sprintf(fmt,"FATAL (%d): %s", err, format);
  vprintf(fmt, args);
  abort();
  va_end(args);
  return;
};

/* is in math nowadays (or just on OS X ?
int round (double x) {
  return (int) (x + sgn(x) * 0.5);
}
*/

int sgn (double x) {
  if (x > 0) {
    return 1;
  } else if (x == 0) {
    return 0;
  } else {
    return -1;
  };
};
  
real ***grid3_alloc(int nx, int ny, int nz) {
  /* allocate a 3D array, size nx*ny*nz which is referenced through pointers.
     By using snew() the memory is already 0-initialised.
     range is g[0..nx-1,0..ny-1,0..nz-1].
     x: rows
     y: columns
     z: depth
     This is a trimmed down version of "Numerical Recipes in C", p944 , f3tensor() 

     This memory/pointer stuff is mind boggling: Draw it on
     paper. Basically, g[] is a vector, pointing to the rows of an
     array g[][], nx*ny (allocated as one block). These entries again
     point to rows in the 'float' chunk of memory nx*ny*nz, g[][][].

     Because we have allocate the memory in blocks we can initialise
     all pointers in a fairly tricky (to me, anyway) fashion.

     The original routine allows for another memory range of the
     indices as well.  

  */
  int  i,j;
  real ***g = NULL;

  /* allocate a row of 'ptr to ptr to rows' (= index by x values), length nx */
  g=(real ***)calloc(nx,sizeof(real **));
  if (!g) fatal_error(-1,"grid3_alloc(): Cannot allocate memory for nx=%d\n",nx);
  
  /* allocate nx 'ptrs to rows' (length ny each) and set ptrs to them */
  g[0]=(real **)calloc(nx*ny,sizeof(real *));
  if (!g[0]) fatal_error(-1,"grid3_alloc(): Cannot allocate memory for nx*ny=%d*%d\n",
		      nx,ny);
    
  /* allocate nx*ny 'rows' (length nz each; these memory locations hold the floats) */
  g[0][0]=(real *)calloc(nx*ny*nz,sizeof(real));
  if (!g[0][0]) fatal_error(-1,"grid3_alloc(): Cannot allocate memory for data "
		      "(nx*ny*nz=%d*%d%d)\n",nx,ny,nz);

  /* initialise all pointers (starting from the already initialised
     initial locations), "leaping" through the array 
  */
  for(j=1; j<ny; j++) g[0][j]=g[0][j-1]+nz;
  for(i=1; i<nx; i++) {
    g[i]=g[i-1]+ny;
    g[i][0]=g[i-1][0]+ny*nz;
    for(j=1; j<ny; j++) g[i][j]=g[i][j-1]+nz;
  }

  /* return ptr to an array of ptr to rows */
  return g;
}

real **grid2_alloc(int nx, int ny) {
  /* allocate a 3D array, size nx*ny*nz which is referenced through pointers.
     By using snew() the memory is already 0-initialised.
     range is g[0..nx-1,0..ny-1].
     x: rows
     y: columns
     This is a trimmed down version of "Numerical Recipes in C", p944 , f3tensor() 
  */
  int  i;
  real **g = NULL;

  /* allocate a row of 'ptr to ptr to rows' (= index by x values), length nx */
  g=(real **)calloc(nx,sizeof(real *));
  if (!g) fatal_error(-1,"grid2_alloc(): Cannot allocate memory for nx=%d\n",nx);
  
  /* allocate nx 'ptrs to rows' (length ny each) and set ptrs to them */
  g[0]=(real *)calloc(nx*ny,sizeof(real));
  if (!g[0]) fatal_error(-1,"grid2_alloc(): Cannot allocate memory for nx*ny=%d*%d\n",
		      nx,ny);

  for(i=1;i<nx;i++) g[i]=g[i-1]+ny;
  
  /* return ptr to an array of ptr to rows */
  return g;
}
