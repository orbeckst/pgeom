/* 	$Id: xgtypes.h,v 1.2 2002/08/27 02:20:51 oliver Exp $	 
	$Log: xgtypes.h,v $
	Revision 1.2  2002/08/27 02:20:51  oliver
	works but contains lots of dirt hacks in cylint.c...
	
	Revision 1.1  2002/08/25 22:19:52  oliver
	- substitued float with real; compiling with -DDOUBLE uses double precision for real, otherwise it is float
	- Boolean is now in the new xgtypes.h header (as is real)
	

	common types
*/
#ifndef XGTYPES_H
#define XGTYPES_H

#include <float.h>
#include <math.h>

#ifndef HAVE_REAL
#define HAVE_REAL
#ifdef DOUBLE
typedef double real;
#define FTOL   sqrt(DBL_EPSILON)
         /* tolerance for minimisation; cannot do better than
            sqrt(machine precision) (see Numerical Recipes) 
	    This is ca 1e-8 for DBL, and 1e-4 for FLOAT
	 */
#else
typedef float real;
#define FTOL   sqrt(FLT_EPSILON)
#endif /* DOUBLE    */
#endif /* HAVE_REAL */

typedef real rvec[3];

#define DIM 3
/* cartesian */
#define XX 0
#define YY 1
#define ZZ 2

/* cylindrical */
#define RAD 0
#define PHI 1
#define ZZZ 2

typedef enum {FALSE, TRUE} Boolean;
#endif /* XGTYPES_H */
