#pragma once
#ifndef COILBASE
typedef struct _POLYNOMINAL POLYNOMINAL;

//#include "protection/machine.h"
#ifdef SUN_OS_OLD
#define _POSIX_SOURCE
#endif

#include <math.h>
//#include <stdlib.h>
//#include "basic/define_real.h"

#endif
#define PI		3.14159265358979323846
#define	TwoPi	6.2831853071796
#define HalfPi	1.5707963267949 
#ifndef COILBASE
#define	DEG_to_RAD 0.0174532925199433
#define    PMU0    ( 4.e-7 * PI )
#define    PNU0    ( 1. / PMU0 ) 
#define    EPS0    8.854187817e-12
#define    LIGHT_V  2.99792458e8
#define    MY_INFINITY   999.
#define    MY_NODEFINE   1.e30 
#endif
#define SQRT(x) 	sqrt( (double) (x) )
#ifndef COILBASE
#define COS(x)  	cos( (double) (x) )
#define SIN(x)  	sin( (double) (x) )
#define EXP(x)  	exp( (double) (x) )
#define ABS  		abs
#define FABS(x) 	fabs( (double) (x) )
#define ATAN2(y, x)   	atan2( (double) (y), (double) (x))
#define ATOF  		atof
#define LOG10(x)        log10( (double) (x))
#endif
#define LOG(x)          log( (double) (x))
#ifndef COILBASE
#define FLOOR(x)          (int) floor( (double) (x))    

#define    NEARY_EQUAL(a, b)  (((b)==0.) ?  (FABS(a)<1.e-18) : ( FABS((a)/(b)-1.) < 1.e-6 ))

#endif
#define maximum(a,b)            (((a) > (b)) ? (a) : (b))
#define minimum(a,b)            (((a) < (b)) ? (a) : (b))
#ifndef COILBASE
struct _POLYNOMINAL
{
	unsigned int order;
	EM_REAL *c;
};

#ifdef  MOTION_GLOBAL_VALUE_MY_MATH
#define GLOBAL  
#define GLOBAL_VAL(v)      = (v)
#else
#define GLOBAL extern
#define GLOBAL_VAL(v)      /* */
#endif
GLOBAL POLYNOMINAL *create_polynominal(unsigned int);
GLOBAL EM_REAL polynominal_get_value(POLYNOMINAL *, EM_REAL);
GLOBAL EM_REAL polymominal_get_derivative(POLYNOMINAL *, EM_REAL);

//GLOBAL int nearly_equal_MY_POINTS(MY_POINT *p1, MY_POINT *p2);


#undef GLOBAL
#undef GLOBAL_VAL

#endif