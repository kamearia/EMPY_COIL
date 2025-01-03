#pragma once

#define PI		3.14159265358979323846
#define	TwoPi	6.2831853071796
#define HalfPi	1.5707963267949 
#define	DEG_to_RAD 0.0174532925199433
#define    PMU0    ( 4.e-7 * PI )
#define    PNU0    ( 1. / PMU0 ) 
#define    EPS0    8.854187817e-12
#define    LIGHT_V  2.99792458e8
#define    MY_INFINITY   999.
#define    MY_NODEFINE   1.e30 
#define SQRT(x) 	sqrt( (double) (x) )
#define COS(x)  	cos( (double) (x) )
#define SIN(x)  	sin( (double) (x) )
#define EXP(x)  	exp( (double) (x) )
//#define ABS  		abs
#define FABS(x) 	fabs( (double) (x) )
#define ATAN2(y, x)   	atan2( (double) (y), (double) (x))
#define ATOF  		atof
#define LOG10(x)        log10( (double) (x))
#define LOG(x)          log( (double) (x))
#define FLOOR(x)          (int) floor( (double) (x))    
#define    NEARY_EQUAL(a, b)  (((b)==0.) ?  (FABS(a)<1.e-18) : ( FABS((a)/(b)-1.) < 1.e-6 ))
#define maximum(a,b)            (((a) > (b)) ? (a) : (b))
#define minimum(a,b)            (((a) < (b)) ? (a) : (b))

#define EM_MAX(a,b) ((a)>(b)? (a):(b))
#define EM_MIN(a,b) ((a)<(b)? (a):(b))

bool EM_IsFinite(double x);
