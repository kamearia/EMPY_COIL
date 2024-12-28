#pragma once
#ifndef COILBASE
#include "basic/define_real.h"
#include "basic/my_math.h"
#include "EM_COIL.h"
#endif
#define ATAN4(PP1, PP2, PP3,PP4) ( atan2(PP4-PP2,1.+PP2*PP4) - atan2(PP3-PP1,1.+PP1*PP3) )
#ifndef COILBASE
/*#define PI 3.1415926535898*/
#endif
#define TwoDegree 0.0349065850398    /* 2 deg */
#define Myu0by4PI  1.e-7
#define ArcIntegralAccuracy 1.e-6

/*
struct field_sum{
	EM_REAL ARHO,ATHWAV,BZWAVE,BRHO,BTHETA;
};
*/
