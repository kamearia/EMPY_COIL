#pragma once
#ifndef COILBASE
#include "parallel/EM_Parallel.h"
#include <limits.h>
#include <float.h>
#endif
#include <complex>
#ifndef COILBASE
#include <assert.h>

#include "protection/user.h"
#include "protection/machine.h"
#include "basic/my_math.h"
#include <stdio.h>

//#define DEBUG_WRITE
#include "basic/EM_DebugWriter.h"
/*
#define SINGLE
*/
//#define DOUBLE
#endif
#ifdef SINGLE
#define EM_REAL float
#else
#define EM_REAL double
#endif
#ifndef COILBASE
#endif
typedef std::complex<EM_REAL> COMPLEX;
#ifndef CODEBASE
#define EM_CLASS
#define EM_DECL
#define EM_EXTERN_DECL

#define EM_PI 3.141592653589793238462643

#define EM_DEGREES_TO_RADIANS EM_PI/180.0
#define EM_RADIANS_TO_DEGREES 180.0/EM_PI
#if !defined(DBL_EPSILON)
#define DBL_EPSILON     2.2204460492503131e-016 /* smallest such that 1.0+DBL_EPSILON != 1.0 */
#endif 
#define EM_ZERO_TOLERANCE  1.0e-12

#define EM_INT_MAX 2147483647
#define EM_INT_MIN -2147483648
#define EM_DBL_MAX		DBL_MAX  // 1.7976931348623158e+308
#define EM_UNSET_VALUE	-1.23432101234321e+308
#define EM_DBL_MIN		DBL_MIN
#define EM_SQRT_EPSILON 1.490116119385000000e-8
#define EM_DEFAULT_ANGLE_TOLERANCE (EM_PI/180.0)
typedef unsigned int UINT;
typedef unsigned short USINT;
#define EM_UINT_UNSET_VALUE UINT_MAX//4,294,967,295
typedef unsigned short int USHRT;  
#define EM_USINT_UNSET_VALUE USHRT_MAX  //65,535
typedef unsigned char UCHAR;
#define EM_UCHAR_UNSET_VALUE UCHAR_MAX  //255

//#define POST_REAL_TIME

#if defined(DBL_EPSILON)
#define EM_EPSILON DBL_EPSILON
#else
#define EM_EPSILON 2.2204460492503131e-16
#endif

#define EM_3fPoint ON_3fPoint
#define EM_TextLog ON_TextLog
#define EM_Workspace ON_Workspace

#define emmalloc onmalloc
#define emfree onfree

template < class U, class V >
struct SameType{
	enum { result =false};
};

template < class T>
struct SameType<T, T>{
	enum { result =true};
};

const unsigned char INVALID_UCHAR=UCHAR_MAX;

typedef int INDEX;
const int INVALID_INDEX=-1;
struct SIGNED_INDEX{
	SIGNED_INDEX(){}
	SIGNED_INDEX(INDEX index, bool sign): index(index), sign(sign){}
	INDEX Index() const { return index;}
	void DebugWrite(FILE *file) const { fprintf(file, "    %d  (%d)\n", index, sign); }
	INDEX index;
	bool sign;
	EM_REAL Dir() const { return sign ? 1.0 : -1.0; }
	static int CompareIndex(const SIGNED_INDEX *a, const SIGNED_INDEX *b) {
		int ia = a->index;
		int ib = b->index;
		if (ia < ib) return -1;
		else if (ia> ib) return 1;
		else return 0;
	}

};

inline INDEX Index(INDEX i) { return i; }
inline bool Sign(INDEX i) { return true; }
inline INDEX Index(SIGNED_INDEX i) { return i.index; }
inline bool Sign(SIGNED_INDEX i) { return i.sign; }

typedef int DOMAIN_INDEX;
const DOMAIN_INDEX INVALID_DOMAIN_INDEX=-1;

typedef INDEX EM_ID;
const EM_ID INVALID_ID=INVALID_INDEX;

//typedef INDEX EM_SIZE_T;
typedef size_t EM_SIZE_T;

#define EM_DELETE(p) if(p){ delete p; p=0; }
#define EM_DELETE_ARRAY(p) if(p){ delete[] p; p=0; }

typedef unsigned int DIRECTION;
enum FACE_DIRECTION { FACE_INWARD=0, FACE_OUTWARD };
enum PositionDependency { POSITION_INDEPEND, POSITION_DEPEND};
enum PopertyDistributionType { REGIONWISE=0, ELEMENTWISE
};

typedef size_t EM_size_t;
typedef long long int EM_long_int;

template <class T>
class EM_Singleton{
protected:
	EM_Singleton(){}
	~EM_Singleton(){}
	EM_Singleton(const EM_Singleton &x){}
	EM_Singleton &operator=(const EM_Singleton &) {return *this; };
//	static T obj;
public:
	static T *p() { static T obj; return &obj; }
//	static T *p()  { return &obj; }
};



//template <class T>
//T EM_Singleton<T>::obj;

#endif



