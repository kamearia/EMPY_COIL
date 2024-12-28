#pragma once
#ifndef COILBASE
#include <math.h>
#include <stdlib.h>
#include "basic/my_math.h"
#ifndef COILBASE
#include "basic/EM_SimpleArray.h"
#endif
#endif
#define EM_MAX(a,b) ((a)>(b)? (a):(b))
#define EM_MIN(a,b) ((a)<(b)? (a):(b))
#ifndef COILBASE
namespace EM_Constant
{
	extern const EM_REAL Pi;					//! ‰~Žü—¦
	extern const EM_REAL DOUBLE_PI;				//! ‰~Žü—¦‚Ì2”{
	extern const EM_REAL MU0;					//! ^‹ó‚Ì“§Ž¥—¦
	extern const EM_REAL C;						//! ^‹ó‚ÌŒõ‘¬
	extern const EM_REAL EPSILON0;				//! ^‹ó‚Ì—U“d—¦
}

template <typename T>
T Conjugate(const T &c);

template <typename T>
EM_REAL Norm2(const T &c);

template <>
inline EM_REAL Conjugate(const EM_REAL &c) {
	return c;
}
template <>
inline EM_REAL Norm2(const EM_REAL &c) {
	return c*c;
}


template <>
inline COMPLEX Conjugate(const COMPLEX &c) {
	return COMPLEX(c.real(), -c.imag());
}

template <>
inline EM_REAL Norm2(const COMPLEX &c) {
	EM_REAL cr = c.real();
	EM_REAL ci = c.imag();
	return cr*cr + ci*ci;
}


template <class T>
int EM_CompareValue(const  T * a, const T * b)
{
	return *a>*b ?
		1 : *a == *b ? 0 : -1;
}

template <typename T>
EM_REAL RealValue(T a) { return a.real(); }

template <>
inline EM_REAL RealValue(EM_REAL a) { return a; }


/*
EM_REAL Conjugate(const EM_REAL &c){
return c;
}
COMPLEX Conjugate(const COMPLEX &c){
return COMPLEX(c.real(), -c.imag());
}
*/
/*
template <class T>
bool EM_Permute( EM_SimpleArray<T> & array,const int* index ){
bool rc = false;
int  m_count=array.Count();
T *m_a=array.Array();
if ( m_a && m_count > 0 && index ) {
int i;
int * buffer = (int*)emmalloc(m_count*sizeof(int));
memcpy( buffer, index, m_count*sizeof(int) );
for (i = 0; i < m_count; i++ ){
int to=i;
int from=buffer[to];
T tmp;
if(to !=from){
memcpy( &tmp, m_a+to, sizeof(T) );
while(from !=i){
memcpy( m_a+to, m_a+from, sizeof(T) );
to=from;
from=buffer[to];
buffer[to]=to;
};
memcpy( m_a+to, &tmp, sizeof(T) );
buffer[to]=to;
}
}
emfree(buffer);
rc = true;
}
return rc;
}
*/

bool EM_IsValid(double x);

int EM_GetNumDigits(int  n);
#endif
bool EM_IsFinite(double x);
#ifndef COILBASE
int BinarySearchInt(const EM_SimpleArray<int> &data, int first, int last, int &search_key);

// returns rank - if rank != 2, system is under determined
// If rank = 2, then solution to 
//
//          a00*x0 + a01*x1 = b0, 
//          a10*x0 + a11*x1 = b1 
//
// is returned
EM_DECL
int EM_Solve2x2(
	double, double,   // a00 a01 = first row of 2x2 matrix
	double, double,   // a10 a11 = second row of 2x2 matrix
	double, double,   // b0 b1
	double*, double*, // x0, x1 if not NULL, then solution is returned here
	double*           // if not NULL, then pivot_ratio returned here
);

/*
Description:
Use Gauss-Jordan elimination with full pivoting to solve
a system of 3 linear equations and 3 unknowns(x,y,z)

x*row0[0] + y*row0[1] + z*row0[2] = d0
x*row1[0] + y*row1[1] + z*row1[2] = d1
x*row2[0] + y*row2[1] + z*row2[2] = d2

Parameters:
row0 - [in] first row of 3x3 matrix
row1 - [in] second row of 3x3 matrix
row2 - [in] third row of 3x3 matrix
d0 - [in]
d1 - [in]
d2 - [in] (d0,d1,d2) right hand column of system
x_addr - [in] first unknown
y_addr - [in] second unknown
z_addr - [in] third unknown
pivot_ratio - [out] if not NULL, the pivot ration is
returned here.  If the pivot ratio is "small",
then the matrix may be singular or ill
conditioned. You should test the results
before you use them.  "Small" depends on the
precision of the input coefficients and the
use of the solution.  If you can't figure out
what "small" means in your case, then you
must check the solution before you use it.

Returns:
The rank of the 3x3 matrix (0,1,2, or 3)
If EM_Solve3x3() is successful (returns 3), then
the solution is returned in
(*x_addr, *y_addr, *z_addr)
and *pivot_ratio = min(|pivots|)/max(|pivots|).
If the return code is < 3, then (0,0,0) is returned
as the "solution".

See Also:
EM_Solve2x2
EM_Solve3x2
EM_Solve4x4
*/
EM_DECL
int EM_Solve3x3(
	const double row0[3],
	const double row1[3],
	const double row2[3],
	double d0,
	double d1,
	double d2,
	double* x_addr,
	double* y_addr,
	double* z_addr,
	double* pivot_ratio
);

// Random double Value -EM_RAND_MIN <-> EM_RAND_MAX.
//   Initially call srand( (unsigned) seed) ).
#define EM_RAND_MIN -1.e-6
#define EM_RAND_MAX 1.e-6
EM_DECL
double EM_Rand();


class EM_MathFunction {
public:
	EM_MathFunction(int dim);
	EM_MathFunction(int dim, double *f, double **df) :
		dim(dim), f(f), df(df) {}
	~EM_MathFunction();
	int GetDimension() const { return dim; }
	//	virtual double *Func(const double *x)=0;
	//	virtual double **DerivFunc(const double *x)=0;
	virtual void GetValues(const double *x, double *f, double **df) = 0;

	/* Solve Nonlinear Equation
	input x: initial solution
	eps: absolute error criterion
	output: x:result
	return iteration numbers
	*/
	virtual int Solve(double *x, double eps);
	virtual void LinearSolve(double *x, double *f, double **df) = 0;

protected:
	int dim;
	double *f;
	double **df;
};
#endif