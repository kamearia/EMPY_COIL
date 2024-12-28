#pragma once

#ifndef COILBASE
#define _USE_MATH_DEFINES
#include "basic/EM_Headers.h"
#include <memory.h>
#include <cmath>
#endif
#include "basic/EM_Math.h"
#ifndef COILBASE
const EM_REAL EM_Constant::Pi        = M_PI;					//! â~é¸ó¶
const EM_REAL EM_Constant::DOUBLE_PI = 2. * EM_Constant::Pi;	//! â~é¸ó¶ÇÃ2î{
const EM_REAL EM_Constant::MU0       = 4.e-7 * M_PI;			//! ê^ãÛÇÃìßé•ó¶
const EM_REAL EM_Constant::C         = 2.99792458e8;			//! ê^ãÛÇÃåıë¨
const EM_REAL EM_Constant::EPSILON0  = 1. / EM_Constant::MU0 / EM_Constant::C / EM_Constant::C;		//! ê^ãÛÇÃóUìdó¶

/*
template <class T> 
bool EM_Permute( ON_SimpleArray<T> & array,const int* index ){
  bool rc = false;
  int  m_count=array.Count();
  T *m_a=array.Array();
  if ( m_a && m_count > 0 && index ) {
    int i;
    T* buffer = (T*)onmalloc(m_count*sizeof(buffer[0]));
    memcpy( buffer, m_a, m_count*sizeof(T) );
    for (i = 0; i < m_count; i++ )
      memcpy( m_a+i, buffer+index[i], sizeof(T) ); // must use memcopy and not operator=
    onfree(buffer);
    rc = true;
  }
  return rc;
}
*/
/*
EM_REAL Conjugate(const EM_REAL &c){
	return c;
}
COMPLEX Conjugate(const COMPLEX &c){
	return COMPLEX(c.real(), -c.imag());
}
*/

bool EM_IsValid(double x)
{
  return (x != EM_UNSET_VALUE && EM_IsFinite(x) );
}
#endif
#include <float.h>
bool EM_IsFinite(double x)
{
//#ifdef WIN95
  return (_finite(x)?true:false);
//#else
 // return (std::isfinite(x)?true:false);
//#endif

}


int EM_GetNumDigits(int  n)
{
	if (n < 0)
		n *= -1;
	int ndigit = 1;
	int powered = 10;
	while (powered <= n) {
		ndigit++;
		powered *= 10;
	}

	return ndigit;
}
#ifndef COILBASE

int BinarySearchInt(const EM_SimpleArray<int> &data, int first, int last, int &search_key)
{
	int index;

	if (first > last)
		index = -1;

	else
	{
		int mid = (first + last) / 2;

		if (search_key == data[mid])
			index = mid;
		else

			if (search_key < data[mid])
				index = BinarySearchInt(data, first, mid - 1, search_key);
			else
				index = BinarySearchInt(data, mid + 1, last, search_key);

	} 
	return index;
}


int
EM_Solve2x2( double m00, double m01, double m10, double m11, double d0, double d1,
                  double* x_addr, double* y_addr, double* pivot_ratio)
/* Solve a 2x2 system of linear equations
 *
 * INPUT:
 *   m00, m01, m10, m11, d0, d1
 *      coefficients for the 2x2 the linear system:
 *   x_addr, y_addr
 *      addresses of doubles
 *   pivot_ratio
 *      address of double
 * OUTPUT:
 *   ON_Solve2x2() returns rank (0,1,2)
 *
 *   If ON_Solve2x2() is successful (return code 2), then
 *   the solution is returned in {*x_addr, *y_addr} and
 *   *pivot_ratio = min(|pivots|)/max(|pivots|).
 *
 *   WARNING: If the pivot ratio is small, then the matrix may
 *   be singular or ill conditioned.  You should test the results
 *   before you use them.
 *
 * COMMENTS:
 *      The system of 2 equations and 2 unknowns (x,y),
 *         m00*x + m01*y = d0
 *         m10*x + m11*y = d1,
 *      is solved using Gauss-Jordan elimination
 *      with full pivoting.
 * EXAMPLE:
 *      // Find the intersection of 2 2D lines where
 *      // P0, P1  are points on the lines and
 *      // D0, D1, are nonzero directions
 *      rc = ON_Solve2x2(D0[0],-D1[0],D0[1],-D1[1],P1[0]-P0[0],P1[1]-P0[1],
 *                       &x, &y,&pivot_ratio);
 *      switch(rc) {
 *      case  0: // P0 + x*D0 = P1 + y*D1 = intersection point
 *        if (pivot_ratio < 0.001) {
 *          // small pivot ratio - test answer before using ...
 *        }
 *        break;
 *      case -1: // both directions are zero!
 *        break;
 *      case -2: // parallel directions
 *        break;
 *      }
 *
 * REFERENCE:
 *      STRANG
 *
 * RELATED FUNCTIONS:
 *      ON_Solve3x2(), ON_Solve3x3
 */
{
  int i = 0;
  double maxpiv, minpiv;
  double x = fabs(m00);
  double y = fabs(m01); if (y > x) {x = y; i = 1;}
  y = fabs(m10); if (y > x) {x = y; i = 2;}
  y = fabs(m11); if (y > x) {x = y; i = 3;}
  *pivot_ratio = *x_addr = *y_addr = 0.0;
  if (x == 0.0) 
    return 0; // rank = 0
  minpiv = maxpiv = x;
  if (i%2) {
    {double* tmp = x_addr; x_addr = y_addr; y_addr = tmp;}
    x = m00; m00 = m01; m01 = x;
    x = m10; m10 = m11; m11 = x;
  }
  if (i > 1) {
    x = d0; d0 = d1; d1 = x;
    x = m00; m00 = m10; m10 = x;
    x = m01; m01 = m11; m11 = x;
  }
  
  x = 1.0/m00;
  m01 *= x; d0 *= x;
  if (m10 != 0.0) {m11 -= m10*m01; d1 -= m10*d0;}

  if (m11 == 0.0) 
    return 1; // rank = 1

  y = fabs(m11);
  if (y > maxpiv) maxpiv = y; else if (y < minpiv) minpiv = y;
  
  d1 /= m11;
  if (m01 != 0.0)
    d0 -= m01*d1;

  *x_addr = d0;
  *y_addr = d1;
  *pivot_ratio = minpiv/maxpiv;
  return 2;  
}

int
EM_Solve3x3(const double row0[3], const double row1[3], const double row2[3],
                double d0, double d1, double d2,
                double* x_addr, double* y_addr, double* z_addr,
                double* pivot_ratio)
{
  /* Solve a 3x3 linear system using Gauss-Jordan elimination 
   * with full pivoting.
   */
  int i, j;
  double* p0;
  double* p1;
  double* p2;
  double x, y, workarray[12], maxpiv, minpiv;

  const int sizeof_row = 3*sizeof(row0[0]);

  *pivot_ratio = *x_addr = *y_addr = *z_addr = 0.0;
  x = fabs(row0[0]); i=j=0;
  y = fabs(row0[1]); if (y>x) {x=y;j=1;}
  y = fabs(row0[2]); if (y>x) {x=y;j=2;}
  y = fabs(row1[0]); if (y>x) {x=y;i=1;j=0;}
  y = fabs(row1[1]); if (y>x) {x=y;i=1;j=1;}
  y = fabs(row1[2]); if (y>x) {x=y;i=1;j=2;}
  y = fabs(row2[0]); if (y>x) {x=y;i=2;j=0;}
  y = fabs(row2[1]); if (y>x) {x=y;i=2;j=1;}
  y = fabs(row2[2]); if (y>x) {x=y;i=2;j=2;}
  if (x == 0.0) 
    return 0;
  maxpiv = minpiv = fabs(x);
  p0 = workarray;
  switch(i) {
  case 1: /* swap rows 0 and 1 */
    memcpy(p0,row1,sizeof_row); p0[3] = d1; p0 += 4;
    memcpy(p0,row0,sizeof_row); p0[3] = d0; p0 += 4;
    memcpy(p0,row2,sizeof_row); p0[3] = d2;
    break;
  case 2: /* swap rows 0 and 2 */
    memcpy(p0,row2,sizeof_row); p0[3] = d2; p0 += 4;
    memcpy(p0,row1,sizeof_row); p0[3] = d1; p0 += 4;
    memcpy(p0,row0,sizeof_row); p0[3] = d0;
    break;
  default:
    memcpy(p0,row0,sizeof_row); p0[3] = d0; p0 += 4;
    memcpy(p0,row1,sizeof_row); p0[3] = d1; p0 += 4;
    memcpy(p0,row2,sizeof_row); p0[3] = d2;
    break;
  }
  switch(j) {
  case 1: /* swap columns 0 and 1 */
    p0 = x_addr; x_addr = y_addr; y_addr = p0;
    p0 = &workarray[0]; 
    x = p0[0]; p0[0]=p0[1]; p0[1]=x; p0 += 4;
    x = p0[0]; p0[0]=p0[1]; p0[1]=x; p0 += 4;
    x = p0[0]; p0[0]=p0[1]; p0[1]=x;
    break;
  case 2: /* swap columns 0 and 2 */
    p0 = x_addr; x_addr = z_addr; z_addr = p0;
    p0 = &workarray[0]; 
    x = p0[0]; p0[0]=p0[2]; p0[2]=x; p0 += 4;
    x = p0[0]; p0[0]=p0[2]; p0[2]=x; p0 += 4;
    x = p0[0]; p0[0]=p0[2]; p0[2]=x;
    break;
  }

  x = 1.0/workarray[0];
  /* debugger set workarray[0] = 1 */
  p0 = p1 = workarray + 1;
  *p1++ *= x; *p1++ *= x; *p1++ *= x;
  x = -(*p1++);
  /* debugger set workarray[4] = 0 */
  if (x == 0.0) 
    p1 += 3;
  else 
    {*p1++ += x*(*p0++); *p1++ += x*(*p0++); *p1++ += x*(*p0); p0 -= 2;}
  x = -(*p1++);
  /* debugger set workarray[8] = 0 */
  if (x != 0.0)
    {*p1++ += x*(*p0++); *p1++ += x*(*p0++); *p1++ += x*(*p0); p0 -= 2;}

  x = fabs(workarray[ 5]);i=j=0;
  y = fabs(workarray[ 6]);if (y>x) {x=y;j=1;}
  y = fabs(workarray[ 9]);if (y>x) {x=y;i=1;j=0;}
  y = fabs(workarray[10]);if (y>x) {x=y;i=j=1;}
  if (x == 0.0) 
    return 1; // rank = 1;
  y = fabs(x);
  if (y > maxpiv) maxpiv = y; else if (y < minpiv) minpiv = y;
  if (j) {
    /* swap columns 1 and 2 */
    p0 = workarray+1;
    p1 = p0+1;
    x = *p0; *p0 = *p1; *p1 = x; p0 += 4; p1 += 4;
    x = *p0; *p0 = *p1; *p1 = x; p0 += 4; p1 += 4;
    x = *p0; *p0 = *p1; *p1 = x; p0 += 4; p1 += 4;
    p0 = y_addr; y_addr = z_addr; z_addr = p0;
  }

  if (i) {
    /* pivot is in row 2 */
    p0 = workarray+1;
    p1 = p0 + 8;
    p2 = p0 + 4;
  }
  else {
    /* pivot is in row 1 */
    p0 = workarray+1;
    p1 = p0 + 4;
    p2 = p0 + 8;
  }

  /* debugger set workarray[5+4*i] = 1 */
  x = 1.0/(*p1++); *p1++ *= x; *p1 *= x; p1--;
  x = -(*p0++);
  /* debugger set p0[-1] = 0 */
  if (x != 0.0) {*p0++ += x*(*p1++); *p0 += x*(*p1); p0--; p1--;}
  x = -(*p2++);
  /* debugger set p2[-1] = 0 */
  if (x != 0.0) {*p2++ += x*(*p1++); *p2 += x*(*p1); p2--; p1--;}
  x = *p2++;
  if (x == 0.0) 
    return 2; // rank = 2;
  y = fabs(x);
  if (y > maxpiv) maxpiv = y; else if (y < minpiv) minpiv = y;
  /* debugger set p2[-1] = 1 */
  *p2 /= x;
  x = -(*p1++);  if (x != 0.0) *p1 += x*(*p2);
  /* debugger set p1[-1] = 0 */
  x = -(*p0++);  if (x != 0.0) *p0 += x*(*p2);
  /* debugger set p0[-1] = 0 */
  *x_addr = workarray[3];
  if (i) {
    *y_addr = workarray[11];
    *z_addr = workarray[7];
  }
  else {
    *y_addr = workarray[7];
    *z_addr = workarray[11];
  }
  *pivot_ratio = minpiv/maxpiv;
  return 3;
}

double EM_Rand(){
	static double max=EM_RAND_MAX, min=EM_RAND_MIN;
	return ((double)rand()/(double)RAND_MAX)*(max-min)+min;
}

EM_MathFunction::EM_MathFunction(int dim):dim(dim){
/*
	f=new double[dim];
	df=new double*[dim];
	for(int i=0; i<dim; ++i)
		df[i]=new double[dim];
*/
}

EM_MathFunction::~EM_MathFunction(){
/*
	delete[] f;
	for(int i=0; i<dim; ++i) delete[] df[i];
	delete[] df;
*/
}

int EM_MathFunction::Solve(double *x, double eps){
	assert(eps>0.);
	double *y=new double[dim];
	double error=eps*10.;
	int k=0;
	while(error>eps){
		++k;
		GetValues(x, f, df); 
		LinearSolve(y, f, df);
		error=0.;
		for(int i=0; i<dim; ++i){
			double tmp = std::abs(y[i]);
			if(tmp>error) error=tmp;
			x[i] += y[i];
		}
	}
	delete[] y;
	return k;
}



#endif