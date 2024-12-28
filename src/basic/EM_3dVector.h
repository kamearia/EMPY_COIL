////////////////////////////////////////////////////////////////
//
//   defines double precision vector in 3-dimension
//
////////////////////////////////////////////////////////////////
#pragma once

////////////////////////////////////////////////////////////////
//
//   EM_3dVector
//
#include <assert.h>
#include "basic/EM_defines.h"
#include "basic/EM_Math.h"
class EM_3dPoint;
template <typename T>
class EM_3dVectorBase/*:public ON_3dVector*/
{
public:
	T x[3];
	inline EM_3dVectorBase();//:ON_3dVector(){}
	inline ~EM_3dVectorBase(){}
	inline EM_3dVectorBase( const float* v ); //:ON_3dVector(v){}
	inline EM_3dVectorBase( const T* v );//:ON_3dVector(v){}
	inline EM_3dVectorBase(T x, T y, T z);//:ON_3dVector(x,y,z){}
//	inline EM_3dVectorBase( const EM_3dVectorBase &v );//:ON_3dVector(v){}

	inline EM_3dVectorBase( const EM_3dVectorBase<EM_REAL> &v );//:ON_3dVector(v){}
	inline EM_3dVectorBase( const EM_3dVectorBase<COMPLEX> &v );//:ON_3dVector(v){}
	inline EM_3dVectorBase( const EM_3dPoint &p );
  // (double*) conversion operators
	inline operator T*();
	inline operator const T*() const;
	inline T operator*(const EM_3dVectorBase&) const; // inner (dot) product

	inline EM_3dVectorBase  operator-() const;
	inline EM_3dVectorBase& operator*=(T);
	inline EM_3dVectorBase& operator/=(double);
	inline EM_3dVectorBase& operator+=(const EM_3dVectorBase&);
	inline EM_3dVectorBase& operator-=(const EM_3dVectorBase&);
	inline EM_3dVectorBase  operator*(const T & ) const;
	inline EM_3dVectorBase  operator/(double) const;
	inline EM_3dVectorBase  operator+(const EM_3dVectorBase&) const;
	inline EM_3dVectorBase  operator-(const EM_3dVectorBase&) const;
  // index operators mimic T[3] behavior
	inline T& operator[](int);
	inline const T& operator[](int) const;
	inline void Set(double x,double y,double z){
	  this->x[0]=x; this->x[1]=y; this->x[2]=z; 
	}
  /*
  Description
    If any coordinate of a point is ON_UNSET_VALUE,
    then the point is not valid.
  Returns:
    true if the point is valid.
  */
  bool IsValid() const;

  double Length() const;
  double Length2() const;

  // set this vector to be perpendicular to another vector
  bool PerpendicularTo( // Result is not unitized. 
                        // returns false if input vector is zero
        const EM_3dVectorBase& 
        );
  // set this vector to be perpendicular to a plane defined by 3 points
  bool PerpendicularTo(
               // about 3 times slower than
               //    EM_3dVectorBase N = EM_CrossProduct(P1-P0,P2-P0); 
               //    N.Unitize();
               // returns false if points are coincident or colinear
         const EM_3dPoint&, const EM_3dPoint&, const EM_3dPoint& 
         );

	int IsParallelTo( 
      // returns  1: this and other vectors are and parallel
      //         -1: this and other vectors are anti-parallel
      //          0: this and other vectors are not parallel
      //             or at least one of the vectors is zero
      const EM_3dVectorBase& v,
      double angle_tolerance=EM_DEFAULT_ANGLE_TOLERANCE// (default=EN_DEFAULT_ANGLE_TOLERANCE) radians
	  ) const; /*{
		  return EM_3dVectorBase::IsParallelTo(v,angle_tolerance);
	}*/
  bool IsPerpendicularTo(
        // returns true:  this and other vectors are perpendicular
        //         false: this and other vectors are not perpendicular
        //                or at least one of the vectors is zero
        const EM_3dVectorBase& other,                           // other vector     
        double angle_tolerance = EM_DEFAULT_ANGLE_TOLERANCE // optional angle tolerance (radians)
        ) const;
  // Returns:
  //   true if vector is the zero vector.
  bool IsZero() const;
  void Zero(); // set all coordinates to zero;
  bool Unitize();  // returns false if vector has zero length


	inline void Conv(T *a) const{
		a[0]=x[0];
		a[1]=x[1];
		a[2]=x[2];
	}
};

typedef EM_3dVectorBase<double>  EM_3dVector;
typedef EM_3dVectorBase<COMPLEX>  EM_3dVectorComplex;

template <typename T>
inline EM_3dVectorBase<T>::EM_3dVectorBase(){ Zero();}

template <typename T>
inline EM_3dVectorBase<T>::EM_3dVectorBase( const float* v )
{
	assert(v);
//  if (v) {
    x[0] = v[0]; x[1] = v[1]; x[2] = v[2];
//  }
//  else {
//    x[0] = x[1] = x[2] = 0.0;
//  }
}

template <typename T>
inline EM_3dVectorBase<T>::EM_3dVectorBase( const T* v )
{
	assert(v);
 // if (v) {
    x[0] = v[0]; x[1] = v[1]; x[2] = v[2];
//  }
//  else {
//    x[0] = x[1] = x[2] = 0.0;
//  }
}


template <typename T>
inline EM_3dVectorBase<T>::EM_3dVectorBase(T xx, T yy, T zz)
{x[0]=xx;x[1]=yy;x[2]=zz;}
/*
template <typename T>
inline EM_3dVectorBase<T>::EM_3dVectorBase<T>(const EM_3dVectorBase<T>& v)
{x[0]=v.x[0];x[1]=v.x[1];x[2]=v.x[2];}
*/
template <>
inline EM_3dVectorBase<EM_REAL>::EM_3dVectorBase(const EM_3dVectorBase<EM_REAL>& v)
{x[0]=v.x[0];x[1]=v.x[1];x[2]=v.x[2];}
template <>
inline EM_3dVectorBase<COMPLEX>::EM_3dVectorBase(const EM_3dVectorBase<EM_REAL>& v)
{x[0]=v.x[0];x[1]=v.x[1];x[2]=v.x[2];}
template <>
inline EM_3dVectorBase<COMPLEX>::EM_3dVectorBase(const EM_3dVectorBase<COMPLEX>& v)
{x[0]=v.x[0];x[1]=v.x[1];x[2]=v.x[2];}
template <>
inline EM_3dVectorBase<EM_REAL>::EM_3dVectorBase(const EM_3dVectorBase<COMPLEX>& v)
{	assert(false); }


template <typename T>
inline EM_3dVectorBase<T>::operator T*()
{
  return x;
}

template <typename T>
inline EM_3dVectorBase<T>::operator const T*() const
{
  return x;
}

template <typename T>
inline T EM_3dVectorBase<T>::operator*( const EM_3dVectorBase<T>& v ) const
{
  return (x[0]*v.x[0] + x[1]*v.x[1] + x[2]*v.x[2]);
}
/*
template <typename T>
EM_3dVectorBase<T>& EM_3dVectorBase<T>::operator*=(double d)
{
  x *= d;
  y *= d;
  x[2] *= d;
  return *this;
}
*/
template <typename T>
inline EM_3dVectorBase<T>& EM_3dVectorBase<T>::operator*=(T d)
{
  x[0] *= d;
  x[1] *= d;
  x[2] *= d;
  return *this;
}

template <typename T>
inline EM_3dVectorBase<T>& EM_3dVectorBase<T>::operator/=(double d)
{
  x[0] /= d;
  x[1] /= d;
  x[2] /= d;
  return *this;
}

template <typename T>
inline EM_3dVectorBase<T>& EM_3dVectorBase<T>::operator+=(const EM_3dVectorBase<T>& v)
{
  x[0] += v.x[0];
  x[1] += v.x[1];
  x[2] += v.x[2];
  return *this;
}

template <typename T>
inline EM_3dVectorBase<T>& EM_3dVectorBase<T>::operator-=(const EM_3dVectorBase<T>& v)
{
  x[0] -= v.x[0];
  x[1] -= v.x[1];
  x[2] -= v.x[2];
  return *this;
}

template <typename T>
inline EM_3dVectorBase<T> EM_3dVectorBase<T>::operator+(const EM_3dVectorBase<T>& v) const
{
  return EM_3dVectorBase<T>(x[0]+v.x[0], x[1]+v.x[1], x[2]+v.x[2]);
}

template <typename T>
inline EM_3dVectorBase<T> EM_3dVectorBase<T>::operator-(const EM_3dVectorBase<T>& v) const
{
  return EM_3dVectorBase<T>(x[0]-v.x[0], x[1]-v.x[1], x[2]-v.x[2]);
}

template <typename T>
inline EM_3dVectorBase<T> EM_3dVectorBase<T>::operator*(const  T &d ) const
{
  return EM_3dVectorBase<T>(x[0]*d,x[1]*d,x[2]*d);
}

template <typename T>
inline EM_3dVectorBase<T> EM_3dVectorBase<T>::operator/( double d ) const
{
  return EM_3dVectorBase<T>(x[0]/d,x[1]/d,x[2]/d);
}

template <typename T>
inline EM_3dVectorBase<T> EM_3dVectorBase<T>::operator-() const
{
  return EM_3dVectorBase<T>(-x[0],-x[1],-x[2]);
}

template <typename T>
inline const T & EM_3dVectorBase<T>::operator[](int i) const
{
  return x[i];
}

template <typename T>
inline T & EM_3dVectorBase<T>::operator[](int i)
{
  return x[i];
}

template <typename T>
bool EM_3dVectorBase<T>::IsZero() const
{
  return (x[0]==0.0 && x[1]==0.0 && x[2]==0.0);
}

template <typename T>
void EM_3dVectorBase<T>::Zero()
{
  x[0] = x[1] = x[2] = 0.0;
}
/*
template <>
void EM_3dVectorBase<COMPLEX>::Zero()
{
  x = x[1] = x[2] = COMPLEX(0.0,0.);
}
*/
template <typename T>
bool EM_3dVectorBase<T>::IsValid() const
{
  return ( EM_IsValid(x[0]) && EM_IsValid(x[1]) && EM_IsValid(x[2]) ) ? true : false;
}

template <typename T>
bool EM_3dVectorBase<T>::PerpendicularTo( const EM_3dVectorBase<T>& v )
{
  //bool rc = false;
  int i, j, k; 
  double a, b;
  k = 2;
  if ( fabs(v.x[1]) > fabs(v.x[0]) ) {
    if ( fabs(v.x[2]) > fabs(v.x[1]) ) {
      // |v.x[2]| > |v.x[1]| > |v.x|
      i = 2;
      j = 1;
      k = 0;
      a = v.x[2];
      b = -v.x[1];
    }
    else if ( fabs(v.x[2]) >= fabs(v.x[0]) ){
      // |v.x[1]| >= |v.x[2]| >= |v.x[0]|
      i = 1;
      j = 2;
      k = 0;
      a = v.x[1];
      b = -v.x[2];
    }
    else {
      // |v.x[1]| > |v.x| > |v.x[2]|
      i = 1;
      j = 0;
      k = 2;
      a = v.x[1];
      b = -v.x[0];
    }
  }
  else if ( fabs(v.x[2]) > fabs(v.x[0]) ) {
    // |v.x[2]| > |v.x[0]| >= |v.x[1]|
    i = 2;
    j = 0;
    k = 1;
    a = v.x[2];
    b = -v.x[0];
  }
  else if ( fabs(v.x[2]) > fabs(v.x[1]) ) {
    // |v.x| >= |v.x[2]| > |v.x[1]|
    i = 0;
    j = 2;
    k = 1;
    a = v.x[0];
    b = -v.x[2];
  }
  else {
    // |v.x| >= |v.x[1]| >= |v.x[2]|
    i = 0;
    j = 1;
    k = 2;
    a = v.x[0];
    b = -v.x[1];
  }
  double* this_v = &x[0];
  this_v[i] = b;
  this_v[j] = a;
  this_v[k] = 0.0;
  return (a != 0.0) ? true : false;
}

template <typename T>
int EM_3dVectorBase<T>::IsParallelTo( 
      // returns  1: this and other vectors are and parallel
      //         -1: this and other vectors are anti-parallel
      //          0: this and other vectors are not parallel
      //             or at least one of the vectors is zero
      const EM_3dVectorBase<T>& v,
      double angle_tolerance // (default=EM_DEFAULT_ANGLE_TOLERANCE) radians
      ) const
{
  int rc = 0;
  const double ll = Length()*v.Length();
  if ( ll > 0.0 ) {
    const double cos_angle = (x[0]*v.x[0] + x[1]*v.x[1] + x[2]*v.x[2])/ll;
    const double cos_tol = cos(angle_tolerance);
    if ( cos_angle >= cos_tol )
      rc = 1;
    else if ( cos_angle <= -cos_tol )
      rc = -1;
  }
  return rc;
}

template <typename T>
bool EM_3dVectorBase<T>::IsPerpendicularTo(
      // returns true:  this and other vectors are perpendicular
      //         false: this and other vectors are not perpendicular
      //                or at least one of the vectors is zero
      const EM_3dVectorBase<T>& v,
      double angle_tolerance // (default=EM_DEFAULT_ANGLE_TOLERANCE) radians
      ) const
{
  bool rc = false;
  const double ll = Length()*v.Length();
  if ( ll > 0.0 ) {
    if ( fabs((x[0]*v.x[0] + x[1]*v.x[1] + x[2]*v.x[2])/ll) <= sin(angle_tolerance) )
      rc = true;
  }
  return rc;
}

template <typename T>
double EM_3dVectorBase<T>::Length2() const{
	return x[0]*x[0]+x[1]*x[1]+x[2]*x[2];
}

template <typename T>
double EM_3dVectorBase<T>::Length() const
{
	return sqrt(x[0]*x[0]+x[1]*x[1]+x[2]*x[2]);
#if 0
  double len;
  double fx = fabs(x[0]);
  double fy = fabs(x[1]);
  double fz = fabs(x[2]);
  if ( fy >= fx && fy >= fz ) {
    len = fx; fx = fy; fy = len;
  }
  else if ( fz >= fx && fz >= fy ) {
    len = fx; fx = fz; fz = len;
  }

  // 15 September 2003 Dale Lear
  //     For small denormalized doubles (positive but smaller
  //     than DBL_MIN), some compilers/FPUs set 1.0/fx to +INF.
  //     Without the EM_DBL_MIN test we end up with
  //     microscopic vectors that have infinte length!
  //
  //     This code is absolutely necessary.  It is a critical
  //     part of the bug fix for RR 11217.
  if ( fx > EM_DBL_MIN ) 
  {
    len = 1.0/fx;
    fy *= len;
    fz *= len;
    len = fx*sqrt(1.0 + fy*fy + fz*fz);
  }
  else if ( fx > 0.0 && EM_IsFinite(fx) )
    len = fx;
  else
    len = 0.0;

  return len;
#endif
}

template <typename T>
bool EM_3dVectorBase<T>::Unitize()
{
  // 15 September 2003 Dale Lear
  //     Added the EM_DBL_MIN test.  See EM_3dVectorBase<T>::Length()
  //     for details.
  bool rc = false;
  double d = Length();
  if ( d > EM_DBL_MIN )
  {
    d = 1.0/d;
    x[0] *= d;
    x[1] *= d;
    x[2] *= d;
    rc = true;
  }
  else if ( d > 0.0 && EM_IsFinite(d) )
  {
    // This code is rarely used and can be slow.
    // It multiplies by 2^1023 in an attempt to 
    // normalize the coordinates.
    // If the renormalization works, then we're
    // ok.  If the renormalization fails, we
    // return false.
    EM_3dVectorBase<T> tmp;
    tmp.x[0] = x[0]*8.9884656743115795386465259539451e+307;
    tmp.x[1] = x[1]*8.9884656743115795386465259539451e+307;
    tmp.x[2] = x[2]*8.9884656743115795386465259539451e+307;
    d = tmp.Length();
    if ( d > EM_DBL_MIN )
    {
      d = 1.0/d;
      x[0] = tmp.x[0]*d;
      x[1] = tmp.x[1]*d;
      x[2] = tmp.x[2]*d;
      rc = true;
    }
    else
    {
      x[0] = 0.0;
      x[1] = 0.0;
      x[2] = 0.0;
    }
  }
  else
  {
    x[0] = 0.0;
    x[1] = 0.0;
    x[2] = 0.0;
  }

  return rc;
}



///////////////////////////////////////////////////////////////
//
// EM_3dVector utilities
//
///////////////////////////////////////////////////////////////
//
// common vectors
//


template <typename T>
EM_3dVectorBase<T> operator*(const T &d, const EM_3dVectorBase<T>& v)
{
  return EM_3dVectorBase<T>(d*v.x[0],d*v.x[1],d*v.x[2]);
}

EM_3dVectorBase<COMPLEX> operator*(const COMPLEX &d, const EM_3dVectorBase<EM_REAL>& v);
EM_3dVectorBase<COMPLEX> operator*( const EM_3dVectorBase<EM_REAL>& v, const COMPLEX &d);

template <typename T>
double EM_DotProduct( const EM_3dVectorBase<T>& a , const EM_3dVectorBase<T>& b )
{
  // inner (dot) product between 3d vectors
  return (a.x[0]*b.x[0] + a.x[1]*b.x[1] + a.x[2]*b.x[2]);
}


template <typename T>
EM_3dVectorBase<T> EM_CrossProduct( const EM_3dVectorBase<T>& a , const EM_3dVectorBase<T>& b )
{
  return EM_3dVectorBase<T>(a.x[1]*b.x[2] - b.x[1]*a.x[2], a.x[2]*b.x[0] - b.x[2]*a.x[0], a.x[0]*b.x[1] - b.x[0]*a.x[1] );
}
/*
template <typename T>
EM_3dVectorBase<T> EM_CrossProduct( const EM_3dVectorBase<T>& a , const EM_3dVectorBase<EM_REAL>& b )
{
  return EM_3dVectorBase<T>(a.x[1]*b.x[2] - b.x[1]*a.x[2], a.x[2]*b.x - b.x[2]*a.x, a.x*b.x[1] - b.x*a.x[1] );
}


EM_3dVectorBase<EM_REAL> EM_CrossProduct( const EM_3dVectorBase<EM_REAL>& a , const EM_3dVectorBase<EM_REAL>& b )
{
  return EM_3dVectorBase<EM_REAL>(a.x[1]*b.x[2] - b.x[1]*a.x[2], a.x[2]*b.x - b.x[2]*a.x, a.x*b.x[1] - b.x*a.x[1] );
}
*/
bool 
EM_IsOrthogonalFrame( // true if X, Y, Z are nonzero and mutually perpindicular
    const EM_3dVector&, // X
    const EM_3dVector&, // Y
    const EM_3dVector&  // Z 
    );
bool 
EM_IsOrthonormalFrame( // true if X, Y, Z are orthogonal and unit length
    const EM_3dVector&, // X
    const EM_3dVector&, // Y
    const EM_3dVector&  // Z 
    );
bool 
EM_IsRightHandFrame( // true if X, Y, Z are orthonormal and right handed
    const EM_3dVector&, // X
    const EM_3dVector&, // Y
    const EM_3dVector&  // Z 
    );

#define EM_UNSET_VECTOR EM_UNSET_VECTOR

