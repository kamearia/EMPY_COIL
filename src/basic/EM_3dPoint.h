////////////////////////////////////////////////////////////////
//
//   defines double precision pointin 3-dimension
//
////////////////////////////////////////////////////////////////
#pragma once

//#include "basic/EM_defines.h"
#include <math.h>
#include "basic/EM_3dVector.h"
#include "basic/EM_Math.h"

bool EM_IsValid(double x);
////////////////////////////////////////////////////////////////
//
//   EM_3dPoint
//
class EM_CLASS EM_3dPoint//: public ON_3dPoint
{
public:
	typedef double value_type; 
  double x, y, z;
    EM_3dPoint(){} //:ON_3dPoint(){}
	inline EM_3dPoint(double x,double y,double z);//:ON_3dPoint(x,y,z){}
	inline EM_3dPoint( const float* p );//:ON_3dPoint(p){}
	inline EM_3dPoint( const double* p );//:ON_3dPoint(p){}
	EM_3dPoint(const EM_3dPoint &p);//:ON_3dPoint(p){}
	EM_3dPoint(const EM_3dVector& );    // from 3d vector

  // (double*) conversion operators
  inline operator double*();
  inline operator const double*() const;

  inline EM_3dPoint& operator=(const double*); // point = double[3] support
  inline EM_3dPoint& operator=(const EM_3dVector&);
  inline EM_3dPoint& operator*=(double);
  inline EM_3dPoint& operator/=(double);
  inline EM_3dPoint& operator+=(const EM_3dPoint&);
  inline EM_3dPoint& operator-=(const EM_3dVector&);
  inline EM_3dPoint  operator*(double) const;
  inline EM_3dPoint  operator/(double) const;
  inline EM_3dPoint  operator+(const EM_3dPoint&) const;
  inline EM_3dPoint  operator+(const EM_3dVector&) const;
  inline EM_3dVector operator-(const EM_3dPoint&) const;
  inline bool operator==(const EM_3dPoint&) const;
  // index operators mimic double[3] behavior
  inline double& operator[](int);
  inline double operator[](int) const;

  // set 3d point value
  inline void Set(double x,double y,double z);
  inline double DistanceTo( const EM_3dPoint& ) const;
  inline double MaximumCoordinate() const; // absolute value of maximum coordinate
  inline void Zero(); // set all coordinates to zero;
  inline EM_3dPoint ToCylindrical() const;  // (r, theta, z)
  bool Coincide(const EM_3dPoint &p) const;
  bool CoincideX(const EM_3dPoint &p) const;
  /*
  Description
    If any coordinate of a point is EM_UNSET_VALUE,
    then the point is not valid.
  Returns:
    true if the point is valid.
  */
  inline bool IsValid() const;


	int CompareZXY(const EM_3dPoint *p, double tolerance=0.) const{
		if(z < p->z-tolerance) return -1;
		if(z > p->z+tolerance) return 1;
		if(y < p->y-tolerance) return -1;
		if(y > p->y+tolerance) return 1;
		if(x < p->x-tolerance) return -1;
		if(x > p->x+tolerance) return 1;
		return 0;
	}

	void Write(class EM_BinaryFile &ffile) const;
	void Read(EM_BinaryFile &ffile);

	inline void Conv(EM_REAL *a) const{
		a[0]=x;
		a[1]=y;
		a[2]=z;
	}

	double distance( const EM_3dPoint &node)
	{
		double dx = x - node.x;
		double dy = y - node.y;
		double dz = z - node.z;
	return EM_MAX(fabs(dx),EM_MAX(fabs(dy),fabs(dz)));

	}

	bool operator <(const EM_3dPoint &left) const {
		if(x < left.x) return true;
		if(x==left.x){
			if(y< left.y) return true;
			if(y==left.y){
				if(z<left.z) return true;
			}
		}
		return false;
	}
#if 0
	void DebugWrite(FILE *f) const; /*{
		fprintf(ffile_check->file, " %15.5e , %15.5e ,%15.5e ",
			x, y, z);
	}*/
#endif
	static bool Test();
};

inline EM_3dPoint::EM_3dPoint(double xx,double yy,double zz) // :x(xx),y(yy),z(zz)
	{x=xx;y=yy;z=zz;}

inline EM_3dPoint::EM_3dPoint( const float* p )
{
  if (p) {
    x = p[0]; y = p[1]; z = p[2];
  }
  else {
    x = y = z = 0.0;
  }
}

inline EM_3dPoint::EM_3dPoint( const double* p )
{
  if (p) {
    x = p[0]; y = p[1]; z = p[2];
  }
  else {
    x = y = z = 0.0;
  }
}

inline EM_3dPoint::EM_3dPoint(const EM_3dPoint& p) // : x(p.x),y(p.y),z(0.0)
{x=p.x;y=p.y;z=p.z;}

inline EM_3dPoint::EM_3dPoint(const EM_3dVector& v) // : x(p.x),y(p.y),z(0.0)
{x=v.x[0];y=v.x[1];z=v.x[2];}

inline EM_3dPoint::operator double*()
{
  return &x;
}

inline EM_3dPoint::operator const double*() const
{
  return &x;
}

inline EM_3dPoint& EM_3dPoint::operator*=(double d)
{
  x *= d;
  y *= d;
  z *= d;
  return *this;
}
inline EM_3dPoint& EM_3dPoint::operator/=(double d)
{
  const double one_over_d = 1.0/d;
  x *= one_over_d;
  y *= one_over_d;
  z *= one_over_d;
  return *this;
}


inline EM_3dPoint& EM_3dPoint::operator+=(const EM_3dPoint& p)
{
  x += p.x;
  y += p.y;
  z += p.z;
  return *this;
}
inline EM_3dPoint& EM_3dPoint::operator-=(const EM_3dVector& v)
{
  x -= v.x[0];
  y -= v.x[1];
  z -= v.x[2];
  return *this;
}
inline EM_3dPoint& EM_3dPoint::operator=(const double* p)
{
  if ( p ) {
    x = p[0];
    y = p[1];
    z = p[2];
  }
  else {
    x = y = z = 0.0;
  }
  return *this;
}
inline EM_3dPoint& EM_3dPoint::operator=(const EM_3dVector& v)
{
  x = v.x[0];
  y = v.x[1];
  z = v.x[2];
  return *this;
}
inline EM_3dPoint EM_3dPoint::operator*( double d ) const
{
  return EM_3dPoint(x*d,y*d,z*d);
}
inline EM_3dPoint EM_3dPoint::operator/( double d ) const
{
  const double one_over_d = 1.0/d;
  return EM_3dPoint(x*one_over_d,y*one_over_d,z*one_over_d);
}
inline EM_3dPoint EM_3dPoint::operator+( const EM_3dPoint& p ) const
{
  return EM_3dPoint(x+p.x,y+p.y,z+p.z);
}
inline EM_3dPoint EM_3dPoint::operator+( const EM_3dVector& v ) const
{
  return EM_3dPoint(x+v.x[0],y+v.x[1],z+v.x[2]);
}

inline bool EM_3dPoint::operator==( const EM_3dPoint& p ) const
{
  return (x==p.x&&y==p.y&&z==p.z)?true:false;
}
inline double EM_3dPoint::operator[](int i) const
{
  return ( (i<=0)?x:((i>=2)?z:y) );
}

inline double& EM_3dPoint::operator[](int i)
{
  double* pd = (i<=0)? &x : ( (i>=2) ?  &z : &y);
  return *pd;
}
inline bool EM_3dPoint::IsValid() const
{
  return (EM_IsValid(x) && EM_IsValid(y) && EM_IsValid(z) ) ? true : false;
}
inline void EM_3dPoint::Set(double xx, double yy, double zz)
{
  x = xx; y = yy; z = zz;
}
inline double EM_3dPoint::MaximumCoordinate() const
{
  double c = fabs(x); if (fabs(y)>c) c=fabs(y); if (fabs(z)>c) c=fabs(z);
  return c;
}

inline void EM_3dPoint::Zero()
{
  x = y = z = 0.0;
}
inline double EM_3dPoint::DistanceTo( const EM_3dPoint& p ) const
{
  return (p - *this).Length();
}

inline EM_3dPoint EM_3dPoint::ToCylindrical() const{
	EM_3dVector cp;
	cp[0]=sqrt(x*x+y*y);
	cp[1]=atan2(y,x);
	cp[2]=z;
	return cp;
}

inline EM_3dPoint operator*(double d, const EM_3dPoint&p){
  return EM_3dPoint(d*p.x,d*p.y,d*p.z);
}

int ComparePosition(const EM_3dPoint *a, const EM_3dPoint *b);
int ComparePositionStrict(const EM_3dPoint *a, const EM_3dPoint *b);
int CompareX(const EM_3dPoint *a, const EM_3dPoint *b);
int CompareY(const EM_3dPoint *a, const EM_3dPoint *b);
int CompareZ(const EM_3dPoint *a, const EM_3dPoint *b);

extern const EM_3dPoint  EM_origin; // (0.0, 0.0, 0.0)
extern const EM_3dPoint  EM_UNSET_POINT; // (ON_UNSET_VALUE,ON_UNSET_VALUE,ON_UNSET_VALUE)

//#define EM_origin ON_origin
//#define EM_UNSET_POINT ON_UNSET_POINT
	
template<>
inline EM_3dVectorBase<double>::EM_3dVectorBase( const EM_3dPoint &p ){
	this->x[0]=p.x;
	this->x[1]=p.y;
	this->x[2]=p.z;
}
template<>
inline EM_3dVectorBase<COMPLEX>::EM_3dVectorBase( const EM_3dPoint &p ){
	assert(false);
}

EM_DECL
inline double EM_TriangleArea(const EM_3dPoint& p1, const EM_3dPoint& p2, const EM_3dPoint& p3){
	EM_3dVector a,b;
	a=p2-p1;
	b=p3-p1;
	return EM_CrossProduct(a,b).Length()/2.;
}

inline EM_3dVector EM_3dPoint::operator-( const EM_3dPoint& p ) const
{
  return EM_3dVector(x-p.x,y-p.y,z-p.z);
}

template <typename T>
bool
EM_3dVectorBase<T>::PerpendicularTo( 
      const EM_3dPoint& P0, const EM_3dPoint& P1, const EM_3dPoint& P2
      )
{
  // Find a the unit normal to a triangle defined by 3 points
  EM_3dVectorBase<T> V0, V1, V2, N0, N1, N2;

  Zero();

  V0 = P2 - P1;
  V1 = P0 - P2;
  V2 = P1 - P0;

  N0 = EM_CrossProduct( V1, V2 );
  if ( !N0.Unitize() )
    return false;
  N1 = EM_CrossProduct( V2, V0 );
  if ( !N1.Unitize() )
    return false;
  N2 = EM_CrossProduct( V0, V1 );
  if ( !N2.Unitize() )
    return false;

  const double s0 = 1.0/V0.Length();
  const double s1 = 1.0/V1.Length();
  const double s2 = 1.0/V2.Length();

  // choose normal with smallest total error
  const double e0 = s0*fabs(EM_DotProduct(N0,V0)) + s1*fabs(EM_DotProduct(N0,V1)) + s2*fabs(EM_DotProduct(N0,V2));
  const double e1 = s0*fabs(EM_DotProduct(N1,V0)) + s1*fabs(EM_DotProduct(N1,V1)) + s2*fabs(EM_DotProduct(N1,V2));
  const double e2 = s0*fabs(EM_DotProduct(N2,V0)) + s1*fabs(EM_DotProduct(N2,V1)) + s2*fabs(EM_DotProduct(N2,V2));

  if ( e0 <= e1 ) {
    if ( e0 <= e2 ) {
      *this = N0;
    }
    else {
      *this = N2;
    }
  }
  else if (e1 <= e2) {
    *this = N1;
  }
  else {
    *this = N2;
  }
  
  return true;
}

