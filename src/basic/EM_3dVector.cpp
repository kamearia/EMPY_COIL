#include "basic/EM_Headers.h"
#include "basic/EM_Math.h"

EM_3dVectorBase<COMPLEX> operator*(const COMPLEX &d, const EM_3dVectorBase<EM_REAL>& v){
  return EM_3dVectorBase<COMPLEX>(d*v.x[0],d*v.x[1],d*v.x[2]);
}
EM_3dVectorBase<COMPLEX> operator*( const EM_3dVectorBase<EM_REAL>& v, const COMPLEX &d){
  return EM_3dVectorBase<COMPLEX>(d*v.x[0],d*v.x[1],d*v.x[2]);
}


#if 0
EM_3dVector::EM_3dVector()
{}
EM_3dVector::EM_3dVector( const float* v )
{
  if (v) {
    x = v[0]; y = v[1]; z = v[2];
  }
  else {
    x = y = z = 0.0;
  }
}
EM_3dVector::EM_3dVector( const double* v )
{
  if (v) {
    x = v[0]; y = v[1]; z = v[2];
  }
  else {
    x = y = z = 0.0;
  }
}

EM_3dVector::EM_3dVector(double xx,double yy,double zz)
{x=xx;y=yy;z=zz;}

EM_3dVector::EM_3dVector(const EM_3dVector& v)
{x=v.x;y=v.y;z=v.z;}

EM_3dVector::operator double*()
{
  return &x;
}

EM_3dVector::operator const double*() const
{
  return &x;
}
double EM_3dVector::operator*( const EM_3dVector& v ) const
{
  return (x*v.x + y*v.y + z*v.z);
}

EM_3dVector& EM_3dVector::operator*=(double d)
{
  x *= d;
  y *= d;
  z *= d;
  return *this;
}
EM_3dVector& EM_3dVector::operator/=(double d)
{
  x /= d;
  y /= d;
  z /= d;
  return *this;
}
EM_3dVector& EM_3dVector::operator+=(const EM_3dVector& v)
{
  x += v.x;
  y += v.y;
  z += v.z;
  return *this;
}

EM_3dVector EM_3dVector::operator+(const EM_3dVector& v) const
{
  return EM_3dVector(x+v.x, y+v.y, z+v.z);
}

EM_3dVector EM_3dVector::operator-(const EM_3dVector& v) const
{
  return EM_3dVector(x-v.x, y-v.y, z-v.z);
}

EM_3dVector EM_3dVector::operator*( double d ) const
{
  return EM_3dVector(x*d,y*d,z*d);
}
EM_3dVector EM_3dVector::operator/( double d ) const
{
  return EM_3dVector(x/d,y/d,z/d);
}
EM_3dVector EM_3dVector::operator-() const
{
  return EM_3dVector(-x,-y,-z);
}
double EM_3dVector::operator[](int i) const
{
  return ( (i<=0)?x:((i>=2)?z:y) );
}

double& EM_3dVector::operator[](int i)
{
  double* pd = (i<=0)? &x : ( (i>=2) ?  &z : &y);
  return *pd;
}

bool EM_3dVector::IsZero() const
{
  return (x==0.0 && y==0.0 && z==0.0);
}
void EM_3dVector::Zero()
{
  x = y = z = 0.0;
}
bool EM_3dVector::IsValid() const
{
  return ( EM_IsValid(x) && EM_IsValid(y) && EM_IsValid(z) ) ? true : false;
}

bool EM_3dVector::PerpendicularTo( const EM_3dVector& v )
{
  //bool rc = false;
  int i, j, k; 
  double a, b;
  k = 2;
  if ( fabs(v.y) > fabs(v.x) ) {
    if ( fabs(v.z) > fabs(v.y) ) {
      // |v.z| > |v.y| > |v.x|
      i = 2;
      j = 1;
      k = 0;
      a = v.z;
      b = -v.y;
    }
    else if ( fabs(v.z) >= fabs(v.x) ){
      // |v.y| >= |v.z| >= |v.x|
      i = 1;
      j = 2;
      k = 0;
      a = v.y;
      b = -v.z;
    }
    else {
      // |v.y| > |v.x| > |v.z|
      i = 1;
      j = 0;
      k = 2;
      a = v.y;
      b = -v.x;
    }
  }
  else if ( fabs(v.z) > fabs(v.x) ) {
    // |v.z| > |v.x| >= |v.y|
    i = 2;
    j = 0;
    k = 1;
    a = v.z;
    b = -v.x;
  }
  else if ( fabs(v.z) > fabs(v.y) ) {
    // |v.x| >= |v.z| > |v.y|
    i = 0;
    j = 2;
    k = 1;
    a = v.x;
    b = -v.z;
  }
  else {
    // |v.x| >= |v.y| >= |v.z|
    i = 0;
    j = 1;
    k = 2;
    a = v.x;
    b = -v.y;
  }
  double* this_v = &x;
  this_v[i] = b;
  this_v[j] = a;
  this_v[k] = 0.0;
  return (a != 0.0) ? true : false;
}

int EM_3dVector::IsParallelTo( 
      // returns  1: this and other vectors are and parallel
      //         -1: this and other vectors are anti-parallel
      //          0: this and other vectors are not parallel
      //             or at least one of the vectors is zero
      const EM_3dVector& v,
      double angle_tolerance // (default=EM_DEFAULT_ANGLE_TOLERANCE) radians
      ) const
{
  int rc = 0;
  const double ll = Length()*v.Length();
  if ( ll > 0.0 ) {
    const double cos_angle = (x*v.x + y*v.y + z*v.z)/ll;
    const double cos_tol = cos(angle_tolerance);
    if ( cos_angle >= cos_tol )
      rc = 1;
    else if ( cos_angle <= -cos_tol )
      rc = -1;
  }
  return rc;
}

bool EM_3dVector::IsPerpendicularTo(
      // returns true:  this and other vectors are perpendicular
      //         false: this and other vectors are not perpendicular
      //                or at least one of the vectors is zero
      const EM_3dVector& v,
      double angle_tolerance // (default=EM_DEFAULT_ANGLE_TOLERANCE) radians
      ) const
{
  bool rc = false;
  const double ll = Length()*v.Length();
  if ( ll > 0.0 ) {
    if ( fabs((x*v.x + y*v.y + z*v.z)/ll) <= sin(angle_tolerance) )
      rc = true;
  }
  return rc;
}
double EM_3dVector::Length() const
{
  double len;
  double fx = fabs(x);
  double fy = fabs(y);
  double fz = fabs(z);
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
}

bool EM_3dVector::Unitize()
{
  // 15 September 2003 Dale Lear
  //     Added the EM_DBL_MIN test.  See EM_3dVector::Length()
  //     for details.
  bool rc = false;
  double d = Length();
  if ( d > EM_DBL_MIN )
  {
    d = 1.0/d;
    x *= d;
    y *= d;
    z *= d;
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
    EM_3dVector tmp;
    tmp.x = x*8.9884656743115795386465259539451e+307;
    tmp.y = y*8.9884656743115795386465259539451e+307;
    tmp.z = z*8.9884656743115795386465259539451e+307;
    d = tmp.Length();
    if ( d > EM_DBL_MIN )
    {
      d = 1.0/d;
      x = tmp.x*d;
      y = tmp.y*d;
      z = tmp.z*d;
      rc = true;
    }
    else
    {
      x = 0.0;
      y = 0.0;
      z = 0.0;
    }
  }
  else
  {
    x = 0.0;
    y = 0.0;
    z = 0.0;
  }

  return rc;
}

bool
EM_3dVector::PerpendicularTo( 
      const EM_3dPoint& P0, const EM_3dPoint& P1, const EM_3dPoint& P2
      )
{
  // Find a the unit normal to a triangle defined by 3 points
  EM_3dVector V0, V1, V2, N0, N1, N2;

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
EM_3dVector operator*(double d, const EM_3dVector& v)
{
  return EM_3dVector(d*v.x,d*v.y,d*v.z);
}
double EM_DotProduct( const EM_3dVector& a , const EM_3dVector& b )
{
  // inner (dot) product between 3d vectors
  return (a.x*b.x + a.y*b.y + a.z*b.z);
}
EM_3dVector EM_CrossProduct( const EM_3dVector& a , const EM_3dVector& b )
{
  return EM_3dVector(a.y*b.z - b.y*a.z, a.z*b.x - b.z*a.x, a.x*b.y - b.x*a.y );
}

#endif
/*
EM_3dVector<COMPLEX> & operator *( const EM_3dVector & v, const COMPLEX &f) {
	return 
}
*/
bool EM_IsOrthogonalFrame( const EM_3dVector& X,  const EM_3dVector& Y,  const EM_3dVector& Z )
{
  // returns true if X, Y, Z is an orthogonal frame
  if (! X.IsValid() || !Y.IsValid() || !Z.IsValid() )
    return false;

  double lx = X.Length();
  double ly = Y.Length();
  double lz = Z.Length();
  if ( lx <=  EM_SQRT_EPSILON )
    return false;
  if ( ly <=  EM_SQRT_EPSILON )
    return false;
  if ( lz <=  EM_SQRT_EPSILON )
    return false;
  lx = 1.0/lx;
  ly = 1.0/ly;
  lz = 1.0/lz;
  double xy = /*(X.x*Y.x + X.y*Y.y + X.z*Y.z)*/ X*Y*lx*ly;
  double yz = /*(Y.x*Z.x + Y.y*Z.y + Y.z*Z.z)*/ Y*Z*ly*lz;
  double zx = /*(Z.x*X.x + Z.y*X.y + Z.z*X.z)*/ Z*X*lz*lx;
  if (    fabs(xy) > EM_SQRT_EPSILON 
       || fabs(yz) > EM_SQRT_EPSILON
       || fabs(zx) > EM_SQRT_EPSILON
     )
  {
    double t = 0.0000152587890625;
    if ( fabs(xy) >= t || fabs(yz)  >= t || fabs(zx) >= t )
      return false;

    // do a more careful (and time consuming check)
    // This fixes RR 22219 and 22276
    EM_3dVector V;
    V = (lx*ly)*EM_CrossProduct(X,Y);
    t = fabs(/*(V.x*Z.x + V.y*Z.y + V.z*Z.z)*/ V*Z*lz);
    if ( fabs(t-1.0) > EM_SQRT_EPSILON )
      return false;

    V = (ly*lz)*EM_CrossProduct(Y,Z);
    t = fabs(/*(V.x*X.x + V.y*X.y + V.z*X.z)*/ V*X*lx);
    if ( fabs(t-1.0) > EM_SQRT_EPSILON )
      return false;

    V = (lz*lx)*EM_CrossProduct(Z,X);
    t = fabs(/*(V.x*Y.x + V.y*Y.y + V.z*Y.z)*/ V*Y*ly);
    if ( fabs(t-1.0) > EM_SQRT_EPSILON )
      return false;
  }
  return true;
}


bool EM_IsOrthonormalFrame( const EM_3dVector& X,  const EM_3dVector& Y,  const EM_3dVector& Z )
{
  // returns true if X, Y, Z is an orthonormal frame
  if ( !EM_IsOrthogonalFrame( X, Y, Z ) )
    return false;
  double x = X.Length();
  if ( fabs(x-1.0) >  EM_SQRT_EPSILON )
    return false;
  x = Y.Length();
  if ( fabs(x-1.0) >  EM_SQRT_EPSILON )
    return false;
  x = Z.Length();
  if ( fabs(x-1.0) >  EM_SQRT_EPSILON )
    return false;

  return true;
}

bool EM_IsRightHandFrame( const EM_3dVector& X,  const EM_3dVector& Y,  const EM_3dVector& Z )
{
  // returns true if X, Y, Z is an orthonormal right hand frame
  if ( !EM_IsOrthonormalFrame(X,Y,Z) )
    return false;
  double x = EM_DotProduct( EM_CrossProduct( X, Y ), Z );
  if ( x <=  EM_SQRT_EPSILON )
    return false;
  return true;
}