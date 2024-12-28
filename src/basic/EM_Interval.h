#pragma once
#include "basic/EM_3dPoint.h"
class EM_File;
////////////////////////////////////////////////////////////////
//
//   ON_Interval
//
template <typename T>
class EM_Interval {
public:
  ////////
  // The default constructor creates the empty set interval (ON_UNSET_VALUE,ON_UNSET_VALUE)
	EM_Interval(){
		m_t[0] = m_t[1] = 0 /*EM_UNSET_VALUE*/; 
	}

	EM_Interval(T t0, T t1){
		Set(t0,t1);
	}

  // Interval = [m_t[0], m_t[1]]
  T m_t[2];

  /*
  Description:
    Sets interval to [t0,t1]
  Parameters:
    t0 - [in]
    t1 - [in]
  See Also:
    ON_Interval::ON_Interval( double, double )
  */
	inline void Set(T t0, T t1) {
		m_t[0] = t0;
		m_t[1] = t1;
	}

	inline void Reset() {
		m_t[0] = INT_MAX;
		m_t[1] = 0;
	}

	bool operator ==(const EM_Interval<T> & other) const {
		if(m_t[0] != other.m_t[0]
		  || m_t[1] != other.m_t[1]) return false;
		return true;
	}

	inline T Mid() const {// returns 0.5*(m_t[0] + m_t[1])
		return (m_t[0]+m_t[1])/2;
	} 
	inline T Length() const {
		return m_t[1]-m_t[0];
//		return ( EM_IsValid(m_t[0]) && EM_IsValid(m_t[1]) ) ? m_t[1]-m_t[0] : 0.0;
	}

	inline bool IsInside(const T &t) const {
		if(t<m_t[0] || t>m_t[1]) return false;
		return true;
	}

	inline void Add(const T &t){
		if(m_t[0]> t) m_t[0]=t;
		if(m_t[1]< t) m_t[1]=t;
	}

	inline void Add(const EM_Interval &other){
		if(other.m_t[0]< m_t[0]) m_t[0]=other.m_t[0];
		if(other.m_t[1]> m_t[1]) m_t[1]=other.m_t[1];
	}

	void Write(EM_File *ffile) const;
/*
	{
		FWRITE(&m_t[0], sizeof(T), 1, ffile);
		FWRITE(&m_t[1], sizeof(T), 1, ffile);
	}
*/
	void Read(EM_File *ffile);
/*
	{
		FREAD(&m_t[0], sizeof(T), 1, ffile);
		FREAD(&m_t[1], sizeof(T), 1, ffile);
	}
*/
	inline T Lower() const { return m_t[0]; }
	inline T Upper() const { return m_t[1]; }
	inline void Shift(T w) { 
		m_t[0] +=w;
		m_t[1] +=w;
	}
	T &operator ++(){ return ++m_t[1]; }
/*
	EM_Interval<T>& operator=( const EM_Interval<T>& src ){
		m_t[0]=src.m_t[0];
		m_t[1]=src.m_t[1];
		return *this;
	}
*/
};

class EM_3dInterval{
public:
	EM_3dInterval(){
		for(int i=0; i<3; ++i)
			m_Interval[i]=EM_Interval<double>(EM_DBL_MAX, EM_DBL_MIN);
	}
	inline EM_3dInterval &Add(const EM_3dPoint &point){
		for(int i=0; i<3; ++i)
			m_Interval[i].Add(point[i]);
		return *this;
	}

	inline EM_3dInterval &Add(const EM_3dInterval & interval){
		for(int i=0; i<3; ++i)
			m_Interval[i].Add(interval.m_Interval[i]);
		return *this;
	}

	inline const EM_3dPoint Upper() const {
		return EM_3dPoint(m_Interval[0].Upper(),
						  m_Interval[1].Upper(),	
						  m_Interval[2].Upper());
	}
	inline const EM_3dPoint Lower() const {
		return EM_3dPoint(m_Interval[0].Lower(),
						  m_Interval[1].Lower(),	
						  m_Interval[2].Lower());
	}
	inline bool IsInside(const EM_3dPoint &p) const{
		for(int i=0; i<3; ++i) 
			if(!m_Interval[i].IsInside(p[i]) ) return false;
		return true;
	}

private:
	EM_Interval<double> m_Interval[3];
};
