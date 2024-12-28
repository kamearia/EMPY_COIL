#pragma once
#include "basic/EM_3dPoint.h"
#include <fem.hpp>
#include <fstream>

using namespace ngfem;
template <class T>
class B_CF : public CoefficientFunction
{
public:
	T* coil;
	B_CF(T* coil)
		: CoefficientFunction(3) {
		this->coil = coil;
	}

	virtual double Evaluate(const BaseMappedIntegrationPoint& mip) const override
	{
		return 0;
	}

	virtual void Evaluate(const BaseMappedIntegrationPoint& mip,
		FlatVector<> result) const override
	{
		auto pnt = mip.GetPoint();
		EM_3dVector v = coil->B_Field(EM_3dPoint(pnt[0], pnt[1], pnt[2]));
		result(0) = v[0];
		result(1) = v[1];
		result(2) = v[2];
	}
};

template <class T>
class A_CF : public CoefficientFunction
{
public:
	T* coil;
	A_CF(T* coil)
		: CoefficientFunction(3) {
		this->coil = coil;
	}

	virtual double Evaluate(const BaseMappedIntegrationPoint& mip) const override
	{
		return 0;
	}

	virtual void Evaluate(const BaseMappedIntegrationPoint& mip,
		FlatVector<> result) const override
	{
		auto pnt = mip.GetPoint();
		EM_3dVector v = coil->A_Field(EM_3dPoint(pnt[0], pnt[1], pnt[2]));
		result(0) = v[0];
		result(1) = v[1];
		result(2) = v[2];
}
};

template <class T>
class Omega_CF : public CoefficientFunction
{
public:
	T* coil;
	Omega_CF(T* coil)
		: CoefficientFunction(3) {
		this->coil = coil;
	}

	virtual double Evaluate(const BaseMappedIntegrationPoint& mip) const override
	{
		auto pnt = mip.GetPoint();
		return coil->Omega_Field(EM_3dPoint(pnt[0], pnt[1], pnt[2]));
	}

};


class EM_COIL_Element{
public:
	EM_COIL_Element() {
	}
	virtual EM_3dVector A_Field(const EM_3dPoint& x) const { return EM_3dPoint(0, 0, 0); }
	virtual EM_3dVector B_Field(const EM_3dPoint& x) const { return EM_3dPoint(0, 0, 0); }
	virtual double Omega_Field(const EM_3dPoint& x) const { return 0; }
	virtual const char* Name() const { return NULL; };

protected:
	static EM_REAL LINE_APP_LIMIT;
	struct field_param {
		EM_REAL AU, AUSQ, AL, ALSQ, ARU, ARL, ZZU, ZUSQ, ZZL, ZLSQ, ARZ[4], RHO,
			RHOSQ, AA, ASQ, TAR, Z1, ZSQ, ARZSQ, AREA;
		EM_REAL DFT;
		int calcA, calcB, symmetricArc;
	};
/*
	std::ofstream* file;
	void OpenFile() {
		file = new std::ofstream("debug_out.txt");
	}
*/
};

#include <list>
class  EM_COIL : public EM_COIL_Element {
public:
	EM_COIL(){}
	virtual ~EM_COIL() {};

	virtual EM_3dVector A_Field(const EM_3dPoint& x) const {
		EM_3dVector A = EM_3dVector(0, 0, 0);
		for (EM_COIL_Element* e : elements) {
			A += e->A_Field(x);
		}
		return A;
	}
	virtual EM_3dVector B_Field(const EM_3dPoint& x) const {
		EM_3dVector B = EM_3dVector(0, 0, 0);
		for (EM_COIL_Element* e : elements) {
			B +=e->B_Field(x);
		}
		return B;
	}

	virtual CoefficientFunction* B() {
		return new B_CF<EM_COIL>(this);
	}
	virtual CoefficientFunction* A() {
		return new A_CF<EM_COIL>(this);
	}

	std::string Add(EM_COIL_Element &e);

private:
	std::list< EM_COIL_Element* > elements;

};

class UNIFORM_FIELD : public EM_COIL_Element
{
public:
	UNIFORM_FIELD(double bx, double by, double bz, int dir);
	virtual const char *Name() const { return "UINF"; }
	EM_3dVector A_Field(const EM_3dPoint &x) const;
	EM_3dVector B_Field(const EM_3dPoint &x) const;
	double Omega_Field(const EM_3dPoint& x) const {
		double mu = 4.e-7 * EM_PI;
		return B_vector* x / mu;
	}

	virtual CoefficientFunction* B() {
		return new B_CF< UNIFORM_FIELD>(this);
	}
	virtual CoefficientFunction* A() {
		return new A_CF<UNIFORM_FIELD>(this);
	}
	virtual CoefficientFunction* Omega() {
		return new Omega_CF<UNIFORM_FIELD>(this);
	}
private:
	EM_3dVector B_vector;
	int    A_direction;
};

class INFINITE_LINE : public EM_COIL_Element {
public:
	INFINITE_LINE(double current, double px, double py, double pz
		, double dx, double dy, double dz);
	virtual const char *Name() const { return "LINE"; }
	EM_3dVector A_Field(const EM_3dPoint &x) const;
	EM_3dVector B_Field(const EM_3dPoint &x) const;
	virtual CoefficientFunction* B() {
		return new B_CF<INFINITE_LINE>(this);
	}
	virtual CoefficientFunction* A() {
		return new A_CF<INFINITE_LINE>(this);
	}

private:
	EM_REAL   current;
	EM_3dPoint  on_point;
	EM_3dVector direction;
	int    A_direction;
};


class EM_COIL_Block : public EM_COIL_Element {
public:
	EM_COIL_Block(){}
	EM_COIL_Block(
		double current,
		const EM_3dPoint& start,
		const EM_3dPoint& end,
		const EM_3dVector& dx,
		const EM_3dVector& dy
	);
	EM_COIL_Block(
		double current,
		double sx, double sy, double sz,
		double ex, double ey, double ez,
		double width,
		double thick
	) : EM_COIL_Block(current,
		EM_3dPoint(sx, sy, sz), EM_3dPoint(ex, ey, ez),
		width, thick) {}

	EM_COIL_Block(
		double current,
		const EM_3dPoint &start,
		const EM_3dPoint &end,
		double width,
		double thick
		);

	virtual const char *Name() const { return "GCE"; }
	void Initialize(); 
	virtual EM_3dVector A_Field(const EM_3dPoint &x) const;
	virtual EM_3dVector B_Field(const EM_3dPoint &x) const;

	virtual CoefficientFunction* B() {
		return new B_CF<EM_COIL_Block>(this);
	}
	virtual CoefficientFunction* A() {
		return new A_CF<EM_COIL_Block>(this);
	}

private:
	void normal_vectors();
	EM_REAL  current;
	EM_3dPoint start_point;
	EM_3dPoint end_point;
	EM_3dVector  half_width_vector_a;
	EM_3dVector   half_width_vector_b;
	EM_REAL length, area, center[3], Exyz[9], S1, S2;

	static EM_REAL inverse_r_integral_line(EM_REAL A[2], EM_REAL B[2], EM_REAL z);
	static EM_REAL inverse_r_integral_rectangle(EM_REAL s1, EM_REAL s2, EM_REAL x, EM_REAL y, EM_REAL z);

	void CalcFields(const EM_REAL *X, EM_REAL *A, EM_REAL *B, field_param *param) const;
	void CalcFields_Waseda(const EM_REAL *X, EM_REAL *A, EM_REAL *B, struct field_param *param) const;

};

class EM_COIL_Arc : public EM_COIL_Element {
public:
	EM_COIL_Arc(){}
	EM_COIL_Arc(
		double current,
		double radius,
		double z,
		double radialWidth,
		double axialWidth,
		double start_angle,
		double end_angle
	);
	EM_COIL_Arc(
		double current,
		double radius,
		const EM_3dPoint& center,
		double radialWidth,
		double axialWidth,
		double start_angle,
		double end_angle
	) :current(current), radius(radius), center(center)
		, axial_width(axialWidth), radial_width(radialWidth)
		, Euler_angle_alpha(0.), Euler_angle_beta(0.)
		, start_angle(start_angle), end_angle(end_angle)
	{}
	void SetCenter(double x, double y, double z);
	virtual const char *Name() const { return "ARC"; }
	virtual EM_3dVector A_Field(const EM_3dPoint& x) const;
	virtual EM_3dVector B_Field(const EM_3dPoint& x) const;
	virtual CoefficientFunction* B() {
		return new B_CF<EM_COIL_Arc>(this);
	}
	virtual CoefficientFunction* A() {
		return new A_CF<EM_COIL_Arc>(this);
	}

private:
	EM_REAL  current;
	EM_3dPoint center;
	EM_REAL  radius;
	EM_REAL  axial_width, radial_width;
	EM_REAL  Euler_angle_alpha, Euler_angle_beta;
	EM_REAL  start_angle, end_angle;

	struct field_sum{
		EM_REAL ARHO, ATHWAV, BZWAVE, BRHO, BTHETA;
	};
	static void segment_integral(EM_REAL PHI1, EM_REAL PHI2, field_param *param, field_sum *sum);
	static void arc_integration(EM_REAL P1, EM_REAL P2, field_param *param, field_sum *sum);
	static void arc_field_integrand(EM_REAL ALPHA, EM_REAL *F, void *param0);
	static void field_ARC_line(EM_REAL PHI1, EM_REAL PHI2, field_param *param, field_sum *sum);
	void CalcFields(const EM_REAL *X, EM_REAL *A, EM_REAL *B, struct field_param *param) const;

	friend class EM_COIL_Loop;
	};

class EM_COIL_Loop : public EM_COIL_Arc {
public:
	const char *Name() const { return "LOOP"; }
	EM_COIL_Loop(){}
	EM_COIL_Loop(EM_REAL current, EM_REAL radius, EM_REAL z, EM_REAL radialWidth, EM_REAL axialWidth);
	virtual EM_3dVector A_Field(const EM_3dPoint& x) const;
	virtual EM_3dVector B_Field(const EM_3dPoint& x) const;

	virtual CoefficientFunction* B() {
		B_CF< EM_COIL_Loop>* tmp=new B_CF< EM_COIL_Loop>(this);
		return tmp;
	}
	virtual CoefficientFunction* A() {
		return new A_CF< EM_COIL_Loop>(this);
	}

private:
	void field_LOOP(EM_REAL R, EM_REAL Z, EM_REAL CR, EM_REAL CZ, EM_REAL CA, EM_REAL CB, EM_REAL CI,
		EM_REAL *PSI, EM_REAL *BR, EM_REAL *BZ, field_param *param) const;
	static void bring0_(double R, double Z, double CR, double CZ, double CI,
		double& BR, double& BZ);
	static double psi_(double R, double Z, double CR, double CZ, double CI);
	static void CELIDD(double KX, double& ZK, double& ZE, int& ILL);

};
