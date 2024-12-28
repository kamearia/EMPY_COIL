#pragma once
#include "basic\EM_3dPoint.h"
//typedef double EM_REAL
#include "fem/EM_Integrator.h"
class EM_COIL_Element_base {
public:
	virtual const char* Name() const { return "ARC"; }
	virtual EM_3dVector A_Field(const EM_3dPoint& x) const = 0;
	virtual EM_3dVector B_Field(const EM_3dPoint& x) const = 0;
};

class EM_COIL_Arc_base : public EM_COIL_Element_base {
public:
	EM_COIL_Arc_base(
		double current,
		const EM_3dPoint& center,
		double radius,
		double width,
		double thick,
		double start,
		double end
	);
	virtual const char* Name() const { return "ARC"; }
	virtual EM_3dVector A_Field(const EM_3dPoint& x) const;
	virtual EM_3dVector B_Field(const EM_3dPoint& x) const;
private:
	void CalcFields(const EM_REAL* X, EM_REAL* A, EM_REAL* B, struct field_param* param) const;
	void EM_COIL_Arc_base::Integrand::Initialize();
	EM_REAL  current;
	EM_3dPoint center;
	EM_REAL  radius;
	EM_REAL  axial_width, radial_width;
	EM_REAL  Euler_angle_alpha, Euler_angle_beta;
	EM_REAL  start_angle, end_angle;

	class Integrand :public EM_Integrand {
	public:
		Integrand(const EM_COIL_Arc* parent, int numOutput)
			:EM_Integrand(3, numOutput), parent(parent) {
			Initialize();
		}
		void Initialize();
		void Next();
		EM_REAL Factor() const { return dphi; }

	protected:
		const EM_COIL_Arc* parent;
		EM_REAL phi0, dphi, sina, cosa, sinb, cosb;
	};
};
class EM_COIL_Loop_base : public EM_COIL_Arc_base {
public:
	EM_COIL_Loop_base(EM_REAL current, EM_REAL radius, EM_REAL z,
		EM_REAL radialWidth, EM_REAL axialWidth);
	const char* Name() const { return "LOOP"; }
	virtual EM_3dVector A_Field(const EM_3dPoint& x) const;
	virtual EM_3dVector B_Field(const EM_3dPoint& x) const;
};