#include <iostream>
#include <string>
#include <fem.hpp>
#include <python_ngstd.hpp>
#include <pybind11/pybind11.h>
#include <pybind11/operators.h>
#include <pybind11/stl.h>
#include <cmath>
#include "source/coil/EM_COIL.h"
#include "basic/EM_3dVector.h"
#include <fstream>

namespace ngfem
{
#if 0
	template <class T> class  EMPY_Tmplate {

	public:
		T* coil;
//		std::ofstream *out;
		EMPY_Tmplate<T>() {}
		EMPY_Tmplate<T>(double , double , double , int);
		EMPY_Tmplate<T>(double, double, double, double, double);
		EMPY_Tmplate<T>(double, double, double, double, double, double);
		EMPY_Tmplate<T>(double, double, double, double, double, double, double);
		EMPY_Tmplate<T>(double, double, double, double, double, double, double, double, double);
//		void OpenFile() {
//			out = new std::ofstream("out.txt");
//		}
		
		class Omega_Field : public CoefficientFunction
		{
		public:
			T* coil;
	
			Omega_Field(T* coil)
				: CoefficientFunction(/*dimension = */ 1) {
				this->coil = coil;
			}

			// ********************* Necessary **********************************
			virtual double Evaluate(const BaseMappedIntegrationPoint& mip) const override
			{
				// maybe do something with mip and trafo?
				auto pnt = mip.GetPoint();
				return coil->Omega_Field(EM_3dPoint(pnt[0], pnt[1], pnt[2]));

				//double mu = 4.e-7 * M_PI;
				//return  B0*pnt[2] / mu;
			}
		};
		class B_Field : public CoefficientFunction
		{
		public:
			//T* coil;
			EM_COIL_Element* coil;
	//		std::ofstream* out;
			B_Field(EM_COIL_Element* coil/*, std::ofstream* out*/)
				: CoefficientFunction(/*dimension = */ 3) {
				this->coil = coil;
//				this->out = out;
//				coil->SetFile(out);
			}

			virtual double Evaluate(const BaseMappedIntegrationPoint& mip) const override
			{ return 0; }

			// ********************* Necessary **********************************
			virtual void Evaluate(const BaseMappedIntegrationPoint& mip,
				FlatVector<> result) const override
			{
				auto pnt = mip.GetPoint();
				EM_3dVector v = coil->B_Field(EM_3dPoint(pnt[0], pnt[1], pnt[2]));
				result(0) = v[0];
				result(1) = v[1];
				result(2) = v[2];
//				*out << "  x: " << pnt[0] << "  y= " << pnt[1] << "  z= " << pnt[2]
//					<< "  bx: " << v[0] << "  by= " << v[1] << "  bz= " << v[2] 
//					<<std::endl ;
//				//*out<<std::fflush:
			}

		};

		class A_Field : public CoefficientFunction
		{
		public:
			T* coil;
			A_Field(T* coil)
				: CoefficientFunction(/*dimension = */ 3) {
				this->coil = coil;
			}
			
			virtual double Evaluate(const BaseMappedIntegrationPoint& mip) const override { return 0; }

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
	public:
		CoefficientFunction* B() {
			return new B_Field(coil/*, out */);
		}
		CoefficientFunction* A() {
			return new A_Field(coil);
		}
		CoefficientFunction* Omega() {
			return new Omega_Field(coil);
		}


	};
	
	template<> 
	EMPY_Tmplate< UNIFORM_FIELD>::
	EMPY_Tmplate<UNIFORM_FIELD>(double bx, double by, double bz, int dir) {
		coil = new UNIFORM_FIELD(bx, by, bz, dir);
	//	OpenFile();
	}

	template<>
	EMPY_Tmplate<EM_COIL_Block>::
		EMPY_Tmplate<EM_COIL_Block>(double current,
			double sx, double sy, double sz,
			double ex, double ey, double ez,
			double width, double thick) {
		coil = new EM_COIL_Block(current, EM_3dPoint(sx,sy,sz), EM_3dPoint(ex,ey,ez),
			width, thick);
		//	OpenFile();
	}

	template<>
	EMPY_Tmplate< EM_COIL_Arc>::
		EMPY_Tmplate<EM_COIL_Arc>(double current, double radius, double zc, 
			double radialWidth, double axialWidth, double start_angle, double end_angle) {
		coil = new EM_COIL_Arc(current, radius, zc, radialWidth, axialWidth,
			start_angle* EM_DEGREES_TO_RADIANS, end_angle* EM_DEGREES_TO_RADIANS);
		//	OpenFile();
	}

	template<>
	EMPY_Tmplate< EM_COIL_Loop>::
		EMPY_Tmplate<EM_COIL_Loop>(double current, double radius, double z, double radialWidth, double axialWidth) {
		coil = new EM_COIL_Loop(current, radius, z, radialWidth, axialWidth);
	//	OpenFile();
	}
	
	template<>
	EMPY_Tmplate<INFINITE_LINE>::
		EMPY_Tmplate<INFINITE_LINE>(double current, double px, double py, double pz
			, double dx, double dy, double dz) {
		coil = new INFINITE_LINE(current, px, py, pz, dx, dy, dz);
	//		OpenFile();
	}
	
	typedef EMPY_Tmplate< UNIFORM_FIELD> UNIF;
	typedef EMPY_Tmplate< INFINITE_LINE> LINE;
	typedef EMPY_Tmplate< EM_COIL_Block> BLOCK;
	typedef EMPY_Tmplate< EM_COIL_Arc> ARC;
	typedef EMPY_Tmplate< EM_COIL_Loop> LOOP; 

	class LOOP {

	public:
		EM_COIL_Loop* coil;
		LOOP(double current, double radius, double z, double radialWidth, double axialWidth) {
			this->coil = new EM_COIL_Loop(current, radius, z, radialWidth, axialWidth);
		}

		class B_Field : public CoefficientFunction
		{
		public:
			EM_COIL_Loop* coil;
			B_Field(EM_COIL_Loop* coil)
				: CoefficientFunction(/*dimension = */ 3) {
				this->coil = coil;
			}

			virtual double Evaluate(const BaseMappedIntegrationPoint& mip) const override
			{
				return 0;
			}

			// ********************* Necessary **********************************
			virtual void Evaluate(const BaseMappedIntegrationPoint& mip,
				FlatVector<> result) const override
			{
				auto pnt = mip.GetPoint();
				EM_3dVector& v = coil->B_Field(EM_3dPoint(pnt[0], pnt[1], pnt[2]));
				result(0) = v[0];
				result(1) = v[1];
				result(2) = v[2];
			}

		};

		class A_Field : public CoefficientFunction
		{
		public:
			EM_COIL_Loop* coil;
			A_Field(EM_COIL_Loop* coil)
				: CoefficientFunction(/*dimension = */ 3) {
				this->coil = coil;
			}

			virtual double Evaluate(const BaseMappedIntegrationPoint& mip) const override
			{
				return 0;
			}

			// ********************* Necessary **********************************
			virtual void Evaluate(const BaseMappedIntegrationPoint& mip,
				FlatVector<> result) const override
			{
				auto pnt = mip.GetPoint();
				EM_3dVector& v = coil->A_Field(EM_3dPoint(pnt[0], pnt[1], pnt[2]));
				result(0) = v[0];
				result(1) = v[1];
				result(2) = v[2];
			}

		};
		CoefficientFunction* B() {
			return new B_Field(coil);
		}
		CoefficientFunction* A() {
			return new A_Field(coil);
		}
	};
#endif

#ifdef PYTHON
/*
#include <pybind11/stl.h>
#include <pybind11/numpy.h>
#include <pybind11/complex.h>
#include <memory>
*/


namespace py = pybind11;

PYBIND11_MODULE(EMPY_Field, m) {

		py::class_<EM_COIL_Element>(m, "EM_COIL_Element")
			.def(py::init<>())
		;

		py::class_<EM_COIL, EM_COIL_Element>(m, "EM_COIL")
			.def(py::init<>())
			.def("Add", &EM_COIL::Add)
			.def("B", &EM_COIL::B)
			.def("A", &EM_COIL::A)
			;

		py::class_<UNIFORM_FIELD, EM_COIL_Element>(m, "UNIF")
			.def(py::init<double, double, double, int>())
			.def("B", &UNIFORM_FIELD::B)
			.def("A", &UNIFORM_FIELD::A)
			.def("Omega", &UNIFORM_FIELD::Omega)
			;

		py::class_<INFINITE_LINE, EM_COIL_Element>(m, "LINE")
			.def(py::init<double, double, double, double, double, double, double>())
			.def("B", &INFINITE_LINE::B)
			.def("A", &INFINITE_LINE::A)
			;

		py::class_<EM_COIL_Block, EM_COIL_Element>(m, "BLOCK")
			.def(py::init<double, const EM_3dPoint&, const EM_3dPoint&,
						const EM_3dVector&, const EM_3dVector&>())
			.def(py::init<double, const EM_3dPoint&, const EM_3dPoint&, double, double>())
			.def(py::init<double, double, double, double, double, double, double, double, double>())
			.def("B", &EM_COIL_Block::B)
			.def("A", &EM_COIL_Block::A)
			;

		py::class_<EM_COIL_Arc, EM_COIL_Element>(m, "ARC")
			.def(py::init<double, double, const EM_3dPoint&, double, double, double, double>())
			.def(py::init<double, double, double, double, double, double, double>())
			.def("B", &EM_COIL_Arc::B)
			.def("A", &EM_COIL_Arc::A)
			;

		py::class_<EM_COIL_Loop, EM_COIL_Arc, EM_COIL_Element>(m, "LOOP")
			.def(py::init<double, double, double, double, double>())
			.def("B", &EM_COIL_Loop::B)
			.def("A", &EM_COIL_Loop::A)
			;


		py::class_<EM_3dPoint>(m, "EM_3dPoint")
			.def(py::init<double, double, double>())
			;

		py::class_<EM_3dVector>(m, "EM_3dVector")
			.def(py::init<double, double, double>())
//			.def(py::self + py::self)
			;
	}
/*

void ExportMyCoefficient(py::module m)
{
	py::class_<MyCoefficientFunction, shared_ptr<MyCoefficientFunction>, CoefficientFunction>
		(m, "MyCoefficient", "CoefficientFunction that returns x*y.")
		.def(py::init<>())
		;
}
*/
}
#endif