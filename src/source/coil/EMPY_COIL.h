#pragma once
#include "EM_COIL_Element_base.h"
#include <fem.hpp>
#include <python_ngstd.hpp>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <cmath>

namespace ngfem
{
	class EMPY_LOOP {
	public:
		EM_COIL_Loop_base* loop;
		EMPY_LOOP(EM_REAL current, EM_REAL radius, EM_REAL z, EM_REAL radialWidth, EM_REAL axialWidth) {
			loop = new EM_COIL_Loop_base(current, radius, z, radialWidth, axialWidth);
		}
		~EMPY_LOOP() {
			delete loop;
		}
		
		class B_Field : public CoefficientFunction
		{
		public:
			EM_COIL_Loop_base* loop;
			B_Field(EM_COIL_Loop_base* loop)
				: CoefficientFunction(/*dimension = */ 3) {
				this->loop = loop;
			}

			virtual double Evaluate(const BaseMappedIntegrationPoint& mip) const override
			{
				return 0;
			}


			virtual void Evaluate(const BaseMappedIntegrationPoint& mip,
				FlatVector<> result) const override
			{
				auto pnt = mip.GetPoint();
				EM_3dPoint p = EM_3dPoint(pnt[0], pnt[1], pnt[2]);
				EM_3dVector b;
				b = loop->B_Field(p);
				result(0) = b[0];
				result(1) = b[1];
				result(2) = b[2];
			}

		};

	};
}