#ifdef COILBASE
#include "source/coil/EM_COIL.h"
#include "basic/EM_Math.h"
#include "basic/my_math.h"
#else
#include "basic/EM_Headers.h"
#include "domain/EM_DomainVectors.h"
//#include "Header/circuit.h"
#include "mesh/model.h"
#include "basic/fortran_function.h"
#include "io/files.h"
#include "basic/list_int.h"
#include "matrix/vector.h"
#include "source/regularization.h"
#include "post/post_files/output.h"
#include "EM_COIL_Assembly.h"
#include "source/coil/coil_field.h"
#include "source/coil/EM_MeshedCOIL.h"
#include "mesh/EM_MeshSet.h"
#include "boundary/EM_PotentialInterface.h"
#include "source/coil/EM_A_SourcePotentials.h"
#include "potential/EM_PotentialConnection.h"
#include "source/EM_Circuit.h"

extern void reverse_translated_position_by_motion(EM_REAL x[3], EM_REAL xt[3], const MOTION *motion);
//extern EM_REAL inductance_between_COIL_series( EM_COIL *sourceTo, EM_COIL *sourceFrom);
#endif

std::string EM_COIL::Add(EM_COIL_Element& e) {
	elements.push_back(&e);

	std::string str;
	int k = 0;
	for (EM_COIL_Element *e : elements) {
		k = k + 1;
		const char* name = e->Name();
//		return std::string(name);
		str +=std::string(name);
	}

//	std::cout << str << endl;
	return str;
}

EM_REAL EM_COIL_Element::LINE_APP_LIMIT = 1000;
#ifndef COILBASE
EM_REAL EM_COIL::INTEG_DELTA = 1.;
EM_COIL::IntegParam *EM_COIL::INTEG_PARAM = new EM_COIL::IntegParam(1, 5, 5, 3);
const EM_Symmetry *EM_COIL::symmetry = NULL;
EM_REAL EM_COIL::CONDUCTIVITY = 0.;

//int EM_COIL::no_COIL=0;
bool EM_COIL::active = false;

EM_COIL::EM_COIL(int num)
:numElements(num), elements(NULL)
//	,potentials(NULL)
, m_A_sources(NULL)
, A_source(NULL), B_source(NULL)
, B_source_at_nodes(NULL), B_source_at_element_centers(NULL)
, potentials(NULL)
, regularizationVector(NULL)
,selfInductance(0)
{
	type = EXTERNAL_CURRENT_TYPE;
	resistive = 0;
	//	joule_heat=0.;
	circuit_parm = new CIRCUIT_PARM;
	if (num) elements = (EM_COIL_Element **)CALLOC(num, sizeof(EM_COIL_Element *), "EM_COIL:: elements");
}

EM_COIL::~EM_COIL(){
	FREE(circuit_parm, "EM_COIL::circuit_parm");

	if (elements) {
		for (int i = 0; i<numElements; ++i)
			delete elements[i];
		FREE(elements, "EM_COIL:: elements");
	}

	if (A_source) FREE(A_source, "EM_COI::A_source");
	if (m_A_sources) delete m_A_sources;
	//	if(regularizationVector) delete regularizationVector;
}



void EM_COIL::Read(EM_TextFile *ffile) {
	TITLERD(ffile_control_input, "*COIL   * SERIES_NO  *  TIME_ID * NO_ELEMENTS * MOTIN_ID * IN_ROTOR *\n",
		"%s %d %d %d %d %d %d", ND7(namfld, &i1, &i2, &i3, &i4, &i5, &i6),
		(i1 >= 1 && i2 >= 0 && i3 >= 1 && i4 >= 0 && (i5 == 0 || i5 == 1) && i6 >= 0));
	int series_no = i1;
	int time_no = i2;
	//	no_COIL+=i3;
#ifdef SUPERPRO_PROTECT
	if (code_protection.version == VERSION_DEMO && no_COIL >COIL_NO_LIMIT_IN_DEMO) {
		fprintf(Stderr, "Protection key not found.\nWithout key, EMSolution allows to run only in demo version (max.20 COIL/max.10 NETWORK, until %d:%d:%d).\nFor a complete version of EMSolution, please contact us or your distributor.\n\n",
			YEAR_MAX_DEMO, MONTH_MAX_DEMO, DAY_MAX_DEMO);
		EXIT(1);
	}
#endif
	numElements = i3;
	if (CALCULATION::p()->COIL_MOTION == ON) motion.id = i4;
	else motion.id = 0;
	id = series_no;
	time_evolution.no = time_no;
	in_rotor = i5;

	if (i6 == 0) potential_region_no = 1;
	else potential_region_no = i6;

	static int coordinate_id;
	static EM_Coordinate *coordinate = EM_Coordinate::Global;
	symmetry = SYMMETRY::p();

	/*	-----  numElements (loop)  -----      */
	elements = (EM_COIL_Element **)CALLOC(numElements, sizeof(EM_COIL_Element *), "EM_COIL:: elements");
	int j = 0;
	while (j<numElements) {

		CARDRD_TITLE(ffile,
			"%s",
			namfld,
			(1));
		if (equal_string(namfld, "COORDID")) {
			TITLERD(ffile, "*COORDID * COORDINATE_ID *\n",
				"%s %d", ND2(namfld, &i1),
				(i1 >= 0));
			coordinate_id = i1;
			coordinate = EM_Coordinate::coordinate_pointer(coordinate_id);
		}
		else if (equal_string(namfld, "LINE_APPROX")){
			TITLERD(ffile, "*LINE_APPROX * LIMIT *\n",
				"%s %lf", ND2(namfld, &r1),
				(r1 >= 0.));
			if (r1) LINE_APP_LIMIT = r1; /* default: 1000 */
		}
		else if (equal_string(namfld, "DELTA")){
			TITLERD(ffile, "*DELTA   * DELTA(m) *\n",
				"%s %lf", ND2(namfld, &r1),
				(r1>0.));
			INTEG_DELTA = r1; /* default: 1.0 */
		}
		else if (equal_string(namfld, "INTEG")){
			TITLERD(ffile, "*INTEG   * NDIV *  INTX *  INTY  * INTZ  *\n",
				"%s %d  %d  %d  %d", ND5(namfld, &i1, &i2, &i3, &i4),
				(i1 >= 0 && i2 >= 0 && i3 >= 0 && i4 >= 0));
			delete INTEG_PARAM;
			INTEG_PARAM = new IntegParam(i1, i2, i3, i4);
		}
		else if (equal_string(namfld, "SYMMETRY")){
			TITLERD(ffile, "*SYMMETRY * XSYM * YSYM * ZSYM  **\n",
				"%s %d  %d  %d", ND4(namfld, &i1, &i2, &i3),
				(i1 == 0 || i1 == 1 || i1 == -1)
				&& (i2 == 0 || i2 == 1 || i2 == -1)
				&& (i3 == 0 || i3 == 1 || i3 == -1));
			symmetry = new EM_Symmetry(i1, i2, i3);
		}
		else if (equal_string(namfld, "CONDUCTIVITY")){
			TITLERD(ffile, "*SYMMETRY * (S/m)  **\n",
				"%s %lf", ND2(namfld, &r1),
				(r1 >= 0.));
			CONDUCTIVITY = r1;
		}
		else {

			bool integralVolume = false;
			EM_COIL_Element *coil = elements[j];
			if (equal_string(namfld, "UNIF")) {
				coil = new UNIFORM_FIELD;
			}
			else if (equal_string(namfld, "LINE")) {
				coil = new INFINITE_LINE;
			}
			else if (equal_string(namfld, "LOOP")) {
				coil = new EM_COIL_Loop;
			}
			else if (equal_string(namfld, "GCE")) {
				coil = new EM_COIL_Block;
			}
			else if (equal_string(namfld, "ARC")) {
				coil = new EM_COIL_Arc;
			}
			else if (equal_string(namfld, "FGCE")) {
				coil = new FGENERAL_CURRENT_ELEMENT;
			}
			else if (equal_string(namfld, "FARC")) {
				coil = new FARC_ELEMENT;
			}
			else if (equal_string(namfld, "CONSTA")) {
				coil = new CONST_A_FIELD;
			}
			else if (equal_string(namfld, "MESH")){
				coil = new EM_MeshedCOIL;
			}
			else if (equal_string(namfld, "DIPO")) {
				coil = new DIPOLE_FIELD;
			}
			else if (equal_string(namfld, "LOOP-")) {
				coil = new EM_COIL_Loop;
				integralVolume = true;
			}
			else if (equal_string(namfld, "GCE-")) {
				coil = new EM_COIL_Block;
				integralVolume = true;
			}
			else if (equal_string(namfld, "ARC-")) {
				coil = new EM_COIL_Arc;
				integralVolume = true;

			}
			else if (equal_string(namfld, "MESH-")) {
				coil = new EM_MeshedCOIL;
				integralVolume = true;
			}
			else if (equal_string(namfld, "RECT")) {
				coil = new EM_RectangleCOIL;

			}
			else if (equal_string(namfld, "RACE")) {
				coil = new EM_RaceTrackCOIL;

			}
			else {
				card_error(ffile);
			}

			coil->Read(ffile);
			if (integralVolume) coil->ReadIntegralParameters(ffile);
			coil->SetConductivity(CONDUCTIVITY);

			//			*(external_field->typs+ j) = field_id;
			elements[j] = coil;
			coil->SetCoordinate(coordinate);
			++j;
		}
	}

	CARDRD(ffile, "*END of COIL *\n", "%s", namfld, (equal_string(namfld, "END")));
}

void EM_COIL::IntegParam::Read(EM_TextFile *ffile){
	CARDRD(ffile, "* NDIV * INT_X * INT_Y * INT_Z *\n",
		"%d %d %d %d", ND4(&i1, &i2, &i3, &i4),
		((i1>0 || i1 == -1) && i2 >= 0 && i3 >= 0 && i4 >= 0));

	if (i1 == -1){
		*this = *INTEG_PARAM;
	}
	else {
		integral_volume_division = i1;
		intx = i2;
		inty = i3;
		intz = i4;
	}
}

void EM_COIL::Initialize(){
	for (int j = 0; j<numElements; ++j){
		elements[j]->Initialize();
	}
}

int EM_COIL::NumIntegralVolume() const{
	int count = 0;
	for (int j = 0; j<numElements; ++j){
		const EM_COIL_Element *elem = elements[j];
		count += elem->NumIntegralVolume();
	}
	return count;
}

EM_REAL EM_COIL::CalcResistance() const{
	EM_REAL resistance = 0.;
	for (int i = 0; i<numElements; ++i){
		if (elements[i]->GetConductivity())
			resistance += elements[i]->CalcResistance();
	}
	return resistance;
}

const EM_REAL EM_COIL::GetNumTurns() const { return Element(0)->GetCurrent(); }

EM_3dVector EM_COIL_Element::Global_A_Field(const EM_3dPoint &x) const{
	EM_3dVector A;
	if (!IsFieldActive()) {
		A.Zero();
		return A;
	}

	const EM_Coordinate *coordinate = GetCoordinate();
	EM_3dPoint xlocal = coordinate->ToLocal(x);

	EM_3dVector alocal = A_Field(xlocal);

	EM_Coordinate::Helper helper(coordinate, xlocal);

	A = coordinate->ToGlobal(alocal, helper);
	return A;
}

EM_3dVector EM_COIL_Element::Global_B_Field(const EM_3dPoint &x) const{
	EM_3dVector A;
	if (!IsFieldActive()) {
		A.Zero();
		return A;
	}

	const EM_Coordinate *coordinate = GetCoordinate();
	EM_3dPoint xlocal = coordinate->ToLocal(x);

	EM_3dVector alocal = B_Field(xlocal);

	EM_Coordinate::Helper helper(coordinate, xlocal);

	A = coordinate->ToGlobal(alocal, helper);
	return A;
}

void  EM_COIL_Element::Global_A_B_Field(const EM_3dPoint &x, EM_3dVector &A, EM_3dVector &B) const{
	if (!IsFieldActive()) {
		A.Zero();
		B.Zero();
		return;
	}

	const EM_Coordinate *coordinate = GetCoordinate();
	EM_3dPoint xlocal = coordinate->ToLocal(x);

	EM_3dVector alocal, blocal;
	A_B_Field(xlocal, alocal, blocal);

	EM_Coordinate::Helper helper(coordinate, xlocal);

	A = coordinate->ToGlobal(alocal, helper);
	B = coordinate->ToGlobal(blocal, helper);
}

EM_3dVector EM_COIL::A_Field(const EM_3dPoint &x, const MOTION *targetMotion) const{
	EM_3dVector A;

	EM_3dPoint x_tmp2 = x;
	if (CALCULATION::p()->COIL_MOTION == ON) move_point_global(x_tmp2, MOTION::GlobalMotion());

	EM_3dVector x_tmp;
	EM_3dVector x_tmp1;

	translated_position_by_motion(x_tmp2, x_tmp1, motion.ptr);

	if(targetMotion)
		reverse_translated_position_by_motion(x_tmp1, x_tmp, targetMotion);
	else
		x_tmp = x_tmp1;

	A.Zero();

	for (int i = 0; i<numElements; ++i) {
		const EM_COIL_Element *elem = elements[i];
		A += elem->Global_A_Field(x_tmp);
	}

	translate_vector_by_motion(A, motion.ptr);

	const MOTION *global_motion = MOTION::GlobalMotion();
	if (CALCULATION::p()->COIL_MOTION == ON && global_motion && global_motion->phiid.ptr)
		rotate_vector(A, global_motion->coordinate, global_motion->cosphi, global_motion->sinphi);

	return A;
}

EM_3dVector EM_COIL::B_Field(const EM_3dPoint &x, const MOTION *targetMotion) const{
	EM_3dVector A;

	EM_3dPoint x_tmp2 = x;
	if (CALCULATION::p()->COIL_MOTION == ON) move_point_global(x_tmp2, MOTION::GlobalMotion());

	EM_3dVector x_tmp1;
	EM_3dVector x_tmp;
	translated_position_by_motion(x_tmp2, x_tmp1, motion.ptr);

	if(targetMotion)
		reverse_translated_position_by_motion(x_tmp1, x_tmp, targetMotion);
	else
		x_tmp = x_tmp1;

	A.Zero();

	for (int i = 0; i<numElements; ++i) {
		const EM_COIL_Element *elem = elements[i];
		A += elem->Global_B_Field(x_tmp);
	}

	translate_vector_by_motion(A, motion.ptr);


	const MOTION *global_motion = MOTION::GlobalMotion();

	if (CALCULATION::p()->COIL_MOTION == ON && global_motion && global_motion->phiid.ptr)
		rotate_vector(A, global_motion->coordinate, global_motion->cosphi, global_motion->sinphi);

	/*
	if(A.Length()>1.) {
	fprintf(ffile_check->file, " point %15.5e  %15.5e  %15.5e  A  %15.5e  %15.5e  %15.5e\n",
	x.x, x.y, x.z, A.x[0], A.x[1], A.x[2]);
	fflush(ffile_check->file);
	}
	*/


	return A;
}

void EM_COIL::A_B_Field(const EM_3dPoint &x, EM_3dVector &A, EM_3dVector &B, const MOTION *targetMotion) const{

	EM_3dPoint x_tmp2 = x;
	if (CALCULATION::p()->COIL_MOTION == ON) move_point_global(x_tmp2, MOTION::GlobalMotion());

	EM_3dVector x_tmp1;
	EM_3dVector x_tmp;
	translated_position_by_motion(x_tmp2, x_tmp1, motion.ptr);

	if (targetMotion)
		reverse_translated_position_by_motion(x_tmp1, x_tmp, targetMotion);
	else
		x_tmp = x_tmp1;

	A.Zero();
	B.Zero();

	for (int i = 0; i<numElements; ++i) {
		const EM_COIL_Element *elem = elements[i];
		EM_3dVector aglobal, bglobal;
		elem->Global_A_B_Field(x_tmp, aglobal, bglobal);

		A += aglobal;
		B += bglobal;
	}

	translate_vector_by_motion(A, motion.ptr);
	translate_vector_by_motion(B, motion.ptr);

	const MOTION *global_motion = MOTION::GlobalMotion();
	if (CALCULATION::p()->COIL_MOTION == ON && global_motion && global_motion->phiid.ptr) {
		rotate_vector(A, global_motion->coordinate, global_motion->cosphi, global_motion->sinphi);
		rotate_vector(B, global_motion->coordinate, global_motion->cosphi, global_motion->sinphi);
	}

}

void EM_COIL::TotalingA() const{

	static_cast<const EM_A_SourceMeshSetPotentials *>(A_SourcePotentials())->m_totalA->Add(*A_SourceVector(), Amplitude());

}


EM_3dVector *EM_COIL::CalcForce(EM_REAL current, EM_3dVector *force, const EM_MeshSet *mesh) const{
	for (int i = 0; i<numElements; ++i){
		const EM_COIL_Element *elem = elements[i];
		if (elem->IsIntegralVolume()){
			elem->CalcForce(current, motion.ptr, mesh, force);
			force += elem->NumIntegralVolume();
		}
	}
	return force;
}

void EM_COIL::PrintForce(FILE *file, int id, const EM_3dVector *force) const{
	for (int k = 0; k<numElements; ++k){
		const EM_COIL_Element *elem = elements[k];
		if (elem->IsIntegralVolume()) {
			elem->PrintForce(file, id, force);
			force += elem->NumIntegralVolume();
		}
	}
}

void EM_COIL::OutputForces(FILE *file, const EM_3dVector *&force, int &elementNo) const{
	for (int k = 0; k<numElements; ++k){
		const EM_COIL_Element *elem = elements[k];
		elem->OutputForces(file, force, elementNo);
	}
}

/****    EM_COIL_Element   ***********************************/
EM_COIL_Element::EM_COIL_Element()
:fieldActive(true), integParam(NULL)
, conductivity(0.)
{
	SetCoordinate(EM_Coordinate::Global);
}

void EM_COIL_Element::SetIntegral(int div, int intx, int inty, int intz){
	if (!div) div = GetDivision(EM_COIL::INTEG_DELTA);
	integParam = new EM_COIL::IntegParam(div, intx, inty, intz);
	fieldActive = false;
}
void EM_COIL_Element::SetIntegral(const EM_COIL::IntegParam *integ){
	SetIntegral(integ->integral_volume_division,
		integ->intx, integ->inty, integ->intz);
}

void EM_COIL_Element::ReadIntegralParameters(EM_TextFile *ffile){
	integParam = new EM_COIL::IntegParam();
	integParam->Read(ffile);
	if (integParam->integral_volume_division == 0)
		integParam->integral_volume_division = GetDivision(EM_COIL::INTEG_DELTA);
	fieldActive = false;
}

void EM_COIL_Element::OutputForces(FILE *file, const EM_3dVector *&force, int &elementNo) const{
	EM_PostFileType &output_file_type = *EM_PostFileType::GetType();
	if (IsIntegralVolume()){
		EM_REAL data[4];
		for (int j = 0; j<NumIntegralVolume(); ++j){
			data[0] = force->x[0];
			data[1] = force->x[1];
			data[2] = force->x[2];
			EM_REAL absF = SQRT((data[0] * data[0] + data[1] * data[1] + data[2] * data[2]));
			if (absF){
				data[3] = absF;
				output_file_type.write_eval_data(file, elementNo++,
					4, data, 8);
			}
			++force;
		}
	}
	else {
		elementNo += NumPlotElements();
	}
}

void EM_COIL_Element::PrintForce(FILE *file, int id, const EM_3dVector *force) const{
	for (int j = 0; j<NumIntegralVolume(); ++j){
		fprintf(file, " %4d  %4s%2d  %13.5e   %13.5e  %13.5e\n",
			id, Name(), j + 1, (*force)[0], (*force)[1], (*force)[2]);
		++force;
	}
}

/*
EM_Integrator *EM_COIL_Element::GetIntegrator() const{
}
*/
/****    CONST_A_FIELD   ***********************************/
void CONST_A_FIELD::Read(EM_TextFile *ffile){
	TITLERD(ffile, "*CONSTA  *   AX(T)    *    AY(T)   *   AZ(T)    *\n",
		"%s %lf %lf %lf", ND4(namfld, &r1, &r2, &r3),
		(1));
	A_vector.Set(r1, r2, r3);
}

EM_3dVector CONST_A_FIELD::A_Field(const EM_3dPoint &x) const{
	return A_vector;
}

EM_3dVector CONST_A_FIELD::B_Field(const EM_3dPoint &x) const{
	return EM_3dVector(0., 0., 0.);
}
#endif

/****    UNIFORM_FIELD   ***********************************/
UNIFORM_FIELD::UNIFORM_FIELD(double bx, double by, double bz, int dir) {
	B_vector.Set(bx, by, bz);
	A_direction = dir;
}

EM_3dVector UNIFORM_FIELD::A_Field(const EM_3dPoint &x) const{
	EM_3dVector A, B;
	B = B_vector;

	switch (A_direction) {
	case 0:
		A = EM_CrossProduct(B, EM_3dVector(x)) / 2.;
		break;

	case 1:
		A[0] = B[1] * x[2];
		A[1] = B[2] * x[0];
		A[2] = B[0] * x[1];
		break;

	case 2:
		A[0] = -B[2] * x[1];
		A[1] = -B[0] * x[2];
		A[2] = -B[1] * x[0];
		break;
	}
	return A;
}

EM_3dVector UNIFORM_FIELD::B_Field(const EM_3dPoint &x) const{
	return B_vector;
}

INFINITE_LINE::INFINITE_LINE(double current, double px, double py, double pz
	, double dx, double dy, double dz) {
	this->current = current;
	A_direction = 0;
	on_point.Set(px, py, pz);
	direction.Set(dx, dy, dz);
	direction.Unitize();
}

EM_3dVector INFINITE_LINE::A_Field(const EM_3dPoint &x) const{
	/*******************************************************************************
	c  calculate magnetic vector potential of infinite line current
	c
	c     (xl,yl,zl) ------>----->--------
	c                   i(a)    (dx,dy,dz)
	c
	c******************************************************************************/
	EM_3dVector A;
	double xx[3], dx[3], aj, r, dd;
	double minimum_r = 1.e-10;
	int i;

	r = 0.;
	dd = 0.;
	for (i = 0; i<3; ++i)
	{
		xx[i] = x[i] - (double)on_point[i];
		dx[i] = direction[i];

		r += xx[i] * xx[i];
		dd += xx[i] * dx[i];
	}
	aj = (double)current * 2.e-7;

	r = SQRT(r - dd*dd);
	if (r< minimum_r) r = minimum_r;

	if (!A_direction)
	{
		r = LOG(r);

		for (i = 0; i<3; ++i)
			A[i] = -r * dx[i] * aj;
	}
	else
	{
		for (i = 0; i<3; ++i)
		{
			xx[i] -= dd*dx[i];
			A[i] = xx[i] / (r*r)*dd *aj;
		}
	}
	return A;
}

EM_3dVector INFINITE_LINE::B_Field(const EM_3dPoint &x) const{
	EM_3dVector B;
	double xx[3], dx[3], cx[3], aj, r2;
	double minimum_r2 = 1.e-20;
	int i, k, l;

	for (i = 0; i<3; ++i)
	{
		xx[i] = x[i] - on_point[i];
		dx[i] = direction[i];
	}
	aj = (double)current * 2.e-7;


	r2 = 0.;
	for (i = 0; i<3; ++i)
	{
		k = i + 1;
		l = i + 2;
		if (k>2) k -= 3;
		if (l>2) l -= 3;
		cx[i] = dx[k] * xx[l] - dx[l] * xx[k];
		r2 += cx[i] * cx[i];
	}
	if (r2<minimum_r2) r2 = minimum_r2;

	for (i = 0; i<3; ++i)
		B[i] = cx[i] / r2 * aj;
	return B;
}
#ifndef COILBASE
void INFINITE_LINE::A_B_Field(const EM_3dPoint &x, EM_3dVector &A, EM_3dVector &B) const{

	double minimum_r = 1.e-10;
	double minimum_r2 = 1.e-20;

	EM_3dVector xx = x - on_point;

	EM_REAL aj = current * 2.e-7;

	EM_3dVector cx = EM_CrossProduct(direction, xx);
	EM_REAL r2 = cx*cx;
	EM_REAL r = SQRT(r2);

	if (r< minimum_r){
		r = minimum_r;
		r2 = minimum_r2;
	}
	B = cx*(aj / r2);

	if (!A_direction)
	{
		r = LOG(r);
		A = direction*(-r*aj);

	}
	else
	{
		EM_REAL dd = xx*direction;
		xx -= dd*direction;
		A = xx*(dd*aj / r2);
	}

}



/****    FGENERAL_CURRENT_ELEMENT   ***********************************/
void FGENERAL_CURRENT_ELEMENT::Read(EM_TextFile *ffile){
	TITLERD(ffile, "*FGCE   * VARIABLE::CURRENT(A)*\n",
		"%s %lf", ND2(namfld, &r1),
		(1));
	current = r1;

	CARDRD(ffile, "*            XS(m)  *   YS(m)  *  ZS(m)  *  XE(m)  *  YE(m)  *  ZE(m)  *\n",
		"%lf %lf %lf %lf %lf %lf", ND6(&r1, &r2, &r3, &r4, &r5, &r6),
		(!(r1 == r4 && r2 == r5 && r3 == r6)));
	start_point.Set(r1, r2, r3);
	end_point.Set(r4, r5, r6);
}

void FGENERAL_CURRENT_ELEMENT::CountMeshData(int &no_nodes, int &no_elements) const{
	no_elements += 1;
	no_nodes += 2;
}
#if 0
EM_3dVector FGENERAL_CURRENT_ELEMENT::A_Field(const EM_3dPoint &x) const{
	EM_3dVector A, B;
	//   double   parameter_of_fgce[7];
	// TODO616
	//   double B[3];
	/*
	int i, ia=1, ib=0;

	for(i=0; i<3; ++i)
	{
	parameter_of_fgce[i]   = start_point[i];
	parameter_of_fgce[i+3] = end_point[i];
	}
	parameter_of_fgce[6] = (double) current;
	*/
	field_param param;
	param.calcA = 1;
	param.calcB = 0;

	CalcFields(x, A, B, &param);
	//	const EM_REAL *X=x;
	//   fldfgce_(X, parameter_of_fgce, &ia, &ib, A, B);
	return A;
}

EM_3dVector FGENERAL_CURRENT_ELEMENT::B_Field(const EM_3dPoint &x) const{
	EM_3dVector B;
	//	double   parameter_of_fgce[7];
	double A[3];
	//   int i, ia=0, ib=1;
	/*
	for(i=0; i<3; ++i)
	{
	parameter_of_fgce[i]   = start_point[i];
	parameter_of_fgce[i+3] = end_point[i];
	}
	parameter_of_fgce[6] = (double) current;
	*/
	field_param param;
	param.calcA = 0;
	param.calcB = 1;
	CalcFields(x, A, B, &param);
	//   fldfgce_(x, parameter_of_fgce, &ia, &ib, A, B);
	return B;
}
#endif
EM_3dVector FGENERAL_CURRENT_ELEMENT::A_Field(const EM_3dPoint &x) const{
	EM_3dVector A, B;
	double   parameter_of_fgce[7];
	// TODO616
	//   double B[3];

	int i, ia = 1, ib = 0;

	for (i = 0; i<3; ++i)
	{
		parameter_of_fgce[i] = start_point[i];
		parameter_of_fgce[i + 3] = end_point[i];
	}
	parameter_of_fgce[6] = (double)current;
	/*
	field_param param;
	param.calcA=1;
	param.calcB=0;

	CalcFields(x, A, B, &param);
	*/
	const EM_REAL *X = x;
	fldfgce_(X, parameter_of_fgce, &ia, &ib, A, B);
	return A;
}
EM_3dVector FGENERAL_CURRENT_ELEMENT::B_Field(const EM_3dPoint &x) const{
	EM_3dVector B;
	double   parameter_of_fgce[7];
	double A[3];
	int i, ia = 0, ib = 1;

	for (i = 0; i<3; ++i)
	{
		parameter_of_fgce[i] = start_point[i];
		parameter_of_fgce[i + 3] = end_point[i];
	}
	parameter_of_fgce[6] = (double)current;
	/*
	field_param param;
	param.calcA=0;
	param.calcB=1;
	CalcFields(x, A, B, &param);
	*/
	fldfgce_(x, parameter_of_fgce, &ia, &ib, A, B);
	return B;
}

/****    FARC_ELEMENT   ***********************************/
void FARC_ELEMENT::Read(EM_TextFile *ffile){
	TITLERD(ffile, "*ARC    *  current(A)  *\n",
		"%s %lf", ND2(namfld, &r1),
		(1));
	current = r1;

	CARDRD(ffile, "*          X(m)  *  Y(m)  *   Z(m)  * RADIUS(m)*\n",
		"%lf %lf %lf %lf", ND4(&r1, &r2, &r3, &r4),
		(r4>0.));

	center.Set(r1, r2, r3);
	radius = r4;

	CARDRD(ffile, "*          ALPHA(deg)* BETA(deg) * PHI1(deg)   * PHI2(deg) *\n",
		"%lf %lf %lf %lf", ND4(&r1, &r2, &r3, &r4),
		(r3< r4));
	Euler_angle_alpha = r1;
	Euler_angle_beta = r2;
	start_angle = r3;
	end_angle = r4;

	Euler_angle_alpha *= PI / 180.;
	Euler_angle_beta *= PI / 180.;
	start_angle *= PI / 180.;
	end_angle *= PI / 180.;
}

void FARC_ELEMENT::CountMeshData(int &no_nodes, int &no_elements) const{
	EM_REAL dt = 40.*PI / 180.;
	EM_REAL phi1 = start_angle;
	EM_REAL phi2 = end_angle;
	do {
		dt /= 2.;
		int no_coil_elements = (int)((phi2 - phi1) / dt + 0.9999);
		dt = (phi2 - phi1) / no_coil_elements;
	} while (no_coil_elements <= MIN_ARC_DIVISIONS);
	no_elements += no_coil_elements;
	no_nodes += no_coil_elements + 1;
}

EM_3dVector FARC_ELEMENT::A_Field(const EM_3dPoint &x) const{
	EM_3dVector A;
	double   parameter_of_farc[9];
	double B[3];
	int i, ia = 1, ib = 0;

	for (i = 0; i<3; ++i)
		parameter_of_farc[i] = (double)center[i];
	parameter_of_farc[3] = (double)radius;
	parameter_of_farc[4] = (double)Euler_angle_alpha;
	parameter_of_farc[5] = (double)Euler_angle_beta;
	parameter_of_farc[6] = (double)start_angle;
	parameter_of_farc[7] = (double)end_angle;
	parameter_of_farc[8] = (double)current;

	fldfarc_(x, parameter_of_farc, &ia, &ib, A, B);
	return A;
}

EM_3dVector FARC_ELEMENT::B_Field(const EM_3dPoint &x) const{
	EM_3dVector B;
	double   parameter_of_farc[9];
	double A[3];
	int i, ia = 0, ib = 1;

	for (i = 0; i<3; ++i)
		parameter_of_farc[i] = (double)center[i];
	parameter_of_farc[3] = (double)radius;
	parameter_of_farc[4] = (double)Euler_angle_alpha;
	parameter_of_farc[5] = (double)Euler_angle_beta;
	parameter_of_farc[6] = (double)start_angle;
	parameter_of_farc[7] = (double)end_angle;
	parameter_of_farc[8] = (double)current;

	ia = 0;
	ib = 1;

	fldfarc_(x, parameter_of_farc, &ia, &ib, A, B);
	return B;
}

void EM_COIL::AssembleLoadVector(EM_MeshSetPotentials *basePotentials, EM_REAL &resistance) {

	if (CALCULATION::p()->COIL_MOTION == ON) {
		if ( !IsIncludedInCircuit()) return;
		if ( motion.id) {
			for (int k = 0; k < 3; k++)
				motion.ptr->xm[k] = 0;
			motion.ptr->cosphi = 1.;
			motion.ptr->sinphi = 0.;
		}
	}

	EM_A_SourceMeshSetPotentials *a_source = static_cast<EM_A_SourceMeshSetPotentials *>(A_SourcePotentials());

//	if (a_source) {
//	a_source->AssembleLoadVector(*this, basePotentials);
	{
		EM_StopWatches::p->Start(60);
		a_source->meshSetPotentials = basePotentials;
		EM_DomainVectors<EM_REAL> *domainCouplings = new EM_DomainVectors<EM_REAL>(a_source->DomainVectorProfile());
		SetLoadVector(domainCouplings);

		CoeffcientVector = new EM_DomainVectors<EM_REAL>(a_source->DomainVectorProfile());
		PARALLEL_FOR_NEW(imesh, a_source->numDomains) {
			const EM_Mesh *mesh = a_source->MeshSet()->Mesh(imesh);
			const EM_MeshPotentials *pot = (*basePotentials)[imesh];
			EM_BlockVector<double> *blockVector = CoeffcientVector->GetBlockVector(imesh);
			for (int i = 0; i<mesh->EdgeCount(); ++i) {
				const EM_EdgeRef & edge = mesh->Edge(i);
				const EM_PotentialRef &potential = pot->Potential(edge);
				if (potential) {
					const VARIABLE *variable = potential.Variable();
					if (a_source->Potential(edge)) {
						(*blockVector)[variable->GetNo()] = variable->coefficient;
					}
				}
			}
		}


		EM_COIL_Regularization *regularization = a_source->Regularization();

		if (regularization->On()) {
			a_source->regPotentials = static_cast<EM_COIL_RegularizationPotentials *>(regularization->Potentials());
			assert(a_source->regPotentials);
			EM_Vector<EM_REAL> *regLoadVector = regularization->GetLoadVector();
		}

		a_source->AssembleHtLoad(*this, *domainCouplings);

		if (regularization->On()) {
			regularization->Regularize(*this);
			regularization->AssembleLoadVector(static_cast<const EM_DomainVectors<double> *>(RegularizationVector()), *domainCouplings);
		}

		BnSourceMeshSetPotentials()->CalcA_SourceVector(*this);

		a_source->meshSetPotentials->SetValueVector(domainCouplings);

		a_source->AssembleAtLoadBySolid(*this, *domainCouplings);

		//	(*domainCouplings) *=EM_Circuit::region_multi;
		EM_StopWatches::p->Stop(60);
	}



		resistance = CalcResistance();
//	}
}

void EM_COIL::RecalculateLoadVector(){

	if (!motion.id) return;

	EM_A_SourceMeshSetPotentials &a_source = static_cast<EM_A_SourceMeshSetPotentials &>(*A_SourcePotentials());

	previousLoadVector = GetLoadVector()->DeepCopy();

	a_source.RecalculateLoadVector(*this);

	(*GetLoadVector()) *= EM_Circuit::region_multi;

}


void EM_FieldSources::SourcePotentialInductances(const EM_MeshSetPotentials &systemPotentials) {

	EM_DomainVectors<EM_REAL> A(&systemPotentials);
	const EM_MeshSet *meshSet = systemPotentials.MeshSet();
	int numDomains = meshSet->Count();

	EM_FieldSources &field_sources = *Analysis().FieldSources();
	//	EM_REAL *inductances=field_sources.inductances;
	int no_data = field_sources.NumSources();
	for (int isrc0 = 0; isrc0<no_data; ++isrc0) {
		const EM_FieldSource *source1 = &field_sources[isrc0];
		if (source1&& source1->type == EXTERNAL_CURRENT_TYPE) {
			int isrc = source1->seq_no;
			const EM_COIL *external1 = static_cast<const EM_COIL *>(source1);

			external1->SourcePotential(systemPotentials, A);

			for (int jsrc0 = 0; jsrc0 <= isrc0; ++jsrc0) {
				const EM_FieldSource *source2 = &field_sources[jsrc0];
				if (source2 && source2->type == EXTERNAL_CURRENT_TYPE) {
					int jsrc = source2->seq_no;
					const EM_COIL *external2 = static_cast<const EM_COIL *>(source2);

					const EM_DomainVectors<EM_REAL> *vec0 = external2->GetLoadVector();
					if (vec0) {
						EM_DomainVectors<EM_REAL> *vec = new EM_DomainVectors<EM_REAL>(*vec0);
						*vec += *external2->GetLoadVector();
						vec->MultiplyElement(*external2->CoeffcientVector);

						if (external1->potential_region_no != external2->potential_region_no) continue;
							//					Inductance(isrc,jsrc) -= *external2->GetLoadVector() * A;
							//COIL_Inductance(isrc, jsrc) -= *vec * A;

						delete vec;
					}
				}
			}
		}
	}
}

void EM_COIL::SourcePotential(const EM_MeshSetPotentials &potentials, EM_DomainVectors<EM_REAL> &A) const{
	const EM_A_SourceMeshSetPotentials *srcA_Potentials = static_cast<const EM_A_SourceMeshSetPotentials *>(A_SourcePotentials());//Set[potential_region_no-1];
	if (!srcA_Potentials) return;
	const EM_MeshSetPotentials &srcPotentials = *srcA_Potentials->BnSourceMeshSetPotentials();
	const EM_MeshSet *meshSet = potentials.MeshSet();
	int numDomains = potentials.NumDomains();
	A.Zero();

	PARALLEL_FOR(imesh, numDomains, parallel for){
		//	for(int imesh=0; imesh<numDomains; ++imesh){
		const EM_Mesh *mesh = meshSet->Mesh(imesh);
		const EM_MeshPotentials &meshPotentials = *potentials[imesh];
		EM_BlockVector<EM_REAL> &v = *A.GetBlockVector(imesh);
		const EM_Vector<EM_REAL> &A_source = *A_SourceVector()->GetBlockVector(imesh)/*->Vector(0)->GetContent()*/;

		for (int i = 0; i<mesh->EdgeCount(); ++i){
			const EM_EdgeRef &edge = mesh->Edge(i);
			if (!edge->OnPotentialInterface()) continue;
			const EM_PotentialRef &potential = meshPotentials.Potential(edge);

			if (potential) {
				EM_REAL a = 0.;
				const EM_PotentialRef &srcPotential = srcPotentials.Potential(edge);
				if (srcPotential){
					for (int j = 0; j<srcPotential->NumVariables(); ++j) {
						const VARIABLE *variable1 = srcPotential.Variable(j);
						int no1 = variable1->GetNo();
						a = A_source[no1];
					}
				}

				for (int j = 0; j<potential->NumVariables(); ++j) {
					const VARIABLE *variable1 = potential.Variable(j);
					if (variable1->GetType() == VARIABLE::A_EDGE) {
						int no1 = variable1->GetNo();
						v[no1] = a;
					}
				}
			}
		}
	}
	A.Accumulated(true);
}


void set_COIL_connection(CONNECTION_VECTOR *connection, const EM_MeshSet *mesh){
	int n, e, i;


	//	v_element=v_elements;
	for (n = 0; n<mesh->SolidCount(); ++n){
		const EM_SolidRef &v_element = mesh->Solid(n);
		if (v_element->contacting_to_potential_interface){
			for (e = 0; e<v_element.NumEdges(); ++e){
				const EM_EdgeRef &edge = v_element.Edge(e);
				assert(false);
				//				add_connection(edge->Potential(), VARIABLE::A_EDGE, connection);
			}
		}
	}
	EM_FieldSources &field_sources = *Analysis().FieldSources();
	for (i = 0; i<field_sources.NumSources(); ++i){
		const EM_FieldSource &source = field_sources[i];
		if (source.type == EXTERNAL_CURRENT_TYPE){
			add_connection(source.circuit_parm->Potential(), VARIABLE::CURRENT, connection);
		}
	}

}

void put_COIL_connection_in_list(EM_PotentialConnection *connections, EM_FieldSources field_sources,
	const EM_MeshSet *mesh){
	LINKED_LIST_INT *lists = connections->Lists();
	CONNECTION_VECTOR *connection = connections->PotentialIndices();
	LINKED_LIST_INT list0(connections->Allocator());
	int i, j, n, init = 1, row;
	const VARIABLE *variable;
	//	LINKED_LIST_INT *list;

	list0.first = NULL;
	/*	list0 =   (LINKED_LIST_INT *) CALLOC(1,  sizeof(LINKED_LIST_INT ), "list0" );*/
	for (i = 0; i<field_sources.NumSources(); ++i){
		const EM_FieldSource &source = field_sources[i];
		if (source.type == EXTERNAL_CURRENT_TYPE){
			const EM_PotentialRef &potential = source.circuit_parm->Potential();
			if (potential){
				for (j = 0; j<potential->NumVariables(); ++j) {
					variable = potential.Variable(j);
					if (variable->GetType() == VARIABLE::CURRENT) {
						row = variable->GetNo();
						LINKED_LIST_INT *list = lists + row;
						if (init){
							connection->clear();
							set_COIL_connection(connection, mesh);
							//							const EM_PotentialRef &potential=source->circuit_parm->Potential();

							/*							list0=list;*/
							for (n = 0; n<connection->size() - 1; ++n){
								put_sorted(&list0, (*connection)[n]);
							}
							init = 0;
						} /*else {*/
						if (!list->first){
							copy_linkable_int2(list, &list0, row + 1);  // XXNO
						}
						/*						}*/
					}
				}
			}
		}
	}
	/*	FREE(list0, "list0");*/
	/*	free_list(list0);*/
}
#endif




