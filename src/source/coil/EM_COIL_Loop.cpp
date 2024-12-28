#ifdef COILBASE
#include "source/coil/EM_COIL.h"
#else
#include "basic/EM_defines.h"

#endif
#include "source/coil/coil_field.h"
#include "basic/my_math.h"

/****    EM_COIL_Loop   ***********************************/
EM_COIL_Loop::EM_COIL_Loop(EM_REAL current, EM_REAL radius, EM_REAL z, EM_REAL radialWidth, EM_REAL axialWidth)
{
	this->current = current;
	this->center=EM_3dPoint(0.,0.,z);
	this->radius=radius;
	radial_width=radialWidth;
	axial_width=axialWidth;
	Euler_angle_alpha=0.;
	Euler_angle_beta=0.;
	start_angle=-PI;
	end_angle=PI;

}
#ifndef COILBASE
void EM_COIL_Loop::Read(EM_TextFile *ffile){
    TITLERD(ffile, "*LOOP  * VARIABLE::CURRENT(A) * RADIUS(m)* CENTER_Z(m)* RADIAL_W(m)*  AXIAL_W(m) *\n",
             "%s %lf %lf %lf %lf %lf", ND6(namfld, &r1,&r2, &r3, &r4, &r5),
            ( r2>0. && r4>=0. && r5>=0.) );
	current  = r1;
	radius   = r2;
	center = EM_3dPoint(0.,0.,r3);
	radial_width = r4;
	axial_width = r5;
	Euler_angle_alpha=0.;
	Euler_angle_beta=0.;
	start_angle=-PI;
	end_angle=PI;
}

EM_REAL EM_COIL_Loop::GetLength() const{
	return 2.*PI*radius;
}
#endif
EM_3dVector EM_COIL_Loop::A_Field( const EM_3dPoint &x) const{
	EM_3dVector A;
	double minimum_R = 1.e-10;

	double aicoil= current;
    double rcoil = radius;
    double zcoil = center.z;
    double acoil = radial_width;
    double bcoil=  axial_width;

    double r2=x[0]*x[0]+x[1]*x[1];
	double r=SQRT(r2);
    if(r< minimum_R)
		A.Zero();
    else {
		double z=x[2];
		double p;
        if(!acoil || !bcoil)    
			p = psi_(r,z,rcoil,zcoil, aicoil);
        else {
			field_param param;
			double br, bz;
			param.calcA=1;
			param.calcB=0;
			field_LOOP(r, z, rcoil, zcoil, acoil, bcoil, aicoil, &p, &br, &bz, &param);
		}

         p /= r2;
         A[0] = -p *x[1];
         A[1] =  p *x[0];
         A[2] = 0.;
      }
	  return A;
}

EM_3dVector EM_COIL_Loop::B_Field( const EM_3dPoint &x) const{
	EM_3dVector B;
    double minimum_R = 1.e-10;


	double aicoil= current;
    double rcoil = radius;
    double zcoil = center.z;
    double acoil = radial_width;
    double bcoil=  axial_width;

    double r=SQRT(x[0]*x[0]+x[1]*x[1]);

    double z=x[2];
	double br,bz;
    if(!acoil || !bcoil)    
		bring0_(r,z,rcoil,zcoil, aicoil, br, bz);
    else {
		field_param param;
		param.calcA=0;
		param.calcB=1;
		double p;
		field_LOOP(r, z, rcoil, zcoil, acoil, bcoil, aicoil, &p, &br, &bz, &param);
	}

    if(r< minimum_R){
         B[0]=0.;
         B[1]=0.;
         B[2]= bz;
     }
     else{
         B[0] =  br *x[0] /r;
         B[1] =  br *x[1] /r;
         B[2] =  bz;
     }
	 return B;
}
#ifndef COILBASE
void EM_COIL_Loop::A_B_Field(const EM_3dPoint &x, EM_3dVector &A, EM_3dVector &B) const{
	double minimum_R = 1.e-10;

	double aicoil= current;
    double rcoil = radius;
    double zcoil = center.z;
    double acoil = radial_width;
    double bcoil=  axial_width;

	double p, br, bz;
    double r2=x[0]*x[0]+x[1]*x[1];
	double r=SQRT(r2);

	double z=x[2];
    if(!acoil || !bcoil) {  
		p = psi_(&r,&z,&rcoil,&zcoil, &aicoil);
		bring0_(&r,&z,&rcoil,&zcoil, &aicoil, &br, &bz);
	}
    else {
		field_param param;
		param.calcA=1;
		param.calcB=1;
		field_LOOP(r, z, rcoil, zcoil, acoil, bcoil, aicoil, &p, &br, &bz, &param);
	}

	if(r< minimum_R){
		A.Zero();
		B[0]=0.;
		B[1]=0.;
		B[2]= bz;
	}
	else{
		p /= r2;
		A[0] = -p *x[1];
		A[1] =  p *x[0];
		A[2] = 0.;

		B[0] =  br *x[0] /r;
		B[1] =  br *x[1] /r;
		B[2] =  bz;
	}

}

#endif
extern int PerfectEllepticIntegral(EM_REAL k2, EM_REAL kd2, EM_REAL *K, EM_REAL *E);
#define ZOT 1.e-10

void EM_COIL_Loop::field_LOOP(EM_REAL R, EM_REAL Z, EM_REAL CR, EM_REAL CZ, EM_REAL CA, EM_REAL CB, EM_REAL CI, 
			EM_REAL *PSI, EM_REAL *BR, EM_REAL *BZ, field_param *param) const{

	EM_REAL XI, A, AMRSQ, TUR, DFT, RZSQ, ASRZS, APR, APRSQ, RAZSQ, CASQ, CPSQ, RAZ, RAZ3;
	EM_REAL FEI, EEI, ALXU, ALXL;
	EM_REAL C, area2;
	int i,J;
/*	static EM_REAL ZOT=1.e-10;*/
	field_sum fsum;
	EM_REAL AU,AUSQ,AL,ALSQ,ARU,ARL,ZZU,ZUSQ,ZZL,ZLSQ,ARZ[4],RHO,
		RHOSQ,AA,ASQ,TAR,Z1,ZSQ,ARZSQ,AREA;

    AREA=CA*CB;
    XI=CI/AREA*1.e-7;  /* myu0*I/(4*PI) */
    RHO=R;
    Z1=Z-CZ;
    RHOSQ=RHO*RHO;
    ZSQ=Z1*Z1;

    if(ZOT*ZSQ >= RHOSQ) {
		RHOSQ=0.;
		RHO=0.;
	} 
    A=CR;
    ASQ=A*A;
      AMRSQ=(A-RHO)*(A-RHO);
      TUR=RHO+RHO;
      TAR=A*TUR;
      ZSQ=Z1*Z1;
      ARZSQ=AMRSQ+ZSQ;
      DFT=LINE_APP_LIMIT*maximum(CA,CB)/2.;
	  DFT= DFT*DFT;

	param->RHO=RHO;
	param->RHOSQ=RHOSQ;
	param->ZSQ=ZSQ;
	param->ASQ=ASQ;
	param->AA=CR;
	param->ARZSQ=ARZSQ;
	param->TAR=TAR;
	param->AREA=AREA;
	param->Z1=Z1;

	if(ARZSQ >= DFT) {
		RZSQ=RHOSQ+ZSQ;
		ASRZS=ASQ+RZSQ;
		APR=A+RHO;
		APRSQ=APR*APR;
		RAZSQ=APRSQ+ZSQ;
		CASQ=(TAR+TAR)/RAZSQ;
		CPSQ=ARZSQ/RAZSQ;
		area2=AREA+AREA;
		if(  PerfectEllepticIntegral(CASQ, CPSQ, &FEI, &EEI)  ){
			RAZ=sqrt(RAZSQ);
			if(param->calcA)
				fsum.ATHWAV=area2*(RAZ/RHO)*(FEI*(ASRZS/RAZSQ)-EEI);
			if(param->calcB) {
				C=area2/RAZ;
				fsum.BZWAVE=C*(FEI+((ASQ-RZSQ)/ARZSQ)*EEI);
				fsum.BRHO=C*(Z1/RHO)*(EEI*(ASRZS/ARZSQ)-FEI);
			}
		} else {
			RAZ=sqrt(RAZSQ);
			RAZ3=RAZ*RAZ*RAZ;
			C=area2*PI*ASQ/RAZ3;
			fsum.ATHWAV=C*0.5*RHO;
/*			ATHWAV=-CONS4*AREA*(A/sqrt(RAZSQ)); */
			if(param->calcB) {
/*				C=-ATHWAV/RAZSQ; */
				fsum.BRHO=C*1.5*RHO*Z1/RAZSQ;
				fsum.BZWAVE=C;
			}
		}
	} else {
		AU=A+CA/2.;
		AUSQ=(AU-RHO)*(AU-RHO);
		AL=A-CA/2.;
		ALSQ=(AL-RHO)*(AL-RHO);
		ARU=TUR*AU;
		ARL=TUR*AL;
		ZZU=Z1-CB/2.;
		ZUSQ=ZZU*ZZU;
		ZZL=Z1+CB/2.;
		ZLSQ=ZZL*ZZL;
		ARZ[0]=ALSQ+ZLSQ;
		ARZ[1]=ALSQ+ZUSQ;
		ARZ[2]=AUSQ+ZLSQ;
		ARZ[3]=AUSQ+ZUSQ;
/*		AA=A;*/

		if(RHO) {
			fsum.ATHWAV=0.;
			fsum.BRHO=0.;
			fsum.BZWAVE=0.;
			param->symmetricArc=1;

			param->AL=AL;
			param->AU=AU;
			param->ALSQ=ALSQ;
			param->AUSQ=AUSQ;
			param->ZLSQ=ZLSQ;
			param->ZUSQ=ZUSQ;
			param->ZZL=ZZL;
			param->ZZU=ZZU;
			for(i=0; i<4; ++i) param->ARZ[i]=ARZ[i];
			param->ARU=TUR*AU;
			param->ARL=TUR*AL;

			arc_integration(TwoDegree,HalfPi, param, &fsum);
			arc_integration(0.,TwoDegree, param, &fsum);

			if((ARZSQ+TAR) >= DFT)  {
				AA=A;
				field_ARC_line(HalfPi,PI, param, &fsum);
			}
			else 
				arc_integration(HalfPi,PI, param, &fsum);

		} else {
			fsum.ATHWAV=0.;
			if(param->calcB) {
				fsum.BRHO=0.;
				for( J=0; J<4; ++J) {
					ARZ[J]=sqrt(ARZ[J]);
				}
				ALXU=log((AU+ARZ[3])/(AL+ARZ[1]));
				ALXL=log((AU+ARZ[2])/(AL+ARZ[0]));
				fsum.BZWAVE=TwoPi*(ZZL*ALXL-ZZU*ALXU);
			}
		}
	}
	if(param->calcB) {
		*BZ=XI*fsum.BZWAVE;
		*BR=XI*fsum.BRHO;
	}
	if(param->calcA) {
	  fsum.ATHWAV*=XI;
      *PSI=fsum.ATHWAV*R;
	}
}
#include <math.h>
#define ZERO(X) ( X + 1.e0) == 1.
//C% BRING0
void EM_COIL_Loop::bring0_(double R, double Z, double CR, double CZ, double CI,
	double& BR, double& BZ)
{
	//	PI = 3.1415926535898;
	double AA = 4.e-7 * CI;
	//C
	if (ZERO(R)) {
		double CR2 = CR * CR;
		double S = CR2 + (CZ - Z) * (CZ - Z);
		S = sqrt(S) * S;
		BR = 0.;
		BZ = AA * PI / 2. / S * CR2;
		return;
	}
	//C
	double S = CR * CR + R * R + (CZ - Z) * (CZ - Z);
	double P = 2.0 * CR * R;
	double RK2 = 2.0 * P / (S + P);
	double RK = sqrt(RK2);
	double ELPK, ELPE;
	int IER;
	CELIDD(RK2, ELPK, ELPE, IER);
	//C
		if (IER != 0) {
			//		WRITE(6, 1000) IER, R, Z, CR, CZ
			BR = 0.;
			BZ = 0.;
			return;
		}
	//C
	AA = AA / (2.0 * sqrt(S + P));
	BR = AA * (Z - CZ) * (-ELPK + S / (S - P) * ELPE) / R;
	BZ = AA * (ELPK - (S - 2.0 * CR * CR) / (S - P) * ELPE);
	return;
	//1000 FORMAT(1H, '***** ERROR IN CELID(BRING) IER=', I5, ' R=', D15.7,
	//	*' Z=', D15.7, ' CR=', D15.7, ' CZ=', D15.7)
}

//% PSI
//====================================================================== =
double EM_COIL_Loop::psi_(double R, double Z, double CR, double CZ, double CI)
//C---------------------------------------------------------------------- -
{

	if (ZERO(R)) {
		double PSI = 0.e0;
		return PSI;
	}
	//C
	double UM = 4.0E-7 * CI;
	double S = CR * CR + R * R + (CZ - Z) * (CZ - Z);
	double P = CR * R;
	double RK2 = 4.0E0 * P / (S + 2.0 * P);
	double RK = sqrt(RK2);
	double ELPK, ELPE;
	int IER;
	CELIDD(RK2, ELPK, ELPE, IER);
	//C
	if (IER != 0) {
		//        WRITE(6, 1000) IER, R, Z, CR, CZ
		double PSI = 0.;
		return PSI;
	}
	//C---------------------------------------------------------------------- -
	double PSI = UM / RK * sqrt(P) * ((1.0 - RK2 / 2.0) * ELPK - ELPE);
	return PSI;
	//C---------------------------------------------------------------------- -
	//C---------------------------------------------------------------------- -
	//        1000 FORMAT(1H, '***** ERROR IN CELID(PSI) IER=', I5, ' R=', D15.7,
	//            *' Z=', D15.7, ' CR=', D15.7, ' CZ=', D15.7)
}

//C% CELIDD
//C====================================================================== =
void EM_COIL_Loop::CELIDD(double KX, double& ZK, double& ZE, int& ILL) {
	//C------------------------------------------------------------ 1985 / 03 / 26
	//C---------------------------------------------------------------------- -
	double A0 = 1.38629436112;
	double A1 = 0.09666344259;
	double A2 = 0.03590092383;
	double A3 = 0.03742563713;
	double A4 = 0.01451196212;

	double B0 = 0.5;
	double B1 = 0.12498593597;
	double B2 = 0.06880248576;
	double B3 = 0.03328355346;
	double B4 = 0.00441787012;

	double C1 = 0.44325141463;
	double C2 = 0.06260601220;
	double C3 = 0.04757383546;
	double C4 = 0.01736506451;

	double D1 = 0.24998368310;
	double D2 = 0.09200180037;
	double D3 = 0.04069697526;
	double D4 = 0.00526449639;

	//C---------------------------------------------------------------------- -
	if ((KX + 1. > 1.) || (KX + 1. < 2.)) {
		double X1 = 1. - KX;
		double X2 = X1 * X1;
		double X3 = X1 * X2;
		double X4 = X2 * X2;
		double XI = 1. / X1;
		double XL = log(XI);
		double AA = A0 + A1 * X1 + A2 * X2 + A3 * X3 + A4 * X4;
		double BB = B0 + B1 * X1 + B2 * X2 + B3 * X3 + B4 * X4;
		ZK = AA + BB * XL;
		//C---------------------------------------------------------------------- -
		double C0 = 1.;
		double CC = C0 + C1 * X1 + C2 * X2 + C3 * X3 + C4 * X4;
		double DD = D1 * X1 + D2 * X2 + D3 * X3 + D4 * X4;
		ZE = CC + DD * XL;
		ILL = 0;
		return;
		//C---------------------------------------------------------------------- -
	}
	ZK = 0.;
	ZE = 0.;
	ILL = 1;
	return;
}
