#include <assert.h>
#include "basic/EM_defines.h"
#include "source/coil/EM_COIL.h"
#include "source/coil/coil_field.h"
#include "basic/EM_Math.h"


/****    EM_COIL_Block   ***********************************/
EM_COIL_Block::EM_COIL_Block(
	double current,
	const EM_3dPoint& start,
	const EM_3dPoint& end,
	const EM_3dVector& dx,
	const EM_3dVector& dy
) : current(current), start_point(start), end_point(end) ,
	half_width_vector_a(dx), half_width_vector_b(dy){
	Initialize();
}

EM_COIL_Block::EM_COIL_Block( 
		double current,
		const EM_3dPoint &start,
		const EM_3dPoint &end,
		double width,
		double thick
		): current(current), start_point(start), end_point(end){

//	assert(start.z==end.z);

	EM_3dVector dir=end_point-start_point;
	dir.Unitize();
	EM_3dVector  ez(0.,0.,1.);
	half_width_vector_a=EM_CrossProduct(ez,dir)*(width/2.);
	half_width_vector_b=ez*(thick/2.);
}
#ifndef COILBASE
void EM_COIL_Block::Read(EM_TextFile *ffile){
    TITLERD(ffile, "*GCE    * VARIABLE::CURRENT(A)*\n",
             "%s %lf", ND2(namfld, &r1),
             ( 1 ) );
	current  = r1;

    CARDRD(ffile, "*            XS(m)  *   YS(m)  *  ZS(m)  *  XE(m)  *  YE(m)  *  ZE(m)  * \n",
           "%lf %lf %lf %lf %lf %lf", ND6(&r1,&r2, &r3, &r4, &r5, &r6),
           (!(r1==r4 && r2==r5 && r3==r6) ) );
	start_point.Set(r1,r2,r3);
	end_point.Set(r4,r5,r6);

	CARDRD(ffile, "*            W1X(m) *   W1Y(m) *  W1Z(m) *  W2X(m) *  W2Y(m) *  W2Z(m) *\n",
           "%lf %lf %lf %lf %lf %lf", ND6(&r1,&r2, &r3, &r4, &r5, &r6),
           ( r1*r1+r2*r2+r3*r3 !=0. && r4*r4+r5*r5+r6*r6 !=0. ) );

	half_width_vector_a.Set(r1,r2,r3);
	half_width_vector_b.Set(r4,r5,r6);

}

EM_REAL EM_COIL_Block::GetLength() const{
	EM_3dVector d=end_point-start_point;
	return d.Length();
}

EM_REAL EM_COIL_Block::CalcResistance() const{
	EM_REAL s=half_width_vector_a.Length() * half_width_vector_b.Length() *4.;
	EM_REAL j=current;
	EM_REAL l=(end_point-start_point).Length();
//	EM_REAL resistance=l*j*j/(s*conductivity);
	EM_REAL resistance = 0.;
	EM_REAL denominator = s * conductivity;
	if (denominator) {
		resistance = l * j * j / denominator;
	}	return resistance;
}


#endif
EM_3dVector EM_COIL_Block::A_Field( const EM_3dPoint &x) const{
	EM_3dVector A, B;
   field_param param;

   param.calcA=1;
   param.calcB=0;
#ifdef WASEDA_METHOD
   CalcFields_Waseda(x, A, B, &param);
#else
    CalcFields(x, A, B, &param);
#endif
	return A;
}

EM_3dVector EM_COIL_Block::B_Field( const EM_3dPoint &x) const{
	EM_3dVector A, B;
   field_param param;
   param.calcA=0;
   param.calcB=1;

#ifdef WASEDA_METHOD
  CalcFields_Waseda(gx, A, B, &param);
#else
    CalcFields(x, A, B, &param);

#endif
	return B;
} 
#ifndef COILBASE
void EM_COIL_Block::A_B_Field(const EM_3dPoint &x, EM_3dVector &A, EM_3dVector &B) const{

   field_param param;

   param.calcA=1;
   param.calcB=1;

#ifdef WASEDA_METHOD
  CalcFields_Waseda(gx, A, B, &param);
#else
    CalcFields(x, A, B, &param);
#endif
} 



/*	EM_REAL ADJ1=1.e-6;
	EM_REAL ADJ2=1.e-3;*/
#endif
#define ADJ1 0.
#define ADJ2 0.
#define FEPS 1.e-50
#define LOGDIV( A, B)   ( (fabs(A)>FEPS && fabs(B)>FEPS )? log((A)/(B)) : 0. )
#define DIV( A, B)   ( fabs(B)>FEPS ? (A)/(B) : 0. )

void EM_COIL_Block::Initialize(){
	EM_REAL *start=start_point;
	EM_REAL *end=end_point;
	EM_REAL X21, Y21, Z21;


    X21=end[0]-start[0];
    Y21=end[1]-start[1];
    Z21=end[2]-start[2];
    length=sqrt(X21*X21+Y21*Y21+Z21*Z21);
    Exyz[6]=X21/length;
    Exyz[7]=Y21/length;
    Exyz[8]=Z21/length;
	normal_vectors();
    area=4.*S1*S2;
	center[0]=.5*(start[0]+end[0]);
	center[1]=.5*(start[1]+end[1]);
	center[2]=.5*(start[2]+end[2]);
}


void EM_COIL_Block::normal_vectors(){
	EM_REAL *a=half_width_vector_a;
	EM_REAL *b=half_width_vector_b;
	int I;
	static EM_REAL ZOT=1.e-7;
      S1=a[0]*Exyz[6]+a[1]*Exyz[7]+a[2]*Exyz[8];
      Exyz[0]=a[0]-S1*Exyz[6];
      Exyz[1]=a[1]-S1*Exyz[7];
      Exyz[2]=a[2]-S1*Exyz[8];
      S1=sqrt(Exyz[0]*Exyz[0]+Exyz[1]*Exyz[1]+Exyz[2]*Exyz[2]);
      Exyz[0] /= S1;
      Exyz[1] /= S1;
      Exyz[2] /= S1;
      S1 -= minimum(ADJ1,ADJ2* S1);
      Exyz[3]=Exyz[7]*Exyz[2]-Exyz[1]*Exyz[8];
      Exyz[4]=Exyz[0]*Exyz[8]-Exyz[6]*Exyz[2];
      Exyz[5]=Exyz[6]*Exyz[1]-Exyz[0]*Exyz[7];
      S2=fabs(b[0]*Exyz[3]+b[1]*Exyz[4]+b[2]*Exyz[5]);
      S2 -= minimum(ADJ1,ADJ2* S2);
	for(I=0; I<9; ++I)
      if(fabs(Exyz[I]) < ZOT) Exyz[I]=0.;
}

void  EM_COIL_Block::CalcFields(const EM_REAL *X, EM_REAL *A, EM_REAL *B, field_param *param) const{
/*
C
C***********************************************************************
C     CALCULATE GCE FIELD  (GENERAL VARIABLE::CURRENT ELEMENT)
C      INPUT  X(1)  CALCULATION POINT X EM_Coordinate (M)
C             X[1]                    Y
C             X[2]                    Z
C             PARGCE(1)   POINT 1     X  (M)
C             PARGCE(2)   POINT 1     Y  (M)
C             PARGCE[2]   POINT 1     Z  (M)
C             PARGCE[3]   POINT 2     X  (M) FROM POINT 1
C             PARGCE[4]   POINT 2     Y  (M) FROM POINT 1
C             PARGCE[5]   POINT 2     Z  (M) FROM POINT 1
C             PARGCE[6]   POINT 3     X  (M) FROM POINT 1
C             PARGCE[7]   POINT 3     Y  (M) FROM POINT 1
C             PARGCE[8]   POINT 3     Z  (M) FROM POINT 1
C             PARGCE[9]  POINT 4     X  (M)
C             PARGCE(11)  POINT 4     Y  (M)
C             PARGCE(12)  POINT 4     Z  (M)
C             PARGCE(13)  XIGCE          (A) VARIABLE::CURRENT
C
C                           *********************
C                        ***               ***  *
C                     ***      4        ***     *
C                  ***        /]     ***        *
C               ***          /-]- ***           *
C            *********************              *
C            *                 ] *              *
C            *                 ] *              *
C            *                 ] *              *
C            *                 ] *              *
C            *                 ] *              *
C            *                 ] *            ***
C            *                 ] *         ***
C            *                 1 *------ 3
C            *              /    *   ***
C            *           /       ***
C            ********* 2 *********
C
C
C             calcA         =0       NO CALCULATION OF VECTER POTENTIAL A
C                        =1       CALCULATION OF VECTER POTENTIAL A
C             calcB         =0       NO CALCULATION OF MAGNETIC INDUCTION B
C                        =1       CALCULATION OF MAGNETIC INDUCTION B
C
C
C      OUTPUT A(I) I=1,3  X,Y,Z COMPONENT OF VECTER POTENTIAL
C             B(I) I=1,3  X,Y,Z COMPONENT OF MAGNETIC INDUCTION (T)
C
C***********************************************************************
*/

	EM_REAL  Myu0J_by_4Pi, Myu0I_by_4Pi;
	EM_REAL X0G, Y0G, Z0G, X0, Y0, Z0, ZH, ZU, ZUSQ, ZL, ZLSQ, XSYS0, RHO1, RHO2, DMIN;
	EM_REAL RHODOT, RHO12, DD, AT, BI, BXI, BYI, XU, XUSQ, XL, XLSQ, YU, YUSQ, YL, YLSQ, XSYS[4];
	EM_REAL XZ[4], YZ[4], D[8], XPD[8], YPD[8], ZPD[4], PP[4], AXL, AXU, AYL, AYU;
	EM_REAL RX[4], RY[4], RZ[4], AZL, AZU;
	int I;
	static EM_REAL ZOT=1.e-5;
    ZH=.5*length;
    X0G=X[0]-center[0];
    Y0G=X[1]-center[1];
    Z0G=X[2]-center[2];

    X0=Exyz[0]*X0G+Exyz[1]*Y0G+Exyz[2]*Z0G;
    Y0=Exyz[3]*X0G+Exyz[4]*Y0G+Exyz[5]*Z0G;
    Z0=Exyz[6]*X0G+Exyz[7]*Y0G+Exyz[8]*Z0G;

    ZU=Z0-ZH;
    ZUSQ=ZU*ZU;
    ZL=Z0+ZH;
    ZLSQ=ZL*ZL;
    XSYS0=X0*X0+Y0*Y0;
    RHO1=sqrt(XSYS0+ZLSQ);
    RHO2=sqrt(XSYS0+ZUSQ);

	if(ZL > 0.) {
		if(ZU < 0.) DMIN=sqrt(XSYS0);
		else		DMIN=RHO2;
	} else			DMIN=RHO1;

	if(DMIN >= LINE_APP_LIMIT*maximum(S1,S2)) {
		RHODOT=XSYS0+ZU*ZL;
		RHO12=RHO1*RHO2;
        Myu0I_by_4Pi=current*1.e-7;
		if(param->calcA) {
			DD=RHO2/RHO1;
			if(DD > 1.) DD=1./DD;
			if((1.-RHODOT/RHO12) >= ZOT) 
				DD=fabs((RHO1-ZL)/(RHO2-ZU));
			AT=Myu0I_by_4Pi*log(DD);
		}
		if(param->calcB){
			BI=(RHO1+RHO2)/(RHO12*(RHO12+RHODOT));
			BI *= -Myu0I_by_4Pi * length;
			BXI=BI*Y0;
			BYI=BI*X0;
		}
	} else {
		Myu0J_by_4Pi=current*1.e-7/area;
		XU=X0-S1;
		XUSQ=XU*XU;
		XL=X0+S1;
		XLSQ=XL*XL;
		YU=Y0-S2;
		YUSQ=YU*YU;
		YL=Y0+S2;
		YLSQ=YL*YL;
		XSYS[0]=XLSQ+YLSQ;
		XSYS[1]=XUSQ+YLSQ;
		XSYS[2]=XLSQ+YUSQ;
		XSYS[3]=XUSQ+YUSQ;
		XZ[0]=XL*ZL;
		XZ[1]=XU*ZL;
		XZ[2]=XL*ZU;
		XZ[3]=XU*ZU;
		YZ[0]=YL*ZL;
		YZ[1]=YU*ZL;
		YZ[2]=YL*ZU;
		YZ[3]=YU*ZU;

		for( I=0; I<4; ++I){
			D[I]=sqrt(XSYS[I]+ZLSQ);
			D[I+4]=sqrt(XSYS[I]+ZUSQ);
		}

		for(I=0; I<7; I +=2){
			XPD[I]=D[I]+XL;
			XPD[I+1]=D[I+1]+XU;
		}

		for(I=0; I<5; I +=4){ 
			YPD[I ]=D[I ]+YL;
			YPD[I+1]=D[I+1]+YL;
			YPD[I+2]=D[I+2]+YU;
			YPD[I+3]=D[I+3]+YU;
		}

		for(I=0; I<4; ++I)
			ZPD[I]=(D[I+4]+ZU)/(D[I]+ZL);

		PP[0]=DIV(YL*D[0]+XSYS[0], XZ[0]);
		PP[1]=DIV(YU*D[2]+XSYS[2], XZ[0]);
		PP[2]=DIV(YL*D[4]+XSYS[0], XZ[2]);
		PP[3]=DIV(YU*D[6]+XSYS[2], XZ[2]);
		AYL=XL*ATAN4(PP[0],PP[1],PP[2],PP[3]);

		PP[0]=DIV(YL*D[1]+XSYS[1], XZ[1]);
		PP[1]=DIV(YU*D[3]+XSYS[3], XZ[1]);
		PP[2]=DIV(YL*D[5]+XSYS[1], XZ[3]);
		PP[3]=DIV(YU*D[7]+XSYS[3], XZ[3]);
		AYU=XU*ATAN4(PP[0],PP[1],PP[2],PP[3]);

		PP[0]=DIV(XL*D[0]+XSYS[0], YZ[0]);
		PP[1]=DIV(XU*D[1]+XSYS[1], YZ[0]);
		PP[2]=DIV(XL*D[4]+XSYS[0], YZ[2]);
		PP[3]=DIV(XU*D[5]+XSYS[1], YZ[2]);
		AXL=YL*ATAN4(PP[0],PP[1],PP[2],PP[3]);

		PP[0]=DIV(XL*D[2]+XSYS[2], YZ[1]);
		PP[1]=DIV(XU*D[3]+XSYS[3], YZ[1]);
		PP[2]=DIV(XL*D[6]+XSYS[2], YZ[3]);
		PP[3]=DIV(XU*D[7]+XSYS[3], YZ[3]);
		AXU=YU*ATAN4(PP[0],PP[1],PP[2],PP[3]);

		if(param->calcA){
			RY[0]=LOGDIV(YPD[2], YPD[0]);
			RY[1]=LOGDIV(YPD[3], YPD[1]);
			RY[2]=LOGDIV(YPD[6], YPD[4]);
			RY[3]=LOGDIV(YPD[7], YPD[5]);

			for(I=0; I<4; ++I){
				RX[I]=LOGDIV(XPD[2*I+1], XPD[2*I]);

				RZ[I]= ZPD[I] ? log(ZPD[I]): 0.;
			}

			AT=ZU*((YU*RX[3]-YL*RX[2])+(XU*RY[3]-XL*RY[2]))
					-ZL*((YU*RX[1]-YL*RX[0])+(XU*RY[1]-XL*RY[0]))
					+(YU*(XU*RZ[3]-XL*RZ[2])-YL*(XU*RZ[1]-XL*RZ[0]))
					-.5*((YU*AXU-YL*AXL)+(XU*AYU-XL*AYL));

			PP[0]=DIV(ZLSQ+XL*XPD[0], YZ[0]);
			PP[1]=DIV(ZLSQ+XU*XPD[1], YZ[0]);
			PP[2]=DIV(ZLSQ+XL*XPD[2], YZ[1]);
			PP[3]=DIV(ZLSQ+XU*XPD[3], YZ[1]);
			AZL=ZLSQ*ATAN4(PP[0],PP[1],PP[2],PP[3]);

			PP[0]=DIV(ZUSQ+XL*XPD[4], YZ[2]);
			PP[1]=DIV(ZUSQ+XU*XPD[5], YZ[2]);
			PP[2]=DIV(ZUSQ+XL*XPD[6], YZ[3]);
			PP[3]=DIV(ZUSQ+XU*XPD[7], YZ[3]);
			AZU=ZUSQ*ATAN4(PP[0],PP[1],PP[2],PP[3]);

			AT -= .5*(AZU-AZL);
			AT *= Myu0J_by_4Pi;

			if(param->calcB) {
				BXI=(ZL*(RX[1]-RX[0])-ZU*(RX[3]-RX[2]))
					-(XU*(RZ[3]-RZ[1])-XL*(RZ[2]-RZ[0]));
				BYI=(ZL*(RY[1]-RY[0])-ZU*(RY[3]-RY[2]))
					-(YU*(RZ[3]-RZ[2])-YL*(RZ[1]-RZ[0]));

				BXI=Myu0J_by_4Pi*(BXI+(AXU-AXL));
				BYI=Myu0J_by_4Pi*(BYI+(AYU-AYL));
			}
		} else 
		if(param->calcB) {
			RX[0]=LOGDIV((XPD[3]*XPD[0]), (XPD[2]*XPD[1]));
			RX[1]=LOGDIV((XPD[7]*XPD[4]), (XPD[6]*XPD[5]));
			RY[0]=LOGDIV((YPD[3]*YPD[0]), (YPD[2]*YPD[1]));
			RY[1]=LOGDIV((YPD[7]*YPD[4]), (YPD[6]*YPD[5]));
			RZ[0]=LOGDIV(ZPD[2], ZPD[0]);
			RZ[1]=LOGDIV(ZPD[3], ZPD[1]);
			RZ[2]=LOGDIV(ZPD[1], ZPD[0]);
			RZ[3]=LOGDIV(ZPD[3], ZPD[2]);
			BXI=(ZL*RX[0]-ZU*RX[1])-(XU*RZ[1]-XL*RZ[0]);
			BYI=(ZL*RY[0]-ZU*RY[1])-(YU*RZ[3]-YL*RZ[2]);
			BXI=Myu0J_by_4Pi*(BXI+(AXU-AXL));
			BYI=Myu0J_by_4Pi*(BYI+(AYU-AYL));
		}

	}
	if(param->calcB){
		B[0]=Exyz[0]*BXI-Exyz[3]*BYI;
		B[1]=Exyz[1]*BXI-Exyz[4]*BYI;
		B[2]=Exyz[2]*BXI-Exyz[5]*BYI;
	}
	if(param->calcA){
		A[0]=-AT*Exyz[6];
		A[1]=-AT*Exyz[7];
		A[2]=-AT*Exyz[8];
	}
}

#define ZOT 1.e-5
void EM_COIL_Block::CalcFields_Waseda(const EM_REAL *X, EM_REAL *A, EM_REAL *B, struct field_param *param) const{
/*
C
C***********************************************************************
C     CALCULATE GCE FIELD  (GENERAL VARIABLE::CURRENT ELEMENT)
C      INPUT  X(1)  CALCULATION POINT X EM_Coordinate (M)
C             X[1]                    Y
C             X[2]                    Z
C             PARGCE(1)   POINT 1     X  (M)
C             PARGCE(2)   POINT 1     Y  (M)
C             PARGCE[2]   POINT 1     Z  (M)
C             PARGCE[3]   POINT 2     X  (M) FROM POINT 1
C             PARGCE[4]   POINT 2     Y  (M) FROM POINT 1
C             PARGCE[5]   POINT 2     Z  (M) FROM POINT 1
C             PARGCE[6]   POINT 3     X  (M) FROM POINT 1
C             PARGCE[7]   POINT 3     Y  (M) FROM POINT 1
C             PARGCE[8]   POINT 3     Z  (M) FROM POINT 1
C             PARGCE[9]  POINT 4     X  (M)
C             PARGCE(11)  POINT 4     Y  (M)
C             PARGCE(12)  POINT 4     Z  (M)
C             PARGCE(13)  XIGCE          (A) VARIABLE::CURRENT
C
C                           *********************
C                        ***               ***  *
C                     ***      4        ***     *
C                  ***        /]     ***        *
C               ***          /-]- ***           *
C            *********************              *
C            *                 ] *              *
C            *                 ] *              *
C            *                 ] *              *
C            *                 ] *              *
C            *                 ] *              *
C            *                 ] *            ***
C            *                 ] *         ***
C            *                 1 *------ 3
C            *              /    *   ***
C            *           /       ***
C            ********* 2 *********
C
C
C             calcA         =0       NO CALCULATION OF VECTER POTENTIAL A
C                        =1       CALCULATION OF VECTER POTENTIAL A
C             calcB         =0       NO CALCULATION OF MAGNETIC INDUCTION B
C                        =1       CALCULATION OF MAGNETIC INDUCTION B
C
C
C      OUTPUT A(I) I=1,3  X,Y,Z COMPONENT OF VECTER POTENTIAL
C             B(I) I=1,3  X,Y,Z COMPONENT OF MAGNETIC INDUCTION (T)
C
C***********************************************************************
*/
	EM_REAL  Myu0J_by_4Pi, Myu0I_by_4Pi;
	EM_REAL X0G, Y0G, Z0G, X0, Y0, Z0, ZH, ZU, ZUSQ, ZL, ZLSQ, XSYS0, RHO1, RHO2, DMIN;
	EM_REAL RHODOT, RHO12, DD, AT, BI, BXI, BYI, XU,  XL, YU, YL;

	EM_REAL Agsum, x, y, z;
	EM_REAL gxl, gxu, gyl, gyu, gzl, gzu;

    ZH=.5*length;
    X0G=X[0]-center[0];
    Y0G=X[1]-center[1];
    Z0G=X[2]-center[2];

    X0=Exyz[0]*X0G+Exyz[1]*Y0G+Exyz[2]*Z0G;
    Y0=Exyz[3]*X0G+Exyz[4]*Y0G+Exyz[5]*Z0G;
    Z0=Exyz[6]*X0G+Exyz[7]*Y0G+Exyz[8]*Z0G;

    ZU=Z0-ZH;
    ZUSQ=ZU*ZU;
    ZL=Z0+ZH;
    ZLSQ=ZL*ZL;
    XSYS0=X0*X0+Y0*Y0;
    RHO1=sqrt(XSYS0+ZLSQ);
    RHO2=sqrt(XSYS0+ZUSQ);

	if(ZL > 0.) {
		if(ZU < 0.) DMIN=sqrt(XSYS0);
		else		DMIN=RHO2;
	} else			DMIN=RHO1;

	if(DMIN >= LINE_APP_LIMIT*maximum(S1,S2)) {
		RHODOT=XSYS0+ZU*ZL;
		RHO12=RHO1*RHO2;
        Myu0I_by_4Pi=current*1.e-7;
		if(param->calcA) {
			DD=RHO2/RHO1;
			if(DD > 1.) DD=1./DD;
			if((1.-RHODOT/RHO12) >= ZOT) 
				DD=fabs((RHO1-ZL)/(RHO2-ZU));
			AT=Myu0I_by_4Pi*log(DD);
		}
		if(param->calcB){
			BI=(RHO1+RHO2)/(RHO12*(RHO12+RHODOT));
			BI *= -Myu0I_by_4Pi * length;
			BXI=BI*Y0;
			BYI=BI*X0;
		}
	} else {
		Myu0J_by_4Pi=current*1.e-7/area;
		XU=X0-S1;
		XL=X0+S1;
		YU=Y0-S2;
		YL=Y0+S2;

		Agsum=0.;

		x=Z0;
		y=Y0;
		z=-XL;
		gxl=inverse_r_integral_rectangle(ZH, S2, x, y, z);
		Agsum += z*gxl;

		x=Z0;
		y=-Y0;
		z=XU;
		gxu=inverse_r_integral_rectangle(ZH, S2, x, y, z);
		Agsum += z*gxu;

		x=Z0;
		y=-X0;
		z=-YL;
		gyl=inverse_r_integral_rectangle(ZH, S1, x, y, z);
		Agsum += z*gyl;

		x=Z0;
		y=X0;
		z=YU;
		gyu=inverse_r_integral_rectangle(ZH, S1, x, y, z);
		Agsum += z*gyu;

		if(param->calcA) {
			x=X0;
			y=-Y0;
			z=-ZL;
			gzl=inverse_r_integral_rectangle(S1, S2, x, y, z);
			Agsum += z*gzl;

			x=X0;
			y=Y0;
			z=ZU;
			gzu=inverse_r_integral_rectangle(S1, S2, x, y, z);
			Agsum += z*gzu;

			AT=0.5*Myu0J_by_4Pi*Agsum;
		}

		if(param->calcB) {
			BXI=Myu0J_by_4Pi*(gyl-gyu);
			BYI=Myu0J_by_4Pi*(gxl-gxu);
		}

	}
	if(param->calcB){
		B[0]=Exyz[0]*BXI-Exyz[3]*BYI;
		B[1]=Exyz[1]*BXI-Exyz[4]*BYI;
		B[2]=Exyz[2]*BXI-Exyz[5]*BYI;
	}
	if(param->calcA){
		A[0]=-AT*Exyz[6];
		A[1]=-AT*Exyz[7];
		A[2]=-AT*Exyz[8];
	}
}

EM_REAL EM_COIL_Block::inverse_r_integral_rectangle(EM_REAL s1, EM_REAL s2, EM_REAL x, EM_REAL y, EM_REAL z){
	EM_REAL A[2], B[2];
	EM_REAL sum=0.;
	A[0]=s1-x;
	A[1]=-s2-y;
	B[0]=s1-x;
	B[1]=s2-y;
	sum +=inverse_r_integral_line(A, B, z);

	A[0]=B[0];
	A[1]=B[1];
	B[0]=-s1-x;
	B[1]=s2-y;
	sum +=inverse_r_integral_line(A, B, z);


	A[0]=B[0];
	A[1]=B[1];
	B[0]=-s1-x;
	B[1]=-s2-y;
	sum +=inverse_r_integral_line(A, B, z);

	A[0]=B[0];
	A[1]=B[1];
	B[0]=s1-x;
	B[1]=-s2-y;
	sum +=inverse_r_integral_line(A, B, z);

	return sum;
}

EM_REAL EM_COIL_Block::inverse_r_integral_line(EM_REAL A[2], EM_REAL B[2], EM_REAL z) {
	double	ax, ay, AB, dA=0., dB=0., d=0., rA, rB;
	double	ff, log_term, ang, gmt;
/*	static EM_REAL LIMIT=1.e-10;*/

	if(z<0.) z=-z;

	ax = B[0] - A[0];
	ay = B[1] - A[1];
	AB = sqrt(ax*ax + ay*ay);

	if (AB != 0.0) {
		dA = (A[0]*ax + A[1]*ay)/AB;
		dB = (B[0]*ax + B[1]*ay)/AB;
		d = (A[1]*ax - A[0]*ay)/AB;
	}
	else {
#ifndef COILBASE
		fprintf(Stderr, "Edge of element is ZERO\n");
		EXIT(1);
#else
		exit(1);
#endif

	}

	rA = sqrt(A[0]*A[0]+A[1]*A[1]+z*z);
	rB = sqrt(B[0]*B[0]+B[1]*B[1]+z*z);

#if 0

	if (dA > 0.0) {
		log_term = log((rB+dB)/(rA+dA));
	} else {
		ff = fabs(rB-dB);

		if (ff > LIMIT)
			log_term = log((rA-dA)/ff);
		else
			log_term = 0.0;		/* meaningless value */
	}

#else

        if (dA > 0.0) {
                log_term = log((rB+dB)/(rA+dA));
        } else if ((ff = rB-dB) <= 0.0) {
                if (dA*dB <= 0.0)
                        log_term = 0.0;               /* meaningless value */
                else
                        log_term = log(dA/dB);
        } else {
                log_term = log((rA-dA)/ff);
        }

#endif

	if (d != 0.0)
		ang = atan2(d,dA)-atan2(d,dB)-atan2(rA*d,dA*z)+atan2(rB*d,dB*z);
	else
		ang = 0.0;

	gmt = -d*log_term + z*ang;
	return gmt;
}
#ifndef COILBASE
void EM_COIL_Block::CountMeshData(int &no_nodes, int &no_elements) const{
	no_elements +=1;
	no_nodes +=8;
}

void EM_COIL_Block::GetMeshData(MeshElementType &type, EM_3dPoint *coil_nodes, ElementNodes *coil_element_nodes, 
		int &no_coil_nodes, int &no_coil_elements) const {

	no_coil_elements=NumIntegralVolume();
	if(!no_coil_elements) no_coil_elements=1;
    type=hexa;


   EM_3dVector a = half_width_vector_a;
   EM_3dVector b = half_width_vector_b;
   EM_3dPoint x0= start_point;
   EM_3dVector d= (end_point- x0)/ no_coil_elements;

   coil_nodes[0]=x0 -a -b;
   coil_nodes[1]=x0 +a -b;
   coil_nodes[2]=x0 +a +b;
   coil_nodes[3]=x0 -a +b;


   no_coil_nodes=0;
   for(int i=0; i<no_coil_elements; ++i){
	   for(int j=0; j<4; ++j) {
		   coil_element_nodes[i].nodes[j]=no_coil_nodes+j;
		   coil_element_nodes[i].nodes[j+4]=no_coil_nodes+4+j;
		   coil_nodes[no_coil_nodes+4+j]=coil_nodes[no_coil_nodes+j]+d;
	   }
	   no_coil_nodes +=4;
	}
	no_coil_nodes +=4;

}
#endif