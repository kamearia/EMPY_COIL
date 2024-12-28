#ifdef COILBASE
#include "source/coil/EM_COIL.h"
#include "basic/my_math.h"
#else
#include "basic/EM_defines.h"
#endif
#include "source/coil/coil_field.h"
#define FEPS 1.e-50
#define LOGDIV( A, B)   ( (fabs(A)>FEPS && fabs(B)>FEPS )? log((A)/(B)) : 0. )
#define EPSINT 1.e-6
#define KMAX 7
#ifndef COILBASE
extern void reverse_translated_position_by_motion(EM_REAL x[3], EM_REAL xt[3], const MOTION *motion);
#endif
EM_COIL_Arc::EM_COIL_Arc(
		double current,
		double radius,
		double z,
		double radialWidth,
		double axialWidth,
		double start_angle,
		double end_angle
	):current(current), radius(radius)
	 ,axial_width(axialWidth), radial_width(radialWidth)
	 ,Euler_angle_alpha(0.), Euler_angle_beta(0.)
	 ,start_angle(start_angle), end_angle(end_angle)
{
	center = EM_3dPoint(0, 0, z);
}
void EM_COIL_Arc::SetCenter(double x, double y, double z) {
	center = *new EM_3dPoint(x, y, z);
}

#ifndef COILBASE
void EM_COIL_Arc::Read(EM_TextFile *ffile){
    TITLERD(ffile, "*ARC    *  current(A)  *\n",
            "%s %lf", ND2(namfld, &r1),
            ( 1 ) );
	current = r1;

    CARDRD(ffile, "*          X(m)  *  Y(m)  *   Z(m)  * RADIUS(m)* AXIAL_W(m)*RADIAL_W(m)*\n",
           "%lf %lf %lf %lf %lf %lf", ND6(&r1,&r2, &r3, &r4, &r5, &r6),
           (r4>0. && r5>0. && r6>0.) );
	center.Set(r1,r2,r3);
	radius = r4;
	axial_width = r5;
	radial_width = r6;

	CARDRD(ffile, "*          ALPHA(deg)* BETA(deg) * PHI1(deg)   * PHI2(deg) *\n",
           "%lf %lf %lf %lf %d", ND4(&r1,&r2, &r3, &r4),
           (r3< r4 ) );

	Euler_angle_alpha = r1;
	Euler_angle_beta = r2;
	start_angle = r3;
	end_angle = r4;

	Euler_angle_alpha *= PI/180.;
	Euler_angle_beta *= PI/180.;
	start_angle  *= PI/180.;
	end_angle *= PI/180.;
}

EM_REAL EM_COIL_Arc::GetLength() const{
	return radius*(end_angle-start_angle);
}

void EM_COIL_Arc::CountMeshData(int &no_nodes, int &no_elements) const{
	EM_REAL dt = 40.*PI/180.;
	EM_REAL phi1 = start_angle;
	EM_REAL phi2 = end_angle;

	int no_coil_elements;
	do {
		dt /= 2.;
		no_coil_elements =(int) ((phi2-phi1)/dt + 0.9999);
		dt = (phi2-phi1)/ no_coil_elements;
	}
	while(no_coil_elements <= MIN_ARC_DIVISIONS);

	no_elements += no_coil_elements;
	no_nodes += 4*(no_coil_elements+1);
}

EM_REAL EM_COIL_Arc::CalcResistance() const{
	EM_REAL s=axial_width*radial_width;
	EM_REAL j=current;
	EM_REAL l=(end_angle-start_angle)*radius;
//	EM_REAL resistance=l*j*j/(s*conductivity);
	EM_REAL resistance = 0.;
	EM_REAL denominator = s * conductivity;
	if (denominator) {
		resistance = l * j * j / denominator;
	}
	return resistance;
}
#endif
EM_3dVector EM_COIL_Arc::A_Field( const EM_3dPoint &x) const{
	EM_3dVector A, B;
	field_param param;
   
	param.calcA =1;
	param.calcB =0;
	CalcFields (x, A, B, &param);
	return A;
} 

EM_3dVector EM_COIL_Arc::B_Field( const EM_3dPoint &x) const{
	EM_3dVector A, B;
   field_param param; 
   param.calcA =0;
   param.calcB =1;
   CalcFields (x, A, B, &param);
   return B;
}
#ifndef COILBASE
void EM_COIL_Arc::A_B_Field(const EM_3dPoint &x, EM_3dVector &A, EM_3dVector &B) const{
   field_param param;  
   param.calcA =1;
   param.calcB =1;
   CalcFields (x, A, B, &param);
} 

int EM_COIL_Arc::NumPlotElements() const {
	EM_REAL dt = 10. / 180.*PI;
	EM_REAL arcAngle = end_angle - start_angle;
	int numElems = ceil(arcAngle / dt);
	return numElems;
}


extern void rotation_by_Euler_angles(EM_REAL alpha, EM_REAL beta, EM_REAL x0[3], EM_REAL x[3]);
void  EM_COIL_Arc::GetMeshData(MeshElementType &type, EM_3dPoint *coil_nodes, ElementNodes *coil_element_nodes, 
		int &no_coil_nodes, int &no_coil_elements) const {
	type=hexa;
   EM_REAL phi, phi1, phi2, a, s12, s22, alpha, beta, cosp, sinp;
   int i, k, ixyz;
   EM_REAL r, z, x[3], x1[3];
   no_coil_nodes=0;
 
   phi1 = start_angle;
   phi2 = end_angle;
	no_coil_elements=NumIntegralVolume();
	if(!no_coil_elements){
		no_coil_elements=NumPlotElements();
	}
	EM_REAL dt=(phi2-phi1)/no_coil_elements;

   a = radius;
   s12 = axial_width/2.;
   s22 = radial_width/2;
   alpha = Euler_angle_alpha ;
   beta  = Euler_angle_beta ;

   for(i=0; i<= no_coil_elements; ++i)
   {
      phi = phi1 + dt* i ;
      cosp = COS(phi);
      sinp = SIN(phi);

      if(i<no_coil_elements)
         for(k=0; k<8 ; ++ k)
            coil_element_nodes[i].nodes[k] = no_coil_nodes + k;

      for(k=0; k<4; ++k)
	  {
         switch(k)
         {
             case 0:
                r = a-s22;
                z = -s12;
                break;   
             case 1:
                r = a-s22;
                z =  s12;
                break;   
             case 2:
                r = a+s22;
                z =  s12;
                break;   
             case 3:
                r = a+s22;
                z = -s12;
                break;   
         }

         x[0] = r*cosp;
         x[1] = r*sinp;
         x[2] = z;

         rotation_by_Euler_angles(alpha, beta, x, x1);
         for(ixyz=0; ixyz<3; ++ixyz)
            coil_nodes[no_coil_nodes][ixyz] = center[ixyz] + x1[ixyz];

         ++ no_coil_nodes;
      }
   }
}
#endif
void EM_COIL_Arc::CalcFields(const EM_REAL *X, EM_REAL *A, EM_REAL *B, struct field_param *param) const {
/*
C
C***********************************************************************
C     CALCULATE ARC FIELD
C      INPUT  X(1)  CALCULATION POINT X EM_Coordinate (M)
C             X(2)                    Y
C             X(3)                    Z
C             PARARC(1)   CENTER  XARC   (M)
C             PARARC(2)   CENTER  YARC   (M)
C             PARARC(3)   CENTER  ZARC   (M)
C             PARARC(4)   RADIUS  AARC   (M)
C             PARARC(5)           ALPARC (RAD) ROTATION TO Z-AXIS(FIRST)
C             PARARC(6)           BETARC (RAD) ROTATION TO Z-AXIS(SECOND
C             PARARC(7)           APHI1  (RAD) INIATIAL ANGLE
C             PARARC(8)           APHI2  (RAD) END ANGLE
C             PARARC(9)           S1ARC2 (M)   AXIAL WIDTH
C             PARARC(10)          S2ARC2 (M)   RARIAL WIDTH
C             PARARC(11)          XIARC  (A) VARIABLE::CURRENT
C             calcA      =0       NO CALCULATION OF VECTER POTENTIAL A
C                        =1       CALCULATION OF VECTER POTENTIAL A
C             calcB      =0       NO CALCULATION OF MAGNETIC INDUCTION B
C                        =1       CALCULATION OF MAGNETIC INDUCTION B
C
C      OUTPUT A(I) I=1,3  X,Y,Z COMPONENT OF VECTER POTENTIAL
C             B(I) I=1,3  X,Y,Z COMPONENT OF MAGNETIC INDUCTION (T)
C
C**********************************************************************
C

*/
	EM_REAL XARC, YARC, ZARC, AARC, ALPARC, BETARC, APHI1, APHI2, S1ARC2, S2ARC2, XIARC;
	EM_REAL a, b;
	EM_REAL XZERO, YZERO, ZZERO, X1, Y1, zc;
	EM_REAL SINEA, COSINA, SINEB, COSINB, COSA, SINA;
	EM_REAL TWAVE, C1, C2, AMRSQ, TUR, DMIN, PHI1, PHI2;
	EM_REAL RL, RU, ART, BRT, ALXU, ALXL, XI, BX1, BY1, T, AX1, AY1, AYC;
	int AB;
	int i, J;
	static EM_REAL ZOT=1.e-10;
//	static EM_REAL PHI0=0.001; 
	static EM_REAL PHI0=0.; 
	field_sum fsum;
	EM_REAL AU,AUSQ,AL,ALSQ,ZZU,ZUSQ,ZZL,ZLSQ,ARZ[4],rp,
		RHOSQ,rc,ASQ,TAR,ZSQ,ARZSQ,AREA;

    XARC  =center[0];
    YARC  =center[1];
    ZARC  =center[2];
    AARC  =radius;
    ALPARC=Euler_angle_alpha;
    BETARC= Euler_angle_beta;
    APHI1 =start_angle;
    APHI2 =end_angle;
    S1ARC2=axial_width;
    S2ARC2=radial_width;
/*    S1ARC2=S1ARC2-min(ADJ1,ADJ2*S1ARC2);
    S2ARC2=S2ARC2-min(ADJ1,ADJ2*S2ARC2); */
    XIARC =current;
    b=S1ARC2/2.;
    a=S2ARC2/2.;
    XIARC = XIARC/(S1ARC2*S2ARC2);
	XI=XIARC*Myu0by4PI;

    fsum.ATHWAV=0.;
    fsum.ARHO=0.;
    fsum.BTHETA=0.;
    fsum.BRHO=0.;
    fsum.BZWAVE=0.;
    XZERO=X[0]-XARC;
    YZERO=X[1]-YARC;
    ZZERO=X[2]-ZARC;						/*  zp-zc */
    if((ALPARC == 0.) && (BETARC == 0.)) {
		X1=XZERO;
		Y1=YZERO;
		zc=ZZERO;
		AB=1;
	} else {
		SINEA=sin(ALPARC);
		COSINA=cos(ALPARC);
		SINEB=sin(BETARC);
		COSINB=cos(BETARC);
/*
		X1=XZERO*COSINA+YZERO*SINEA;
		EM_REAL Yt=-XZERO*SINEA+YZERO*COSINA;
		Y1=Yt*COSB + ZZERO*SINEB;
		zc =-Yt*SINEB+ZZERO*COSINB;
*/
		X1=XZERO*COSINA+YZERO*SINEA;
		zc=XZERO*SINEA-YZERO*COSINA;
		Y1=ZZERO*SINEB-zc*COSINB;
		zc =zc*SINEB+ZZERO*COSINB;

		AB=0;
	}

	RHOSQ=X1*X1+Y1*Y1;
    ZSQ=zc*zc;
    if(ZOT*ZSQ < RHOSQ) {
		TWAVE=atan2(Y1,X1);
		rp=sqrt(RHOSQ);
		C1=X1/rp;
		C2=Y1/rp;
	} else {
      TWAVE=0.;
      C1=1.;
      C2=0.;
      RHOSQ=0.;
      rp=RHOSQ;
	}

    rc=AARC;
    ASQ=rc*rc;
    AMRSQ=(rc-rp)*(rc-rp);
    TUR=rp+rp;					/* 2*rp   */
    TAR=rc*TUR;
    ARZSQ=AMRSQ+ZSQ;
    PHI1=APHI1-TWAVE;
    PHI2=APHI2-TWAVE;
	
	if( PHI2 < -PI+PHI0) {
		PHI1=PHI1+TwoPi;
		PHI2=PHI2+TwoPi;
	} else {
		if( PHI1 > PI-PHI0) {
			PHI1=PHI1-TwoPi;
			PHI2=PHI2-TwoPi;
		}
	}
//	if(fabs(PHI1) < PHI0) PHI1=0.;
//    if(fabs(PHI2) < PHI0) PHI2=0.;

	DMIN=ARZSQ;
    if(PHI1*PHI2 > 0.)
		DMIN +=TAR*(1.-maximum(cos(PHI1),cos(PHI2)));

	AREA=4.*a*b;
    param->DFT= LINE_APP_LIMIT*maximum(a, b);
	param->DFT = param->DFT*param->DFT;

	param->RHO=rp;
	param->RHOSQ=RHOSQ;
	param->ZSQ=ZSQ;
	param->ASQ=ASQ;
	param->AA=rc;
	param->ARZSQ=ARZSQ;
	param->TAR=TAR;
	param->AREA=AREA;
	param->Z1=zc;

    if( DMIN >= param->DFT/* && RHO*/) {
 		param->symmetricArc=0;

		field_ARC_line(PHI1, PHI2, param, &fsum );

	} else {

		AU=rc+a;
		AUSQ=(AU-rp)*(AU-rp);		/* (rc+a-rp)**2  */
		AL=rc-a;
		ALSQ=(AL-rp)*(AL-rp);		/* (rc-a-rp)**2  */
		ZZU=zc-b;
		ZUSQ=ZZU*ZZU;				/* (zp-zc-b)**2     */
		ZZL=zc+b;
		ZLSQ=ZZL*ZZL;				/* (zp-zc+b)**2     */
		ARZ[0]=ALSQ+ZLSQ;			/* (rc-a-rp)**2+ (zp-zc+b)**2  */
		ARZ[1]=ALSQ+ZUSQ;			/* (rc-a-rp)**2+ (zp-zc-b)**2  */
		ARZ[2]=AUSQ+ZLSQ;			/* (rc+a-rp)**2+ (zp-zc+b)**2  */
		ARZ[3]=AUSQ+ZUSQ;			/* (rc+a-rp)**2+ (zp-zc-b)**2  */


		if(rp ) {
			param->AL=AL;
			param->AU=AU;
			param->ALSQ=ALSQ;
			param->AUSQ=AUSQ;
			param->ZLSQ=ZLSQ;
			param->ZUSQ=ZUSQ;
			param->ZZL=ZZL;
			param->ZZU=ZZU;
			for(i=0; i<4; ++i) param->ARZ[i]=ARZ[i];
			param->ARU=TUR*AU;				/*  2*rp*(rc+a)   */
			param->ARL=TUR*AL;				/*  2*rp*(rc-a)   */

			if(PHI1*PHI2>=0.) {
				param->symmetricArc=0;
				segment_integral(PHI1, PHI2, param, &fsum );
			}
			else{
/*				if(fabs(PHI1+PHI2)<PHI0 ){
					param->symmetricArc=1;
					segment_integral(0., PHI2, param, &fsum );
				} else*/ {
					if( PHI1< -PHI2){
						param->symmetricArc=0;
						segment_integral(PHI1,-PHI2, param, &fsum );
						param->symmetricArc=1;
						segment_integral(0.,PHI2, param, &fsum );
					} else {
						param->symmetricArc=1;
						segment_integral(PHI1,0., param, &fsum );
						param->symmetricArc=0;
						segment_integral(-PHI1, PHI2, param, &fsum );
					}
				}
			}


		} else {

			for(J=0; J<4; ++J) ARZ[J]=sqrt(ARZ[J]);
			RL=ARZ[2]-ARZ[0];
			RU=ARZ[3]-ARZ[1];
			SINA=sin(PHI2)-sin(PHI1);
			COSA=cos(PHI1)-cos(PHI2);

			if(param->calcA) {
				ART=.5*((ZZU*RU-ZZL*RL)+AUSQ*log((ARZ[3]+ZZU)/(ARZ[2]+ZZL))
					-ALSQ*log((ARZ[1]+ZZU)/(ARZ[0]+ZZL)));
				fsum.ARHO=COSA*ART;
				fsum.ATHWAV=-SINA*ART;
			}

			if(param->calcB) {
				BRT=RU-RL;
				fsum.BRHO=SINA*BRT;
				fsum.BTHETA=COSA*BRT;
				ALXU=log((AU+ARZ[3])/(AL+ARZ[1]));
				ALXL=log((AU+ARZ[2])/(AL+ARZ[0]));
				fsum.BZWAVE=(PHI2-PHI1)*(ZZL*ALXL-ZZU*ALXU);
			}
		}
	}

	if(param->calcB) {

		BX1=XI*(fsum.BRHO*C1-fsum.BTHETA*C2);
		BY1=XI*(fsum.BRHO*C2+fsum.BTHETA*C1);
		fsum.BZWAVE *=XI;

		if(AB==0){
			T=BY1*COSINB-fsum.BZWAVE*SINEB;
			B[0]=BX1*COSINA-T*SINEA;
			B[1]=BX1*SINEA+T*COSINA;
			B[2]=BY1*SINEB+fsum.BZWAVE*COSINB;
		} else {
			B[0]=BX1;
			B[1]=BY1;
			B[2]=fsum.BZWAVE;
		}

/*		if(!calcA) return; */
	}

	if(param->calcA) {
		AX1=XI*(fsum.ARHO*C1-fsum.ATHWAV*C2);
		AY1=XI*(fsum.ARHO*C2+fsum.ATHWAV*C1);

		if(AB==0) {
			AYC=AY1*COSINB;
			A[0]=AX1*COSINA-AYC*SINEA;
			A[1]=AX1*SINEA+AYC*COSINA;
			A[2]=AY1*SINEB;
		} else {
			A[0]=AX1;
			A[1]=AY1;
			A[2]=0.;
		}
	}
}

void EM_COIL_Arc::segment_integral(EM_REAL PHI1, EM_REAL PHI2, field_param *param, field_sum *sum){
	EM_REAL DMIN, dphi, PP;
	EM_REAL phiLower, phiUpper;
	int div, n;
	static EM_REAL DPMN=0.0001;

	dphi=PHI2-PHI1;
	if(fabs(dphi)<DPMN) return;
	div=(int) (dphi/HalfPi);
	if(fabs(dphi-div*HalfPi) > DPMN) ++div;
	if(div==0) div=1;
	dphi /=div;
	PP=PHI1;
	for(n=1; n<=div; ++n) {
		phiLower=PP;
		PP=PP+dphi;
		phiUpper=PP;

		if((phiLower-TwoPi)*(phiUpper-TwoPi)<0.){
			phiLower -=TwoPi;
			phiUpper -=TwoPi;
		}
		if(phiLower*phiUpper  > 0.) {
			DMIN=param->ARZSQ+param->TAR*(1.-maximum(cos(phiLower),cos(phiUpper)));
			if(DMIN >= param->DFT) 
				field_ARC_line(phiLower,phiUpper, param, sum);
			else{
				if(phiLower>0.){
					if(phiLower<TwoDegree){
						arc_integration(phiLower, TwoDegree, param, sum);
						phiLower=TwoDegree;
					}
					arc_integration(phiLower,phiUpper, param, sum);
				} else if(phiUpper<0.){
					if(phiUpper>-TwoDegree){
						arc_integration(phiLower, -TwoDegree, param, sum);
						phiLower=-TwoDegree;
					}
					arc_integration(phiLower,phiUpper, param, sum);
				}
			}
		} else {
			if(phiLower<0.){
				if(phiLower<-TwoDegree) {
					arc_integration(phiLower, -TwoDegree, param, sum);
					phiLower=-TwoDegree;
				}
				arc_integration(phiLower, 0., param, sum);
			}
			if(phiUpper>0){
				if(phiUpper>TwoDegree) {
					arc_integration(TwoDegree, phiUpper, param, sum);
					phiUpper=TwoDegree;
				}
				arc_integration(0., phiUpper, param, sum);
			}
		}
	}

}
#ifndef COILBASE
void EM_COIL_Arc::Integrand::Initialize(){
	phi0=parent->start_angle;
	dphi=(parent->end_angle-phi0)/parent->NumIntegralVolume();
	EM_REAL alpha=parent->Euler_angle_alpha;
	EM_REAL beta=parent->Euler_angle_beta;
	cosa=cos(alpha);
	sina=sin(alpha);
	cosb=cos(beta);
	sinb=sin(beta);
}

void EM_COIL_Arc::Integrand::Next(){
	phi0 +=dphi;
}




EM_REAL EM_COIL_Arc::CalcInductance(const EM_COIL *source, const MOTION *targetMotion) const {

	double xmin[3],xmax[3];
	int nx=integParam->intx;
	int ny=integParam->inty;
	int nz=integParam->intz;

	double *intAJ;

	EM_Interval<double> range[3];
	range[0].Set(-0.5, 0.5);
	range[1].Set(-0.5, 0.5);
	range[2].Set( 0., 1.);


	JdotA ja(this, source, targetMotion);
	EM_Integrand *integrand=&ja;
	EM_Integrator *integrator;
	if(nz)
		integrator=new EM_Gauss3dIntegrator(integrand, nx, ny, nz);
	else
		integrator=new EM_Gauss2dAdaptive1dIntegrator(integrand, nx, ny, KMAX, EPSINT);

	int ndiv=NumIntegralVolume();

	EM_REAL inductance=0.;
	for(int i=0; i<ndiv; ++i){
		intAJ=integrator->Calc(range);
			inductance +=intAJ[0];
			ja.Next();
	}
	inductance *=  ja.Factor()*this->current;
	
	return inductance;
}

double * EM_COIL_Arc::JdotA::Func(const double *x){
	EM_3dPoint xl, xg;
	double X0,Y0,Z0,X1,Y1, Y2,Z2;
	double AY1,AZ1, AX2,AY2;
	double r, t;
	r=parent->radius+x[0]*parent->radial_width;
	Z2=x[1]*parent->axial_width;
	t=phi0+x[2]*dphi;

	X1=r*cos(t);
	Y2=r*sin(t);
/*	Z2=z;  */

/*	X1=X2;  */
	Y1=Y2*cosb - Z2*sinb;
	Z0=Y2*sinb + Z2*cosb;

	X0=X1*cosa -Y1*sina;
	Y0=X1*sina +Y1*cosa;
/*	Z0=Z1;  */
	
	xl=EM_3dPoint(X0,Y0,Z0) +parent->center;
	const EM_Coordinate *coord=parent->GetCoordinate();
	xg=coord->ToGlobal(xl);
	EM_3dVector ag;
	ag=source->A_Field(xg, this->targetMotion);
	EM_3dVector al;
	EM_Coordinate::Helper helper(coord, xg);
	al=coord->ToLocal(ag, helper);

	AX2=al[0]*cosa + al[1]*sina; 
	AY1=-al[0]*sina+ al[1]*cosa;
	AZ1=al[2];

/*	AX2=AX1;  */
	AY2=AY1*cosb+AZ1*sinb;
/*	AZ2=-AY1*sinb+AZ1*cosb;  */

/*	aj[0]=(-AX2*sin(t)+AY2*cos(t))*r;  */
	EM_REAL aj=-AX2*Y2+AY2*X1;
#ifdef STATISTIC
	++no_A_evaluation_points;
#endif
	value[0]=aj;
	return value;
}


double * EM_COIL_Arc::JcrossB::Func(const double *x){
	EM_3dVector force;
	double X0,Y0,Z0,X1,Y1, Y2,Z2;
	double JX1, JY1, JZ1, JY2;
	double r, t;
	r=parent->radius+x[0]*parent->radial_width;
	Z2=x[1]*parent->axial_width;
	t=phi0+x[2]*dphi;

	X1=r*cos(t);
	Y2=r*sin(t);

	Y1=Y2*cosb - Z2*sinb;
	Z0=Y2*sinb + Z2*cosb;

	X0=X1*cosa -Y1*sina;
	Y0=X1*sina +Y1*cosa;
	
	EM_3dPoint xl=EM_3dPoint(X0,Y0,Z0) +parent->center;

	const EM_Coordinate *coord=parent->GetCoordinate();
	EM_3dPoint xg;
	xg=coord->ToGlobal(xl);
	EM_3dPoint x_tmp;
	reverse_translated_position_by_motion(xg, x_tmp, motion);
	EM_3dVector B=B_Field(x_tmp);

	JX1=-Y2;
	JY2= X1;
/*	JZ2=0.;  */
	 
/*	JX1=JX2;  */
	JY1=JY2*cosb/* - JZ2*sinb */;
	JZ1=JY2*sinb /*+ JZ2*cosb */;

	EM_3dVector jl;
	jl[0]=JX1*cosa -JY1*sina;
	jl[1]=JX1*sina +JY1*cosa;
	jl[2]=JZ1;

	EM_3dVector jg;
	EM_Coordinate::Helper helper(coord, x_tmp);
	coord->ToGlobal(jl, helper).Conv(jg);
	translate_vector_by_motion(jg, motion);
    force=EM_CrossProduct(jg, B);

#ifdef STATISTIC
	++no_B_evaluation_points;
#endif
	force.Conv(value);
	return value;
}


EM_3dVector * EM_COIL_Arc::CalcForce(EM_REAL current, const MOTION *motion, const EM_MeshSet *mesh, EM_3dVector *force) const{
	double xmin[3],xmax[3];
	int nx=integParam->intx;
	int ny=integParam->inty;
	int nz=integParam->intz;

	EM_Interval<double> range[3];
	range[0].Set(-0.5, 0.5);
	range[1].Set(-0.5, 0.5);
	range[2].Set( 0., 1.);

	int ndiv=integParam->integral_volume_division;
	current *=this->current;
	EM_COIL_Arc::JcrossB jb(this, motion, mesh/*, B_integration_by_v_elements_at_a_point*/);
	EM_Integrand *integrand=&jb;
	EM_Integrator *integrator;
	if(nz)
		integrator=new EM_Gauss3dIntegrator(integrand, nx, ny, nz);
	else
		integrator=new EM_Gauss2dAdaptive1dIntegrator(integrand, nx, ny, KMAX, EPSINT);

	for(int i=0; i<ndiv; ++i, ++force){
		EM_3dVector f=integrator->Calc(range);

		*force += f*(jb.Factor()*current);
		jb.Next();
	}
	return force;
}

#endif
extern int adaptive_1D_integration(int ncal, void func(double, double *, void *), double a, double b, double eps, 
		 double *result, double *work, void *param);

void EM_COIL_Arc::arc_field_integrand(EM_REAL ALPHA, EM_REAL *F, void *param0){
	EM_REAL sinPhi, cosPhi, RSA, rpCosPhi, GXL, GXU, COSA1, ARCL,ARCU, CLSQ, CUSQ;
	EM_REAL RL, RU, AXU, AXL, ALXU, ALXL, ZAZA, RZU, RZL, RLZ, RUZ, BRT, F3, ZLRS, ZURS, ART;
	EM_REAL R[4], PP[4];
	int index;
/*	EM_REAL BUNBO;*/
	field_param *param=(field_param *) param0; 
	EM_REAL AU,AUSQ,AL,ALSQ,ARU,ARL,ZZU,ZUSQ,ZZL,ZLSQ,rp;

	rp=param->RHO;
	AL=param->AL;					/*  rc-a */
	AU=param->AU;					/*  rc+a */
	ARL=param->ARL;					/*  2*rp*(rc-a)   */
	ARU=param->ARU;					/*  2*rp*(rc+a)   */
	ALSQ=param->ALSQ;				/* (rc-a-rp)**2  */
	AUSQ=param->AUSQ;				/* (rc+a-rp)**2  */
	ZLSQ=param->ZLSQ;				/* (zp-zc+b)**2     */
	ZUSQ=param->ZUSQ;				/* (zp-zc-b)**2     */
	ZZL=param->ZZL;					/*  zp-zc+b         */
	ZZU=param->ZZU;					/*  zp-zc-b         */

	sinPhi=sin(ALPHA);
    RSA=rp*sinPhi;
    cosPhi=cos(ALPHA);
    rpCosPhi=rp*cosPhi;
    GXL=AL-rpCosPhi;				/*  rc-a-rp*cos(t)  */
    GXU=AU-rpCosPhi;				/*  rc+a-rp*cos(t)  */
    COSA1=1.-cosPhi;
    if(COSA1 == 0.0) COSA1=ALPHA*ALPHA/2.;
    ARCL=ARL*COSA1;					/*  2*rp*(rc-a)*(1-cos(t))   */
    ARCU=ARU*COSA1;					/*  2*rp*(rc+a)*(1-cos(t))   */
    CLSQ=ARCL+ALSQ;					/*  2*rp*(rc-a)*(1-cos(t)) + (rc-a-rp)**2 = r**2 + rp**2 -2*r*rp*cos(t) /.r->(rc-a)   */
    CUSQ=ARCU+AUSQ;					/*  2*rp*(rc+a)*(1-cos(t)) + (rc+a-rp)**2 = r**2 + rp**2 -2*r*rp*cos(t) /.r->(rc+a)   */
		/*  R[r,z,t,rp,zp] = sqrt( r**2+rp**2-2*r*rp*cos(t)+(zp-z)**2 ) */ 
	R[0]=sqrt(param->ARZ[0]+ARCL);  /*  R[rc-a,zc-b,t,rp,zp] =sqrt( (r-rp)**2+ (zp-z)**2 + 2*rp*r*(1-cos(t)) ) /. {r->rc-a, z->zc-b} */
    R[2]=sqrt(param->ARZ[2]+ARCU);  /*  R[rc+a,zc-b,t,rp,zp] =sqrt( (r-rp)**2+ (zp-z)**2 + 2*rp*r*(1-cos(t)) ) /. {r->rc+a, z->zc-b} */
    RL=R[2]-R[0];
    R[1]=sqrt(param->ARZ[1]+ARCL);  /*  R[rc-a,zc+b,t,rp,zp] =sqrt( (r-rp)**2+ (zp-z)**2 + 2*rp*r*(1-cos(t)) ) /. {r->rc-a, z->zc+b} */
    R[3]=sqrt(param->ARZ[3]+ARCU);  /*  R[rc+a,zc+b,t,rp,zp] =sqrt( (r-rp)**2+ (zp-z)**2 + 2*rp*r*(1-cos(t)) ) /. {r->rc+a, z->zc+b} */
    RU=R[3]-R[1];
    AXU=((GXU*R[3]+CUSQ)+ZUSQ)/R[3]; /* rc+a-rp*cos(t) + R[rc+a,zc+b,t,rp,zp] */
    AXL=((GXL*R[1]+CLSQ)+ZUSQ)/R[1]; /* rc-a-rp*cos(t) + R[rc-a,zc+b,t,rp,zp] */
    ALXU=log(AXU/AXL);
    AXU=((GXU*R[2]+CUSQ)+ZLSQ)/R[2]; /* rc+a-rp*cos(t) + R[rc+a,zc-b,t,rp,zp] */
    AXL=((GXL*R[0]+CLSQ)+ZLSQ)/R[0]; /* rc-a-rp*cos(t) + R[rc-a,zc-b,t,rp,zp] */
    ALXL=log(AXU/AXL);
    ZAZA=ZZU*ALXU-ZZL*ALXL;
    RZU=((ZUSQ+ZZU*R[1])+CLSQ)/R[1]; /* (zp-zc-b) + R[rc-a,zc+b,t,rp,zp]  */
    RZL=((ZLSQ+ZZL*R[0])+CLSQ)/R[0]; /* (zp-zc+b) + R[rc-a,zc-b,t,rp,zp]  */
    RLZ=fabs(RZU/RZL);
    RZU=((ZUSQ+ZZU*R[3])+CUSQ)/R[3]; /* (zp-zc-b) + R[rc+a,zc+b,t,rp,zp]  */
    RZL=((ZLSQ+ZZL*R[2])+CUSQ)/R[2]; /* (zp-zc+b) + R[rc+a,zc-b,t,rp,zp]  */
    RUZ=fabs(RZU/RZL);

	index=0;
    if(param->calcB) {
		PP[0]=(ZZL*GXL)/(RSA*R[0]);  /* (zp-zc+b)*(rc-a-rp*cos(t))/( rp*sin(t)*R[rc-a,zc-b,t,rp,zp]) */
		PP[1]=(ZZU*GXL)/(RSA*R[1]);  /* (zp-zc-b)*(rc-a-rp*cos(t))/( rp*sin(t)*R[rc-a,zc+b,t,rp,zp]) */
		PP[2]=(ZZL*GXU)/(RSA*R[2]);  /* (zp-zc+b)*(rc+a-rp*cos(t))/( rp*sin(t)*R[rc+a,zc-b,t,rp,zp]) */
		PP[3]=(ZZU*GXU)/(RSA*R[3]);  /* (zp-zc-b)*(rc+a-rp*cos(t))/( rp*sin(t)*R[rc+a,zc+b,t,rp,zp]) */
		BRT=(RU-RL)+rpCosPhi*(ALXU-ALXL);

		F[index++]=cosPhi*BRT;
/*		F3=RSA*ATAN4(PP[0],PP[1],PP[2],PP[3])-rpCosPhi*log(RLZ/RUZ)-ZAZA; */
/*  TODO 2008.10.15  */
		F3=RSA*ATAN4(PP[0],PP[1],PP[2],PP[3])-rpCosPhi*LOGDIV(RLZ,RUZ)-ZAZA;
		if(param->symmetricArc==0) {
			F[index++]=sinPhi*BRT;
			F[index++]=F3;
		} else {
			F[index++]=F3;
		}
	}

	if(param->calcA) {

		ZLRS=ZZL*RSA;
		ZURS=ZZU*RSA;
		if(ZLRS){
			PP[0]=(CLSQ+GXL*R[0])/ZLRS;
			PP[2]=(CUSQ+GXU*R[2])/ZLRS;
		} else {
			PP[0]=0.;
			PP[2]=0.;
		}
		if(ZURS){
			PP[1]=(CLSQ+GXL*R[1])/ZURS;
			PP[3]=(CUSQ+GXU*R[3])/ZURS;
		} else {
			PP[1]=0.;
			PP[3]=0.;
		}
/*
Miyata equations
		BUNBO=R[0]*RSA;
		if(BUNBO) PP[0]=(GXL*ZZL)/BUNBO;
		else PP[0]=0.;

		BUNBO=R[2]*RSA;
		if(BUNBO) PP[2]=(GXU*ZZL)/BUNBO;
		else PP[2]=0.;

		BUNBO=R[1]*RSA;
		if(BUNBO) PP[1]=(GXL*ZZU)/BUNBO;
		else PP[1]=0.;
		BUNBO=R[3]*RSA;
		if(BUNBO) PP[3]=(GXU*ZZU)/BUNBO;
		else PP[3]=0.;
*/
		ART=.5*(ZZU*RU-ZZL*RL)
			+rpCosPhi*(ZAZA-RSA*ATAN4(PP[0],PP[1],PP[2],PP[3]))
			+((rpCosPhi*GXU+.5*CUSQ)*log(RUZ)-(rpCosPhi*GXL+.5*CLSQ)*log(RLZ));

		if(param->symmetricArc == 0) {
			F[index++]=sinPhi*ART;
			F[index++]=cosPhi*ART;
		} else
			F[index++]=cosPhi*ART;
	}
}

void EM_COIL_Arc::arc_integration(EM_REAL P1,EM_REAL P2, field_param *param, field_sum *sum) {
    EM_REAL work[655/*=131*5*/], result[5];
	int index, k;

	index=0;
    if(param->calcB) {
		++index;
		if(param->symmetricArc==0) {
			index +=2;
		} else {
			index++;
		}

	} 
	if(param->calcA) {
		if(param->symmetricArc == 0) {
			index+=2;
		} else{
			index++;
		}
	}

	k=adaptive_1D_integration(index, arc_field_integrand, P1, P2, ArcIntegralAccuracy, result,  work, (void *)param);

	index=0;
    if(param->calcB) {
		if(param->symmetricArc==0) {
			sum->BRHO+=result[index++];
			sum->BTHETA+=result[index++];
			sum->BZWAVE+=result[index++];
		} else {
			sum->BRHO+=result[index++]*2.;
			sum->BZWAVE+=result[index++]*2.;

		}
	}
	if(param->calcA) {
		if(param->symmetricArc == 0) {
			sum->ARHO +=result[index++];
			sum->ATHWAV -=result[index++];
		} else{
			sum->ATHWAV-=result[index++]*2.;
		}
	}

}

//extern void  EllepticIntegral(EM_REAL k2, EM_REAL kd2, EM_REAL *K, EM_REAL *E, EM_REAL omega1, EM_REAL omega2);

extern void  EllepticIntegral2(EM_REAL k2, EM_REAL kd2, EM_REAL *K, EM_REAL *E, EM_REAL omega1, EM_REAL omega2);
extern EM_REAL EPS_k2;
void EM_COIL_Arc::field_ARC_line(EM_REAL PHI1, EM_REAL PHI2, field_param *param, field_sum *sum) {
	EM_REAL RZSQ, ASRZS, APR, APRSQ, RAZSQ,   RAZ;
	EM_REAL SNOMG1, CSOMG1, SNCS1, SN1SQ, SNOMG2, CSOMG2, SNCS2, SN2SQ;
	EM_REAL C, D1, D2, SL, D, ARO, G;
	EM_REAL CASQ,CPSQ,OMEGA1,OMEGA2,FEI,EEI;
	EM_REAL RHO,RHOSQ,AA,ASQ,TAR,Z1,ZSQ,ARZSQ,AREA;
	static EM_REAL ZOT=1.e-10;
	EM_REAL fac=1.;
	if(param->symmetricArc==1) fac=2.;

	RHO=param->RHO;
	RHOSQ=param->RHOSQ;
	ZSQ=param->ZSQ;
	ASQ=param->ASQ;
	AA=param->AA;
	ARZSQ=param->ARZSQ;
	TAR=param->TAR;
	AREA=param->AREA;
	Z1=param->Z1;

      RZSQ=RHOSQ+ZSQ;
      ASRZS=ASQ+RZSQ;
      APR=AA+RHO;
      APRSQ=APR*APR;
      RAZSQ=APRSQ+ZSQ;
      CASQ=(TAR+TAR)/RAZSQ;
      CPSQ=ARZSQ/RAZSQ;
      RAZ=sqrt(RAZSQ);
      OMEGA1=.5*(PHI1-PI);
      SNOMG1=sin(OMEGA1);
      CSOMG1=cos(OMEGA1);
      SNCS1=SNOMG1*CSOMG1;
      SN1SQ=SNOMG1*SNOMG1;
      OMEGA2=.5*(PHI2-PI);
      SNOMG2=sin(OMEGA2);
      CSOMG2=cos(OMEGA2);
      SNCS2=SNOMG2*CSOMG2;
      SN2SQ=SNOMG2*SNOMG2;

	EllepticIntegral2(CASQ, CPSQ, &FEI, &EEI, OMEGA1, OMEGA2);

	if(CASQ >= EPS_k2) {
		D1=sqrt(1.-CASQ*SN1SQ);
		D2=sqrt(1.-CASQ*SN2SQ);

		if(param->calcA) {
			C=(RAZ/RHO)*AREA*fac;
			if(param->symmetricArc==0)
				sum->ARHO += C*(D1-D2);
			sum->ATHWAV += C*(FEI*(ASRZS/RAZSQ)-EEI);
		}
		if(param->calcB){
			D1=1./D1;
			D2=1./D2;
			SL=EEI-CASQ*(SNCS2*D2-SNCS1*D1);
			C=AREA/RAZ*fac;
			if(Z1 !=0.) {
				D=(Z1/RHO)*C;
				if(param->symmetricArc==0)
					sum->BTHETA += D*(D1-D2);
				sum->BRHO += D*(SL*(ASRZS/ARZSQ)-FEI);
			}

			if(fabs(1.-RHOSQ/ASQ) > ZOT) SL*=((RZSQ-ASQ)/ARZSQ);
			sum->BZWAVE += C*(FEI-SL);
		}

	} else {
		C=((AA+AA)/RAZ)*AREA*fac;
		ARO=C*(SN2SQ-SN1SQ);
		if(param->calcA) {
			if (param->symmetricArc==0) sum->ARHO +=ARO;
/* TODO 2006 */
			sum->ATHWAV += C*(SNCS1-SNCS2);
/* 1999/7/19 A. Kameari
c		sum->ATHWAV-=C*FEI   */
		}
		if(param->calcB) {
			sum->BZWAVE += (C/RAZSQ)*APR*EEI;
			G=Z1/RAZSQ;
			if(param->symmetricArc==0) sum->BTHETA -= G*ARO;
			sum->BRHO += G*C*(SNCS1-SNCS2);
		}
	}
}
