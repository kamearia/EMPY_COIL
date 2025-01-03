#include "basic/EM_defines.h"
#include "basic/EM_Math.h"
#include <cmath>
#include <cstdlib>

EM_REAL  EPS_k2= 1.e-6/* .44408920985006e-15*/;

int PerfectEllepticIntegral(EM_REAL k2, EM_REAL kd2, EM_REAL *K, EM_REAL *E){

static EM_REAL  EPS_k2= 1.e-6/* .44408920985006e-15*/;
static EM_REAL EPS_kd2= .44408920985006e-15;
static EM_REAL  CK1[10]={.30072519903687e-03, .39684709020990e-02
     , .10795990490592e-01, .10589953620989e-01, .75193867218084e-02
     , .89266462945565e-02, .14942029142282e-01, .30885173001900e-01
     , .96573590301742e-01, .13862943611199e+01 };
static EM_REAL CK2[10] = { .66631752464607e-04, .17216147097987e-02
     , .92811603829686e-02, .20690240005101e-01, .29503729348689e-01
     , .37335546682286e-01, .48827155048118e-01, .70312495459547e-01
     , .12499999999764e+00, .50000000000000e+00 };
static EM_REAL CE1[10] = { .32519201550639e-03, .43025377747931e-02
     , .11785841008734e-01, .11841925995501e-01, .90355277375409e-02
     , .11716766944658e-01, .21836131405487e-01, .56805223329308e-01
     , .44314718058337e+00,0. };
static EM_REAL CE2[10] ={ .72031696345716e-04, .18645379184063e-02
     , .10087958494375e-01, .22660309891604e-01, .32811069172721e-01
     , .42672510126592e-01, .58592707184265e-01, .93749995116367e-01
     , .24999999999746e+00,0. };

/* if k2 < EPS_k2, return 0. Else, return 1; */
/* kd2=1. - k2; */

	int i;
	EM_REAL sumK1, sumK2, sumE1, sumE2;
	EM_REAL log_kd2;

	if( k2 >= EPS_k2) {
		if(kd2) {
			log_kd2=log(kd2);
			if(kd2>= EPS_kd2) {
				sumK1=CK1[0];
				sumK2=CK2[0];
				for( i=1; i<10; ++i){
					sumK1=sumK1*kd2+CK1[i];
					sumK2=sumK2*kd2+CK2[i];
				}

				sumE2=0.;
				sumE1=sumE2;
				for( i=0; i<9; ++i) {
					sumE1=kd2*(sumE1+CE1[i]);
					sumE2=kd2*(sumE2+CE2[i]);
				}

				*K=sumK1-log_kd2*sumK2;
				*E=1.+(sumE1-log_kd2*sumE2);

			} else {
				*K=CK1[9]-log_kd2*CK2[9];
				*E=1.0;
			}
		} else {
			*K=13.;
			*E=1.;
		}
			return 1;
	}
	return 0;
}

#define CAT .90992521
#define ZOT1 .00005
void  EllepticIntegral(EM_REAL k2, EM_REAL kd2, EM_REAL *K, EM_REAL *E, EM_REAL omega1, EM_REAL omega2){
	EM_REAL  k, e, kc, ec, kd2_half, en, en_PI;
	EM_REAL SLRG, SSML;
    EM_REAL a[5], b[4], b2[5], a_b[4];
	EM_REAL t, sin_t[4], cotan_t, cotan2, omega, sign_o, sign_t;
	EM_REAL s;
	int i, n, n_mod4, lmt ;


	if(omega2)  {
		sign_o = omega2>0 ? 1. : -1. ;
		omega=fabs(omega2);
		lmt=1;
	} else {
		*E=0.;
		*K=*E;
		if(omega1 == 0.) return;
		sign_o = omega1>0 ? 1. : -1. ;
		omega=fabs(omega1);
		lmt=0;
	}

	a[0]=1.;
    if(k2 >= CAT) b2[0]=k2;
    else b2[0]=kd2;

	for(i=0; i<4; ++i){
		b[i]=sqrt(b2[i]);
		b2[i+1]=a[i]*b[i];
		a[i+1]=.5*(a[i]+b[i]);
		a_b[i]=a[i+1]*(a[i]-b[i]);
	}
    SLRG=2.*(a_b[1]+2.*a_b[2]);
    kd2_half=.5*kd2;
    SSML=(.5+kd2_half)-(SLRG+8.*a_b[3]);
    SLRG=kd2_half+SLRG;

	while(1) {
		en=(int)(omega/PI+ .500005);
		if(en) {
			en_PI=en*PI;
			omega -= en_PI;
			sign_t = omega>0 ? 1. : -1. ;
			t=fabs(omega);

			if(fabs(t-HalfPi) <= ZOT1) {
				t=0.;
				en_PI += HalfPi*sign_t;
				en += .5*sign_t;
			}

			if(k2 >= CAT) {
				if(kd2 ) {
					k=(.5*en/a[3])*log((a[3]*a[2]*a[1]*a[1])*(256./(kd2*kd2)));
					e=2.*en*a[3]+k*SLRG;
				} else {
					k=26.*en;
					e=en+en;
				}
			} else {
				k=en_PI/a[4];
				e=k*SSML;
			}

			if(t >= ZOT1) {
				kc=k;
				ec=e;
			}
		} else {
			t=omega;
			sign_t=1.;
			ec=0.;
			kc=ec;
		}

		if(!en || t >= ZOT1) {
			sin_t[0]=sin(t);
			cotan_t=cos(t)/sin_t[0];

			if(k2 >= CAT) {
				if(kd2) {
					cotan2=cotan_t*cotan_t;
					for(i=1; i<4; ++i) {
						cotan_t=.5*(cotan_t+sqrt(cotan2+2.*a_b[i-1]));
						sin_t[i]=a[i]/sqrt(a[i]*a[i]+cotan2);
					}
					k=log(a[3]*(1.+sin_t[3])/(sin_t[3]*cotan_t))/a[3];
					e=2.*(b[1]*(sin_t[3]-sin_t[1])+2.*b[2]*(sin_t[3]-sin_t[2]))+b[0]*(sin_t[3]-sin_t[0]);
					e += sin_t[3]+k*SLRG;

				} else {
					k=log(tan(.5*(t+HalfPi)));
					e=sin_t[0];
				}

			} else {
				s=0.;
				n=0;
				for(i=0; i<4; ++i){
					n *=2;
					if(cotan_t <=0.)  ++n;
					cotan_t=.5*(cotan_t-b2[i+1]/cotan_t);
					sin_t[0]=sqrt(a[i+1]*a[i+1]+cotan_t*cotan_t);
					n_mod4=n%4;
					if(n_mod4 == 2 || n_mod4 == 3) sin_t[0]=-sin_t[0];
					s += a_b[i]/sin_t[0];
				}
				if(cotan_t <=0.)  ++n;
				k=(atan(a[4]/cotan_t)+n*PI)/(16.*a[4]);
				e=k*SSML+.5*s;
			}

			k=kc+k*sign_t;
			e=ec+e*sign_t;
		}

		if(!lmt) break;
		*K=k*sign_o;
		*E=e*sign_o;
		sign_o = omega1>0 ? 1. : -1. ;
		omega=fabs(omega1);
		lmt=0;
	}


    *K -= k*sign_o;
    *E -= e*sign_o;
 }

 #define sign(x)  (((x)>=0.)? 1. : -1.)
  void ell_fe (double phi, double m, double *K, double *E) {
  double sgn, a, b, c, twon,  phase, cycle, am, eok, cs, f;
  double s;
  int per;

	if (m>=1.) {
 /*		phi= atanh(sin(phi));*/
		s=sin(phi);
		*K=0.5*log((1.+s)/(1.-s));
		per= (int) floor((phi+ HalfPi)/PI);
		*E= sin(phi)*sign(0.5-std::abs(per%2)) + 2.*per;
		a= 0.;
	} else {          /* compute using arithmetic-geometric mean */
		sgn= sign(phi);
		phi= fabs(phi);
		a= 1.;
		b= sqrt(1.-m);
		twon= 1.0;

		eok= 1.-0.5*m;
		cs= 0.;
		for (;;) {      /* maximum of 8 passes for 64-bit double */
			c= 0.5*(a-b);
			if (!c) break;

			phase= (phi+ HalfPi)/PI;
			cycle= floor(phase);
//			phi*= 1. + 1.e-15*(cycle==phase);
//  2012/4/17  By A. Kameari
			phi*= 1. + 1.e-15*(fabs(cycle-phase)<1.e-10);
			phi+= atan((b/a)*tan(phi)) + PI*cycle;

			cs+= c*sin(phi);
			eok-= twon*c*c;

			twon*= 2.0;
			am= a-c;
			if (am==a) break;
			b= sqrt(a*b);
			a= am;
		}
		*K= phi / (twon*a*sgn);

		f= sgn*phi/(twon*a);
		*E= eok*f + sgn*cs;
	}
}

void EllepticIntegral2(EM_REAL k2, EM_REAL kd2, EM_REAL *K, EM_REAL *E, EM_REAL omega1, EM_REAL omega2){
	EM_REAL K1, E1, K2, E2;
	if(omega2)
		ell_fe (omega2, k2, &K2, &E2);
	else {
		K2=0.;
		E2=0.;
	}
	if(omega2==-omega1){
		*K =K2*2.;
		*E =E2*2.;
	}
	else {
		if(omega1)
			ell_fe (omega1, k2, &K1, &E1);
		else {
			K1=0.;
			E1=0.;
		}
		*K=K2-K1;
		*E=E2-E1;
	}
}
