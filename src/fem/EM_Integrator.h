#pragma once

#include "basic/EM_Interval.h"
#ifndef COILBASE
#include "fem/EM_GaussData.h"
#endif
class EM_Integrand{
public:
	EM_Integrand(int numInput, int numOutput)
		:numInput(numInput), numOutput(numOutput){
		value=new double[numOutput];
	}

	virtual ~EM_Integrand(){
		if(value) delete[] value;
	}
	int NumOutput() const { return numOutput; };

	virtual double *Func(const double * const input)=0;
	int NumInput() const { return numInput; }
protected:


	double *value;
	int numInput, numOutput;
};


class EM_Integrator{
public:
	virtual ~EM_Integrator(){
		if(value) delete[] value;
		if(result) delete[] result;
	}

	int NumOutput() const { return integrand->NumOutput(); }

	virtual double *Func(const double * const input){ return integrand->Func(input); }
	virtual double *Calc(const EM_Interval<double> *range)=0;

protected:
	EM_Integrator(){}
	EM_Integrator(EM_Integrand *integrand):integrand(integrand)
		,numInput(integrand->NumInput())
		,numOutput(integrand->NumOutput())
	{
		value=new double[numOutput];
		result=new double[numOutput];
	}

	EM_Integrand *integrand;
	double *value;
	int numInput, numOutput;
	double *result;
};
#ifndef COILBASE

class EM_GaussIntegrator: public EM_Integrator, protected EM_Gauss1dData{
public:
	EM_GaussIntegrator(EM_Integrand *integrand, int numPoints=3)
		:EM_Integrator(integrand)
		,numPoints(numPoints){
		p=points[numPoints];
		w=weights[numPoints];
	}

	double *Calc(const EM_Interval<double> * range);

private:
	int numPoints;
	const double *w, *p;

	friend class EM_Gauss2dIntegrator;
	friend class EM_Gauss3dIntegrator;
	friend class EM_Gauss2dAdaptive1dIntegrator;

};

class EM_Gauss2dIntegrator: public EM_Integrator{
public:
	EM_Gauss2dIntegrator(EM_Integrand *integrand, int numPointsX=3, int numPointsY=3):
		EM_Integrator(integrand){
		SetPoints(numPointsX, numPointsY); }

	void SetPoints(int numPointsX, int numPointsY){
		this->numPointsX=numPointsX;
		this->numPointsY=numPointsY;
		pX=EM_GaussIntegrator::points[numPointsX];
		wX=EM_GaussIntegrator::weights[numPointsX];
		pY=EM_GaussIntegrator::points[numPointsY];
		wY=EM_GaussIntegrator::weights[numPointsY];
	}
	double *Calc(const EM_Interval<double> * range);
	void SetZ(double z) { this->z=z; }

private:
	int numPointsX, numPointsY;
	const double *wX, *pX, *wY, *pY;
	double z;
};
/*
class EM_3dIntegrator:virtual public EM_Integrator{
public:
	EM_3dIntegrator(EM_Integrand *integrand)
		:EM_Integrator(integrand){}
	virtual ~EM_3dIntegrator(){}
//	EM_3dIntegrator *Create(
};
*/
class EM_Gauss3dIntegrator: public EM_Integrator{
public:
	EM_Gauss3dIntegrator(EM_Integrand *integrand, int numPointsX=3, int numPointsY=3, int numPointsZ=3):
		EM_Integrator(integrand){
		SetPoints(numPointsX, numPointsY, numPointsZ); 
	}
	void SetPoints(int numPointsX, int numPointsY, int numPointsZ){
		this->numPointsX=numPointsX;
		this->numPointsY=numPointsY;
		this->numPointsZ=numPointsZ;
		pX=EM_GaussIntegrator::points[numPointsX];
		wX=EM_GaussIntegrator::weights[numPointsX];
		pY=EM_GaussIntegrator::points[numPointsY];
		wY=EM_GaussIntegrator::weights[numPointsY];
		pZ=EM_GaussIntegrator::points[numPointsZ];
		wZ=EM_GaussIntegrator::weights[numPointsZ];
	}
	double *Calc(const EM_Interval<double> *range);

private:
	int numPointsX, numPointsY, numPointsZ;
	const double *wX, *pX, *wY, *pY, *wZ, *pZ;
};
#endif
class EM_AdaptiveIntegrator: public EM_Integrator{
public:
	class Criterion{
	public:
		Criterion(int n, double eps): n(n), eps(eps) {}
		virtual ~Criterion(){}
		virtual bool Apply(const double *oldval, const double *newval, const double *maxval) const =0;
	protected:
		int n;
		double eps;
	};

	class EachCriterion: public Criterion{
	public:
		EachCriterion(int n, double eps): Criterion(n, eps){}
		bool Apply(const double *oldval, const double *newval, const double *maxval) const{
			for(int i=0; i<n; ++i){
				if(fabs(newval[n]-oldval[n]) > eps*fabs(maxval[n])) return false;
			}
			return true;
		}
	};

	class VectorCriterion: public Criterion{
	public:
		VectorCriterion(int n, double eps): Criterion(n, eps){}
		bool Apply(const double *oldval, const double *newval, const double *maxval) const{
			double sum=0., maxsum=0.;
			for(int i=0; i<n; ++i){
				sum +=(newval[n]-oldval[n])*(newval[n]-oldval[n]);
				maxsum += maxval[n]*maxval[n];
			}
			if(sum > maxsum*eps*eps) return false;
			return true;
		}
	};

	EM_AdaptiveIntegrator(EM_Integrand *integrand, double eps, int kmax=7):
		EM_Integrator(integrand)
		,kmax(kmax<=kMAX? kmax: kmax){
		int n=integrand->NumOutput();
		Set(n);
		criterion=new EachCriterion(n, eps); 
	}

	EM_AdaptiveIntegrator(EM_Integrand *integrand, const Criterion *criterion, int kmax=7):
		EM_Integrator(integrand)
		,kmax(kmax<=kMAX? kmax: kmax), criterion(criterion){
		int n=integrand->NumOutput();
		Set(n);
	}

	void Set(int n){
		rold=new double[n];
		rnew=new double[n];
		vals=new double[n*(int)pow(2.,kmax)];
		maxval=new double[n];
	}

	~EM_AdaptiveIntegrator(){
		if(rold) delete[] rold;
		if(rnew) delete[] rnew;
		if(vals) delete[] vals;
		if(maxval) delete[] maxval;
		if(criterion) delete criterion;
	}

	virtual double *Calc(const EM_Interval<double> * range);
	int GetkResult() const { return kResult; } 

	static const double P[381];
private:

	static const int kMAX=7;
	int kmax, kResult;
	double *rold, *rnew, *vals, *maxval;
	const Criterion *criterion;
	friend class EM_Gauss2dAdaptive1dIntegrator;
	const EM_Interval<double> *range;
};
#ifndef COILBASE
class EM_Gauss2dAdaptive1dIntegrator: public EM_AdaptiveIntegrator{
public:
	EM_Gauss2dAdaptive1dIntegrator(EM_Integrand *integrand, 
		int numPointsX=3, int numPointsY=3, int kmax=6, double eps=1.e-6)
		:EM_AdaptiveIntegrator(integrand, eps, kmax)
//		, eps(eps)
	{
		intergator2d=new EM_Gauss2dIntegrator(integrand,numPointsX, numPointsY);
	}

	~EM_Gauss2dAdaptive1dIntegrator(){
		if(intergator2d) delete intergator2d;
	}
	double *Calc(const EM_Interval<double> *range){
		this->range=range;
		return EM_AdaptiveIntegrator::Calc(range);
	}

	virtual double *Func(const double * const input){
		intergator2d->SetZ(*input);
		return intergator2d->Calc(range);
	}

private:
	void SetPoints(int numPointsX, int numPointsY){
		intergator2d=new EM_Gauss2dIntegrator(integrand,numPointsX, numPointsY);
	}
//	static const double *P;
//	int kmax;
	EM_Gauss2dIntegrator *intergator2d;
//	double *rold, *rnew, *vals;
//	double eps;
};

#endif