#pragma once

#include "AbstractTarget.h"
#include "LoadableTarget.h"
#include "AbstractDFData.h"
#include "DDscatMain.h"
#include "DDscatParameters.h"
#include "Enumerator.h"
#include "DielectricManager.h"
#include "FourArray.h"
#include "FileNamer.h"
#include "Functions.h"
#include "TimerManager.h"

class OutputManager
{
protected:
	AbstractTarget *currentTarget;
	DDscatParameters *param;
	DielectricManager *dielec;
	PeriodicBoundaryFlag jpbc;
	FILE *file8, *file10, *file11, *file12;
	real betmid, betmxd, phimid, phimxd, thtmid, thtmxd;
	real betad, phid, thetad;
	OutputManager(void);
	virtual ~OutputManager(void);

protected:
	static OutputManager *item;

public:
	static OutputManager *GetInstance(void);
	static void Kill(void);

public:
	inline void SetPhid(real p) { phid = p; }
	inline void SetBetad(real b) { betad = b; }
	inline void SetThetad(real t) { thetad = t; }

public:
	void Init(AbstractTarget *currentTarget, DDscatParameters *param, DielectricManager *dielec, PeriodicBoundaryFlag jpbc);
	bool WriteHeadings(int nori);
	void WriteDielec(real wave);
	bool Writefml(int iori, int irad, int iwav, int navg, int ncomp, char *cdescr, real aeff, real ak1, Vect3<real> &akr, 
		real wave, real xx, Vect3<Complex> &cxe01r, Vect3<Complex> &cxe02r, FourArray *cxfData);
	bool Writesca(int itheta, int ibeta, int iphi, bool wantIobin, int iori, int irad, int iwav, int navg, int *itnum, int ncomp, int nori, 
		char *cdescr, real aeff, real ak1, Vect3<real> &ak_tf, real wave, real xx, 
		SumPackage &sumPackage, real *s1111, real *s2121, real *sm, real *smori, Complex *cx1121, Vect3<Complex> &cxe01_tf, 
		Vect3<Complex> &cxe02_tf, FourArray *cxfData, real pyddx, real pzddx, Vect6<real> &xMinmax);
	void Writepol(real aeff, Vect3<real> &akr, real wave, Vect3<Complex> &cxe0r, Matrix *theMatrix, Complex *cxpol, const char *cflpol);
	bool Writebin(bool wantIobin, char *cdescr, int navg, int ncomp, int nori, int iwav, int irad, int iori, real aeff, real wave, 
		real xx, real ak1, Vect3<real>&akr, Vect3<Complex>&cxe01r, Vect3<Complex>&cxe02r, SumPackage &sumPackage, real *g, real *qav, real *sm);
};
