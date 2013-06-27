#ifndef __FFTENGINEGPFAFTNEW_H__
#define __FFTENGINEGPFAFTNEW_H__

#include "AbstractFftEngine.h"

#define II(JB,JA,JSTEPL) { JB = JA + JSTEPL; if (JB < istart) JB += ninc;}

// To replace FftEngineGpfaft totally
class FftEngineGpfaftNew : public AbstractFftEngine
{
protected:
	int nx, ny, nz;
	Complex *trigX, *trigY, *trigZ;
	int mxold, myold, mzold;

public:
	FftEngineGpfaftNew(void);
	virtual ~FftEngineGpfaftNew(void);

protected:
	void LoadTrigs(Complex *cs, int num, Complex *trigs, int kk, real s);

public:
	void Gpfa2f(Complex *ab, Complex *trigs, int inc, int jump, int n, int mm, int lot, FftDirection isign);
	void Gpfa3f(Complex *ab, Complex *trigs, int inc, int jump, int n, int mm, int lot, FftDirection isign);
	void Gpfa5f(Complex *ab, Complex *trigs, int inc, int jump, int n, int mm, int lot, FftDirection isign);

protected:
	void Gpfa3fStep1(Complex *tu, Complex *ab, int j1, int j2, int j3, real c1);
	void Gpfa3fStep3A(Complex *tu, Complex *ab, int j1, int j2);
	void Gpfa3fStep3B(Complex *tu, Complex *ab, int j1, int j2, Complex *cs);
	void Gpfa5fStep1(Complex *tu, Complex *ab, int j1, int j2, int j3, int j4, real c1);
	void Gpfa5fStep2(Complex *tu, Complex *ab, int j1, int j5, int j3, int j7, int j6, int j0, real c2, real c3);
	void Gpfa5fStep2Simple(Complex *tu, Complex *ab, int j1, real c2, real c3);
	void Gpfa5fStep3A(Complex *tu, Complex *ab, int j1, int j2, int j3, int j4);
	void Gpfa5fStep3B(Complex *tu, Complex *ab, int j1, int j2, int j3, int j4, Complex *cs);
};

#endif