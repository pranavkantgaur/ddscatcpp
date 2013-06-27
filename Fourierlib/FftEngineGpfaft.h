#ifndef __FFTENGINEGPFAFT_H__
#define __FFTENGINEGPFAFT_H__

#include "AbstractFftEngine.h"

#define II(JB,JA,JSTEPL) { JB = JA + JSTEPL; if (JB < istart) JB += ninc;}

class FOURIERLIB_API FftEngineGpfaft : public AbstractFftEngine
{
protected:
	Complex *trigX, *trigY, *trigZ;
	unsigned int mxold, myold, mzold;

public:
	FftEngineGpfaft(void);
	virtual ~FftEngineGpfaft(void);

public:
	virtual void Init(unsigned int nX, unsigned int nY, unsigned int nZ);
	virtual void SayHello(char *Buffer);
	void DoFFT(Complex *ccc, unsigned int mx, unsigned int my, unsigned int mz, FftDirection isign);

protected:
	bool Gpfa(Complex *ab, Complex *trigs, int inc, int jump, unsigned int n, int lot, FftDirection isign);
	void Gpfa2f(Complex *ab, Complex *trigs, int inc, int jump, unsigned int n, unsigned int mm, int lot, FftDirection isign);
	void Gpfa3f(Complex *ab, Complex *trigs, int inc, int jump, unsigned int n, unsigned int mm, int lot, FftDirection isign);
	void Gpfa5f(Complex *ab, Complex *trigs, int inc, int jump, unsigned int n, unsigned int mm, int lot, FftDirection isign);
	bool Setgpfa(Complex *trigs, unsigned int n);
	void ElementaryInit(unsigned int &nx, unsigned int nX, Complex *&trig);
	void LoadTrigs(Complex *cs, unsigned int num, Complex *trigs, unsigned int kk, real s);
};

#endif
