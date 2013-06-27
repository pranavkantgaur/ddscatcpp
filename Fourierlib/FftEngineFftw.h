#ifndef __FFTENGINEFFTW_H__
#define __FFTENGINEFFTW_H__

#include "AbstractFftEngine.h"

class FOURIERLIB_API FftEngineFftw : public AbstractFftEngine
{
public:
	FftEngineFftw(void);
	virtual ~FftEngineFftw(void);

public:
	virtual void Init(unsigned int nX, unsigned int nY, unsigned int nZ);
	virtual void SayHello(char *Buffer);
	virtual void DoFFT(Complex *ccc, unsigned int mx, unsigned int my, unsigned int mz, FftDirection isign);
};

#endif
