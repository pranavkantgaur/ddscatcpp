#ifndef __FFTENGINECUDA_H__
#define __FFTENGINECUDA_H__

#include "AbstractFftEngine.h"

class FOURIERLIB_API FftEngineCuda : public AbstractFftEngine
{
public:
	FftEngineCuda(void);
	virtual ~FftEngineCuda(void);

public:
	virtual void Init(unsigned int nX, unsigned int nY, unsigned int nZ);
	virtual void SayHello(char *Buffer);
	virtual void DoFFT(Complex *ccc, unsigned int mx, unsigned int my, unsigned int mz, FftDirection isign);
};

#endif
