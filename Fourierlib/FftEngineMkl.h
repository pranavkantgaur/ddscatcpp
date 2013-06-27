#ifndef __FFTENGINEMKL_H__
#define __FFTENGINEMKL_H__

#include "AbstractFftEngine.h"

class FftEngineMkl : public AbstractFftEngine
{
public:
	FftEngineMkl(void);
	virtual ~FftEngineMkl(void);

public:
	virtual void Init(unsigned int nX, unsigned int nY, unsigned int nZ);
	virtual void SayHello(char *Buffer);
	virtual void DoFFT(Complex *ccc, unsigned int mx, unsigned int my, unsigned int mz, FftDirection isign);
};

#endif
