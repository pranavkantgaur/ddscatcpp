#ifndef __ABSTRACTFFTENGINE_H__
#define __ABSTRACTFFTENGINE_H__

#include "Fourierlib.h"
#include "Definitions.h"
#include "Enumerator.h"
#include "Complex.h"

typedef enum _FftDirection
{
	FftForward = +1, FftBackward = -1
} FftDirection;

class FOURIERLIB_API AbstractFftEngine
{
protected:
	static FftMethod method;
	static AbstractFftEngine *engine;
	AbstractFftEngine(void);
	virtual ~AbstractFftEngine(void);
	unsigned int nCashed, njCashed[3];
	unsigned int nx, ny, nz;

public:
	static AbstractFftEngine *GetEngine(FftMethod method);
	static void ClearEngine(void);
	static void Kill(void);
	static FftMethod GetMethod();

public:
	virtual void Init(unsigned int nx, unsigned int ny, unsigned int nz) = 0;
	virtual void SayHello(char *Buffer) = 0;
	virtual void DoFFT(Complex *ccc, unsigned int mx, unsigned int my, unsigned int mz, FftDirection isign) = 0;

protected:
	unsigned int powint(unsigned int x, unsigned int y);
	bool PrimeDecompose(unsigned int n, unsigned int *nj);
};

#endif
