#include "StdAfx.h"
#include "AbstractFftEngine.h"
#include "FftEngineGpfaft.h"
#include "FftEngineFftw.h"
#include "FftEngineMkl.h"
#include "FftEngineCuda.h"

AbstractFftEngine *AbstractFftEngine::engine = NULL;
FftMethod AbstractFftEngine::method = FftMethod_End;

AbstractFftEngine::AbstractFftEngine(void)
{
	method = FftMethod_End;
	engine = NULL;
	nCashed = 0;
	njCashed[0] = njCashed[1] = njCashed[2] = 0;
	nx = ny = nz = 0;
}

AbstractFftEngine::~AbstractFftEngine(void)
{
	if (engine)
		delete engine;
}

AbstractFftEngine *AbstractFftEngine::GetEngine(FftMethod method)
{
	if (!engine)
	{
		if (method != FftMethod_End)
		{
			switch(method)
			{
			case FftMethod_GPFAFT:
				engine = new FftEngineGpfaft;
				break;
		
			case FftMethod_FFTW21:
				engine = new FftEngineFftw;
				break;

			case FftMethod_FFTMKL:
				engine = new FftEngineMkl;
				break;

			case FftMethod_CUDAFX:
				engine = new FftEngineCuda;
				break;

			default:
				engine = NULL;
				break;
			}
		}		
	}
	return engine;
}

FftMethod AbstractFftEngine::GetMethod()
{
	return method;
}

void AbstractFftEngine::ClearEngine(void)
{
	Kill();
}

void AbstractFftEngine::Kill(void)
{
	if (engine)
	{
		delete engine;
		engine = NULL;
		method = FftMethod_End;
	}	
}

unsigned int AbstractFftEngine::powint(unsigned int x, unsigned int y)
{
	if (y == 1)
		return x;
	if (y < 1)
		return 1;

	int res = 1;
	for(unsigned int i=0; i<y; ++i)
	{
		res *= x;
	}
	return res;
}

//
// Decompose n into factors 2, 3, 5
bool AbstractFftEngine::PrimeDecompose(unsigned int n, unsigned int *nj)
{
	if (n != nCashed)
	{
		int nn = n;
		int ifac = 2;

		for(int ll=0; ll<3; ++ll)
		{
			unsigned int kk = 0;
			while(1)
			{
				if (nn % ifac != 0) 
					break;
				++kk;
				nn /= ifac;
			}
			nj[ll] = kk;
			ifac += (ll + 1);
		}

		if (nn != 1)
			return false;

		nCashed = n;
		memcpy(njCashed, nj, 3*sizeof(unsigned int));
	}
	else
	{
		memcpy(nj, njCashed, 3*sizeof(unsigned int));		
	}
	return true;
}
