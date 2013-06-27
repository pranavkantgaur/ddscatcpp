#include "StdAfx.h"
#include "FftEngineMkl.h"

FftEngineMkl::FftEngineMkl(void)
{
	method = FftMethod_FFTMKL;
}

FftEngineMkl::~FftEngineMkl(void)
{

}

void FftEngineMkl::Init(unsigned int nX, unsigned int nY, unsigned int nZ)
{

}

void FftEngineMkl::SayHello(char *Buffer)
{
	strcat(Buffer, " - using DFTI from Intel Math Kernel Library (MKL)");
}

void FftEngineMkl::DoFFT(Complex *ccc, unsigned int mx, unsigned int my, unsigned int mz, FftDirection isign)
{
	fprintf(stderr, "FftEngineMkl::DoFFT not implemented\n");
	exit(0);
}
