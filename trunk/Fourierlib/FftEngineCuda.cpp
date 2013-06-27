#include "StdAfx.h"
#include "FftEngineCuda.h"

FftEngineCuda::FftEngineCuda(void)
{
	method = FftMethod_CUDAFX;
}

FftEngineCuda::~FftEngineCuda(void)
{

}

void FftEngineCuda::Init(unsigned int nX, unsigned int nY, unsigned int nZ)
{

}

void FftEngineCuda::SayHello(char *Buffer)
{
	strcat(Buffer, " - using CUFFT package from NVidia CUDA");
}

void FftEngineCuda::DoFFT(Complex *ccc, unsigned int mx, unsigned int my, unsigned int mz, FftDirection isign)
{
	fprintf(stderr, "FftEngineCuda::DoFFT not implemented\n");
	exit(0);
}
