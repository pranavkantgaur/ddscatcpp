#include "StdAfx.h"
#include "FftEngineFftw.h"

/* **
GPFAPACK - FORTRAN IMPLEMENTATION OF THE SELF-SORTING IN-PLACE GENERALIZED PRIME FACTOR (COMPLEX) FFT [GPFA]
WRITTEN BY CLIVE TEMPERTON RECHERCHE EN PREVISION NUMERIQUE / ECMWF
THE PACKAGE CONSISTS OF THE SETUP ROUTINE SETGPFA, TOGETHER WITH THE ROUTINES GPFA, GPFA2F, GPFA3F, GPFA5F
** */

FftEngineFftw::FftEngineFftw(void)
{
	method = FftMethod_FFTW21;
}

FftEngineFftw::~FftEngineFftw(void)
{

}

void FftEngineFftw::Init(unsigned int nX, unsigned int nY, unsigned int nZ)
{

}

void FftEngineFftw::SayHello(char *Buffer)
{
	strcat(Buffer, " - using FFTW from Frigo & Johnson");
}

void FftEngineFftw::DoFFT(Complex *ccc, unsigned int mx, unsigned int my, unsigned int mz, FftDirection isign)
{
	fprintf(stderr, "FftEngineFftw::DoFFT not implemented\n");
	exit(0);
}

