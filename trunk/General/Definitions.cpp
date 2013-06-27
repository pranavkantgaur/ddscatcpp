#include "StdAfx.h"
#include "Definitions.h"

#ifdef DOUBLEPRECISION
	const real Pi = 3.14159265358979323846;
	const real Degrad = 180. / Pi;
#else
	const real Pi = 3.1415926f;
	const real Degrad = 180.f / Pi;
#endif

const real TwoPi = Pi + Pi;
const real FourPi = TwoPi + TwoPi;
const int ihuge_ = (int)0x7fffffff;

int nint_(real op)
{
	int x = (int)(Fabs(op) + (real)0.5);
	return (op < (real)0.) ? -x : x;
}

real sign_(int a)
{
	return (a < 0) ? -(real)1. : (real)1.;
}

real sign_(real a)
{
	return (a < (real)0.) ? -(real)1. : (real)1.;
}

int min_(int a, int b)
{
	return (a < b) ? a : b;
}

real min_(real a, real b)
{
	return (a < b) ? a : b;
}

int max_(int a, int b)
{
	return (a > b) ? a : b;
}

real max_(real a, real b)
{
	return (a > b) ? a : b;
}
