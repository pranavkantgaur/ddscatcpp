#ifndef __DEFINITIONS_H__
#define __DEFINITIONS_H__

#include "General.h"

#ifdef DOUBLEPRECISION
	#define real double
	#define realFormat "%lf"
	#include "MathDouble.h"
#else
	#define real float
	#define realFormat "%f"
	#include "MathFloat.h"
#endif

extern const GENERALLIB_API real Pi;
extern const GENERALLIB_API real TwoPi;
extern const GENERALLIB_API real FourPi;
extern const GENERALLIB_API real Degrad;
extern const GENERALLIB_API int ihuge_;

int GENERALLIB_API nint_(real op);
real GENERALLIB_API sign_(int a);
real GENERALLIB_API sign_(real a);
int GENERALLIB_API min_(int a, int b);
real GENERALLIB_API min_(real a, real b);
int GENERALLIB_API max_(int a, int b);
real GENERALLIB_API max_(real a, real b);

#endif // __DEFINITIONS_H__