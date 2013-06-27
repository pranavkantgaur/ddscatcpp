#ifndef __FUNCTIONS_GENERAL_H__
#define __FUNCTIONS_GENERAL_H__

#include "Definitions.h"
#include "Vect3.h"
#include "Complex.h"

// Extracts num real numbers from Buffer according to Format
GENERALLIB_API void ExtractFromInt(const char *Buffer, const char *Format, int num, real *shpar);

// Extracts num real numbers from Buffer according to Format
GENERALLIB_API void ExtractFromReal(const char *Buffer, const char *Format, int num, real *shpar);

// Get nearest 2*3*5 factorizable number
GENERALLIB_API Vect3<int> GetNearestFactoredVector(const Vect3<int> &op);
GENERALLIB_API int GetNearestFactor(int op);

GENERALLIB_API void Wrimsg(const char *csubrt, const char *cmsgnm);
GENERALLIB_API void Errmsg(const char *cstats, const char *csubrt, const char *cmsgnm);

GENERALLIB_API const char *NumberInWords(int number);
GENERALLIB_API int ExtractFirstWord(char *Buffer, char fc = ' ');
GENERALLIB_API void ShiftLeft(char *Buffer);
GENERALLIB_API void RemoveSymbols(char *Buffer, char sym, char toSymbol = ' ');

GENERALLIB_API Complex Scalar(Complex *c, real *d);
GENERALLIB_API Complex Scalar(Complex *c, Complex *d);
GENERALLIB_API void Vector(Complex *c, real *d, Complex *res);
GENERALLIB_API void Vector(Complex *c, Complex *d, Complex *res);

GENERALLIB_API real Ran3(int idum);

#ifndef _WIN32
char *strupr(char *Buffer);
#endif

#endif