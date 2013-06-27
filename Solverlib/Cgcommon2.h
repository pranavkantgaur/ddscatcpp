#ifndef __CGCOMMON2_H__
#define __CGCOMMON2_H__

#include "Solverlib.h"
#include "Complex.h"

SOLVERLIB_API real Scsetrhsstop(Complex *b, Complex *wrk, real epsilon, int *ipar, void (*Preconl)(Complex *, Complex *, int *), real (*Pscnrm)(int, Complex *));
SOLVERLIB_API void Progress(int loglen, int itno, real normres, Complex *x, Complex *res, Complex *trveres);
SOLVERLIB_API real Pscnrm2(int loclen, Complex *u);
SOLVERLIB_API real Scnrm2(int n, Complex *cx, int incx);
SOLVERLIB_API void Pcsum(int isize, Complex *x);

#endif