#ifndef __BLAS_H__
#define __BLAS_H__

#include "Solverlib.h"
#include "Definitions.h"
#include "Complex.h"

void Ccopy(int n, Complex *cx, int incx, Complex *cy, int incy);
Complex Cdotc(int n, Complex *cx, int incx, Complex *cy, int incy);
void Caxpy(int n, Complex ca, Complex *cx, int incx, Complex *cy, int incy);
void Cscal(int n, Complex ca, Complex *cx, int incx);
void Cswap(int n, Complex *cx, int incx, Complex *cy, int incy);
void SOLVERLIB_API Cinit(int n, Complex alpha, Complex *cx, int incx);

#endif // __BLAS_H__