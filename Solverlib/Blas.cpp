#include "StdAfx.h"
#include "Blas.h"

// Copies a vector, x, to a vector, y.
// Jack Dongarra, Linpack, 4/11/78.
void Ccopy(int n, Complex *cx, int incx, Complex *cy, int incy)
{
	if (n <= 0)
		return;
	if ((incx == 1) && (incy == 1))				// Code for both increments equal to 1
	{
		for(int i=0; i<n; ++i)
		{
			cy[i] = cx[i];
		}
	}
	else										// Code for unequal increments or equal increments not equal to 1
	{
		int ix = 1;
		int iy = 1;
		if (incx < 0)
			ix = (-n+1)*incx + 1;
		if (incy < 0) 
			iy = (-n+1)*incy + 1;
		for(int i=0; i<n; ++i)
		{
			cy[iy] = cx[ix];
			ix += incx;
			iy += incy;
		}
	}
}

// Forms the dot product of a vector.
// Jack Dongarra, 3/11/78.
Complex Cdotc(int n, Complex *cx, int incx, Complex *cy, int incy)
{
	Complex ctemp;

	if (n <= 0)
		return ctemp;

	if ((incx == 1) && (incy == 1))
	{
		for(int i=0; i<n; ++i)
		{
			ctemp += cx[i].conjg() * cy[i];
		}
	}
	else
	{
		int ix = 1;
		int iy = 1;
		if (incx < 0) 
			ix = (-n+1)*incx + 1;
		if (incy < 0) 
			iy = (-n+1)*incy + 1;
		for(int i=0; i<n; ++i)
		{
			ctemp += cx[ix].conjg() * cy[iy];
			ix += incx;
			iy += incy;
		}
	}
	return ctemp;
}
//
// Constant CA times a vector CX plus a vector CY, result in CY
// Jack Dongarra, linpack, 3/11/78.
// modified 12/3/93, array(1) declarations changed to array(*)
void Caxpy(int n, Complex ca, Complex *cx, int incx, Complex *cy, int incy)
{
	if (n <= 0) return;
	if (Fabs(ca.re) + Fabs(ca.im) == (real)0.) 
		return;
	if ((incx == 1) && (incy == 1))				// Code for both increments equal to 1
	{
		for(int i=0; i<n; ++i)
		{
			cy[i] += ca * cx[i];
		}
	}
	else										// Code for unequal increments or equal increments not equal to 1
	{
		int ix = 1;
		int iy = 1;
		if (incx < 0) ix = (-n+1) * incx + 1;
		if (incy < 0) iy = (-n+1) * incy + 1;
		for(int i=0; i<n; ++i)
		{
	        cy[iy] += ca * cx[ix];
		    ix += incx;
			iy += incy;
		}
	}
}

//
// Scales a vector by a constant.
// Jack Dongarra, 3/11/78.
// Modified to correct problem with negative increment, 8/21/90.
void Cscal(int n, Complex ca, Complex *cx, int incx)
{
	if (n <= 0) 
		return;
	if (incx == 1)						// Code for increment equal to 1
	{
		for(int i=0; i<n; ++i)
		{
			cx[i] *= ca;
		}
	}
	else								// Code for increment not equal to 1
	{
		int ix = 1;
		if (incx < 0) ix = (-n+1) * incx + 1;
		for(int i=0; i<n; ++i)
		{
			cx[ix] *= ca;
			ix += incx;
		}
	}
}

//
// Interchanges two vectors.
// Jack Dongarra, 3/11/78.
void Cswap(int n, Complex *cx, int incx, Complex *cy, int incy)
{
	if (n <= 0) 
		return;
	if ((incx == 1) && (incy == 1))					// Code for both increments equal to 1
	{
		for(int i=0; i<n; ++i)
		{
			Complex ctemp = cx[i];
			cx[i] = cy[i];
			cy[i] = ctemp;
		}
	}
	else											// Code for unequal increments or equal increments not equal to 1
	{
		int ix = 1;
		int iy = 1;
		if (incx < 0) 
			ix = (-n+1) * incx + 1;
		if (incy < 0) 
			iy = (-n+1) * incy + 1;
		for(int i=0; i<n; ++i)
		{
			Complex ctemp = cx[ix];
			cx[ix] = cy[iy];
			cy[iy] = ctemp;
			ix += incx;
			iy += incy;
		}
	}
}

//
// Initialises a vector x with a scalar alpha.
// Modified from scopy, BLAS Level 1.
// Rudnei Dias da Cunha, 14/6/93.
// copies a vector, x, to a vector, y.
// uses unrolled loops for increments equal to one.
// Jack Dongarra, Linpack, 3/11/78.
void Cinit(int n, Complex alpha, Complex *cx, int incx)
{
	if (n <= 0)
		return;
	if (incx == 1)							// Code for both increments equal to 1, clean-up loop
	{
		int m = n%8;
		if (m)
		{
			for(int i=0; i<m; ++i)
			{
				cx[i] = alpha;				
			}
			if (n < 8) return;
		}
		for(int i=m; i<n; i+=8)
		{
			cx[i] = alpha;
			cx[i+1] = alpha;
			cx[i+2] = alpha;
			cx[i+3] = alpha;
			cx[i+4] = alpha;
			cx[i+5] = alpha;
			cx[i+6] = alpha;
			cx[i+7] = alpha;
		}
	}
	else									// Code for unequal increments or equal increments, not equal to 1
	{
		int ix = 0;
		if (incx < 0) 
			ix = (-n+1)*incx;
		for(int i=0; i<n; ++i)
		{
			cx[ix] = alpha;
			ix += incx;
		}
	}
}
