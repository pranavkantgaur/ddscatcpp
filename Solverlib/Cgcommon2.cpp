#include "StdAfx.h"

#include "Definitions.h"
#include "Complex.h"
#include "Blas.h"
#include "Cgcommon2.h"

real Scsetrhsstop(Complex *b, Complex *wrk, real epsilon, int *ipar, void (*Preconl)(Complex *, Complex *, int *), real (*Pscnrm)(int, Complex *))
{	
	int loclen = ipar[3];
	int stoptype = ipar[8];
	
	switch(stoptype)	
	{	
		case 1:									// ||r||<epsilon or ||Q1r||<epsilon ||x(k)-x(k-1)||<epsilon	
		case 4:	
		case 7:		
			return epsilon;		
			break;	

		case 2:									// ||r||<epsilon||b|| or sqrt(r(Q1r))<epsilon||b|| or ||Q1r||<epsilon||b||	
		case 3:	
		case 5:		
			return epsilon * Pscnrm(loclen, b);
			break;	

		case 6:									// ||Q1r||<epsilon||Q1b||		
			Preconl(b, wrk, NULL);        
			return epsilon * Pscnrm(loclen, wrk);
			break;	

		default:		
			return 0;		
			break;	
	}
}

void Stopcrit(Complex *b, Complex *r, Complex *rtrue, Complex *x, Complex *xold, Complex *wrk, real rhsstop, int cnvrtx, real &exitnorm, int &status, int *ipar, 
	void (*Matvec)(Complex *, Complex *, int *), void (*Tmatvec)(Complex *, Complex *, int *), void (*Preconr)(Complex *, Complex *, int *), 
	void (*Pcsum)(int, Complex *), real (*Pscnrm)(int, Complex *))
{
	const Complex coner((real)1., (real)0.);
	int loclen = ipar[3];
	int precontype = ipar[7];
	int stoptype = ipar[8];

	switch(stoptype)
	{
	case 1:
	case 2:
	case 3:											// Compute true residual if needed
		Ccopy(loclen, b, 1, rtrue, 1);
		if ((precontype == 2) || (precontype == 3))
		{
			Preconr(x, wrk, NULL);
			if (cnvrtx == 1)
			{
				Tmatvec(wrk, xold, NULL);
				Matvec(xold, wrk, NULL);
				Caxpy(loclen, -coner, wrk, 1, rtrue, 1);
			}
			else
			{
				Matvec(wrk, xold, NULL);
				Caxpy(loclen, -coner, xold, 1, rtrue, 1);
			}
		}
        else
		{
			if (cnvrtx == 1)
			{
				Tmatvec(x, xold, NULL);
				Matvec(xold, wrk, NULL);
				Caxpy(loclen, -coner, wrk, 1, rtrue, 1);
			}
			else
			{
				Matvec(x, wrk, NULL);
				Caxpy(loclen, -coner, wrk, 1, rtrue, 1);
			}
		}
//
		if ((stoptype == 1) || (stoptype == 2))					// ||r|| < epsilon or ||r|| < epsilon ||b||
		{
			exitnorm = Pscnrm(loclen, rtrue);
			status = (exitnorm < rhsstop) ? 0 : -99;
		}
		else
		{
			if (stoptype == 3)									// sqrt(r(Q1r))<epsilon||b||
			{
				Complex dots[1];
				Complex temp = Cdotc(loclen, rtrue, 1, r, 1);
				dots[0] = Complex(temp.re, (real)0.);
				Pcsum(1, dots);
				exitnorm = dots[0].sqrt().mod();
				status = (exitnorm < rhsstop) ? 0 : -99;
			}
		}
		break;

	case 4:
	case 5:
	case 6:										// ||Q1r|| < epsilon or ||Q1r|| < epsilon||b|| or ||Q1r|| < epsilon||Q1b
		exitnorm = Pscnrm(loclen, r);
		status = (exitnorm < rhsstop) ? 0 : -99;
		break;

	case 7:										// ||x-x0||<epsilon
		Ccopy(loclen, x, 1, wrk, 1);
		Caxpy(loclen, -coner, xold, 1, wrk, 1);
        exitnorm = Pscnrm(loclen, wrk);
		status = (exitnorm < rhsstop) ? 0 : -99;
		break;

	default:
		break;
	}
}


void Progress(int loglen, int itno, real normres, Complex *x, Complex *res, Complex *trveres)
{
/* **

    const char *Format9000 = " iter= %4d frac.err= %11.7lf\n";
    const char *Format9001 = " iter= %4d frac.err= %11.7lf min.err=%11.7lf rate=%8.6lf\n";

    static real errmin = (real)0.;

      USE DDCOMMON_9,ONLY : ERRSCAL,IDVOUT2,ITERMX,ITERN

! Arguments:
      REAL (WP) :: NORMRES
      INTEGER :: ITNO, LOCLEN
      COMPLEX (WP) :: RES(*), TRUERES(*), X(*)

! Common:
!      INTEGER :: IDVOUT2, ITERMX, ITERN
!      REAL (WP) :: ERRSCAL
!-----------------------------------------------------------------------
!      COMMON /NORMERR/ERRSCAL, IDVOUT2, ITERMX, ITERN
!-----------------------------------------------------------------------

! Local:
      INTEGER :: ITNOL
      REAL (WP) :: ERRMIN, NORMERR, RATE
      SAVE ERRMIN, ITNOL

// Part of PIM
History:
Fortran versions history removed

! Note: when used with STOPTYPE=5, NORMRES=|Ax-b|.
! If quantity ERRSCAL is set to |b| elsewhere, then NORMRES/ERRSCAL
! is a measure of the fractional error in the solution vector x.

** */

	/* **
	++Common9::GetInstance()->Itern();
	normerr = normres / errscal;

	if (itern <= 2)
	{
		real errmin = normerr;
		itnol = itern;
		fprintf(stderr, Format9000, itern, normerr);
	}
	else
	{
		if (normerr < errmin)
		{
			rate = log(errmin / normerr) / (itern - itnol);
            errmin = normerr;
            itnol = itern;
			fprintf(stderr, Format9001, itern, normerr, errmin, rate);
		}
		else
		{
			fprintf(stderr, Format9000, itern, normerr);
		}
	}
	** */
}

real Pscnrm2(int loclen, Complex *u)
{
	return Scnrm2(loclen, u, 1);
}

real Scnrm2(int n, Complex *cx, int incx)
{
/* **
Returns SCNRM2=unitary norm of Complex n-vector stored in CX with
storage increment INCX.
Written to replace SCNRM2 written by C.L.Lawson, as g77 compiler
produces bad code when optimizing Lawson's routine.
In any event, Lawson's routine appears to be unnecessarily complicated
for needs of DDSCAT.

History:
Fortran versions history removed
end history
** */
//
	real sum = (real)0.;
	for(int i=0; i<n; ++i)
	{
		int j = i*incx;
		sum += (cx[j] * cx[j].conjg()).re;
	}

	return Sqrt(sum);
}

void Pcsum(int isize, Complex *x)
{
// ChB: it was empty
}
