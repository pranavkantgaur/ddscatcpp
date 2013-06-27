#include "StdAfx.h"

#include "SolverTangcg.h"
#include "Functions.h"

REGISTER_SOLVER(Tangcg,GPBICG)
SolverTangcg::SolverTangcg(void)
{
	ncompte = itern = 0;
	ndim = maxit = 0;
	xr = NULL;
	steperr = status = 0;
	residu = (real)0.;
	wkFileSize = 12;
}

SolverTangcg::~SolverTangcg(void)
{
	CleanDelete2(xr);
}

void SolverTangcg::SetParameters(int nndim, int mmaxit)
{
	ndim = nndim;
	maxit = mmaxit;
	xr = new Complex[ndim];
	wk = new Complex[wkFileSize * nndim];
}

void SolverTangcg::Tangcg(Complex *xi, Complex *b, real tol, real &tole)
{
/* **
! Interface code to Tang et. al. conjugate gradient
! tol - tolerance to be achived
! maxit - maximum number of iterations
! xi - on input initial guess, on output vector x 
! xr - scratch array
! b - right hand side; i.e. Ax=b
! matvec - external subroutine calculating matrix times vector multiplication Ax
! wrk - scratch array wrk(lda, 12). NOTE that this has to be set properly
! lda - maximum leading dimension
! ndim - actuall dimension of the linear problem
! tole - achived tolerance
! nloop - number of iterations
! ncompte - number of  Ax multiplications (this is where calculations are expensive)
! History
! (PJF) = Piotr Jacek Flatau
! (PCC) = P. C. Chaumet 
! (AR)  = A. Rahmani
! April 2, 2010 Original code by (PCC). (AR)
! May 6, 2010  (PJF)  converted to Fortran90, changed interface,
!                     single/double precision kinds
!-------------------------------------------------------------------------- 
** */
	itern = 0;
	nou = 0;
	ncompte = 0;

	while(1)
	{
		Gpbicg(xi, b, tol);
		if(status < 0)
		{
			fprintf(stderr, "stop nstat %d %d\n", status, steperr);
			return;
		}
		ncompte++;
		Matvec(xi, xr, NULL);
		if(status == 1) 
			break;
	}
//
// Finished iterations
	if (steperr == 0)
	{
		fprintf(stderr, "itern has reached maxit %d %d\n", itern, maxit);
	}
//
// Compute the relative error
	tole = (real)0.;
	for(int i=0; i<ndim; ++i)
	{
		tole += (xr[i] - b[i]).modSquared();
	}
	tole = Sqrt(tole) / norm;
}

/* **
****************************************************************
     Iterative solver GPBICG
****************************************************************
     Authors: P. C. Chaumet and A. Rahmani
     Date: 04/02/2010
     Purpose: iterative solver for linear system Ax=b. There is no
     condition on the matrix. Notice that the product A x is provided
     by the user.
     Reference: if you use this routine in your research, please
     reference, as appropriate: P. C. Chaumet and A. Rahmani, Efficient
     discrete dipole approximation for magneto-dielectric scatterers
     Opt. Lett. 34, 917 (2009). J. Tang, Y. Shen, Y. Zheng, and D. Qiu,
     Coastal Eng. 51, 143 (2004).
     license: GNU GPL
     We cannot guarantee correctness of the programs...
     Main program and example. The main program if provided for testing
     purposes and should be commented out to use only the routine GPBICG
     We want to solve Ax=b
     mat:  matrix A
     lda: Input: leading dimension array of the matrix
     ndim: Input: dimension  of the matrix: ndim.le.lda
     nlar: Input: size of the work vector: nlar.ge.12
     MAXIT: Input: Maximum of iteration for the iterative solver
     NLOOP: Output: number of iteration for the iterative solver Should
     be initialize to 0.
     ncompte: number of Ax products.
     STATUS: Output: if STATUS.lt.0 a problem occured in GPBICG
     STATUS=1 the requested tolerance has been reached or the maximum n
     iterations has been reached
     STEPERR: Output: if STEPERR.lt.0: indicates where the problem
     occurs in GPBICG. STEPERR=0 the maximum number of iterations has
     been reached.  STEPERR=-1 routine completed witho

     tol: Input: tolerance requested by the user. At the end we have:
     r=||Ax-b||/||b||<tol
     b: Input: right-hand side of Ax=b
     norm: Output: norm of b
     xi: Input: initial guess; output:solution of the linear equation
     xr: Input: xr = A xi
     NOU: Local int used by GPBICG. Should be initialized to 0.
     ALPHA,BETA,ETA,DZETA,R0RN: Local complex needs for GPBICG
     WRK: local array used by for GPBICG

History:
 (PJF) = Piotr J. Flatau
 May 6, 2010 - added wrmsg, converted to Fortran 90, added single/double precision kinds
               corrected "40 elseif" statement such that it is on separate lines 
** */
void SolverTangcg::Gpbicg(Complex *xi, Complex *b, real tol)
{
//!     Column index of the various variables stored in WRK array.
//!     1:r0
//!     2:p
//!     3:r
//!     4:y
//!     5:t
//!     6:Ap
//!     7:At
//!     8:u
//!     9:w
//!     10:z
//!     11:x
//!     12:t old

	const int nlar = 12;
	Complex ctmp, ctmp1, ctmp2, ctmp3, ctmp4, ctmp5;
	char cmsgnm[81];
	int ii, jj;
	switch(nou)
	{
	case 1:												// initialize r0=b-Ax0,rOt=rO,p0=v0=d0
		for(ii=jj=0; ii<ndim; ++ii, jj+=nlar)
		{
			wk[jj    ] = b[ii] - xr[ii];
			wk[jj + 1] = czero;
			wk[jj + 2] = wk[jj];
			wk[jj + 4] = czero;
			wk[jj + 8] = czero;
			wk[jj + 7] = czero;
			wk[jj + 9] = czero;
		}
//
// Compute the initial residue
		ctmp = czero;
		for(ii=jj=0; ii<ndim; ++ii, jj+=nlar)
		{
			ctmp += wk[jj] * wk[jj].conjg();
		}
		r0rn = ctmp;
//
// Initialize rho,alpha,w=1,tau=norm,theta=0,eta=0
		beta = czero;
//
// Begin the iteration sequence
		itern = -1;
//
// btd add wrimsg call
		residu = (real)1.;
		sprintf(cmsgnm, "IT=%8d f.err=%10.3e", itern + 1, residu);
		Wrimsg("Gpbicg", cmsgnm);
//
l100:	++itern;
//
// Compute p=r+beta*(p-u)
		for(ii=jj=0; ii<ndim; ++ii, jj+=nlar)
		{
			wk[jj + 1] = wk[jj + 2] + beta * (wk[jj + 1] - wk[jj + 7]);
			xi[ii] = wk[jj + 1];
		}
		nou = 2;
		break;
//
// Compute Ap
	case 2:
		for(ii=jj=0; ii<ndim; ++ii, jj+=nlar)
		{
			wk[jj + 5] = xr[ii];
		}
// 
// Compute alpha=r0r/r0Ap
		ctmp = czero;
		for(ii=jj=0; ii<ndim; ++ii, jj+=nlar)
		{
			ctmp += wk[jj].conjg() * wk[jj + 5];
		}
		if (ctmp == czero)
		{
			status = -1;
			steperr = 1;
			return;
		}
		alpha = r0rn/ctmp;
//
// Compute y=t-r-alpha*w+alpha*Ap et de t=r-alpha AP
		for(ii=jj=0; ii<ndim; ++ii, jj+=nlar)
		{
			wk[jj +  3] = wk[jj + 4] - wk[jj + 2] - alpha * wk[jj + 8] + alpha * wk[jj + 5];
			wk[jj + 11] = wk[jj + 4];
			wk[jj +  4] = wk[jj + 2] - alpha * wk[jj + 5];
			xi[ii] = wk[jj + 4];
		}
		nou = 3;
		break;
//
// Compute At
	case 3:
		for(ii=jj=0; ii<ndim; ++ii, jj+=nlar)
		{
			wk[jj + 6] = xr[ii];
		}
//
// Compute dzeta and eta
		if (itern == 0)
		{
			eta = czero;
			dzeta = czero;
			ctmp = czero;
			for(ii=jj=0; ii<ndim; ++ii, jj+=nlar)
			{
				dzeta += (wk[jj + 6].conjg() * wk[jj + 4]);
				ctmp  += (wk[jj + 6].conjg() * wk[jj + 6]);
			}

			if (ctmp == czero)
			{
				status = -1;
				steperr = 2;
				return;
			}
			dzeta = dzeta / ctmp;
		}
		else
		{
			ctmp1 = czero;
			ctmp2 = czero;
			ctmp3 = czero;
			ctmp4 = czero;
			ctmp5 = czero;

			for(ii=jj=0; ii<ndim; ++ii, jj+=nlar)
			{
				ctmp1 += (wk[jj + 6].conjg() * wk[jj + 6]);
				ctmp2 += (wk[jj + 3].conjg() * wk[jj + 3]);
				ctmp3 += (wk[jj + 6].conjg() * wk[jj + 3]);
				ctmp4 += (wk[jj + 6].conjg() * wk[jj + 4]);
				ctmp5 += (wk[jj + 3].conjg() * wk[jj + 4]);
			}

			ctmp = ctmp1 * ctmp2 - ctmp3 * ctmp3.conjg();

			if (ctmp == czero)
			{
				status = -1;
				steperr = 3;
				return;
			}

			dzeta = (ctmp2*ctmp4 - ctmp5*ctmp3) / ctmp;
			eta   = (ctmp1*ctmp5 - ctmp3.conjg()*ctmp4) / ctmp;
		}
//
// Compute u
		for(ii=jj=0; ii<ndim; ++ii, jj+=nlar)
		{
			wk[jj + 7] = dzeta * wk[jj + 5] + eta * (wk[jj + 11] - wk[jj + 2] + beta * wk[jj + 7]);
		}
// 
// Compute z
		for(ii=jj=0; ii<ndim; ++ii, jj+=nlar)
		{
			wk[jj + 9] = dzeta * wk[jj + 2] + eta * wk[jj + 9] - alpha * wk[jj + 7];
		}
//
// Compute x and r
		residu = (real)0.;
		for(ii=jj=0; ii<ndim; ++ii, jj+=nlar)
		{
			wk[jj + 10] = wk[jj + 10] + alpha * wk[jj + 1] + wk[jj + 9];
			wk[jj +  2] = wk[jj +  4] - eta   * wk[jj + 3] - dzeta * wk[jj + 6];
			residu += wk[jj + 2].modSquared();
		}
		residu = Sqrt(residu) / norm;
//
// Flatau add wrimsg call btd modify
		sprintf(cmsgnm, "IT=%8d f.err=%10.3e", itern + 1, residu);
		Wrimsg("Gpbicg", cmsgnm);
//
		if (residu <= tol)
		{
			status = 1;
			for(ii=jj=0; ii<ndim; ++ii, jj+=nlar)
			{
				xi[ii] = wk[jj + 10];
			}
			nou = 4;
			return;
		}
//	do not need break
	case 4:
//
// Compute beta
		ctmp = czero;
		for(ii=jj=0; ii<ndim; ++ii, jj+=nlar)
		{
			ctmp += (wk[jj].conjg() * wk[jj + 2]);
		}
		if (r0rn == czero)
		{
			status = -1;
			steperr = 4;
			return;
		}
		beta = alpha * ctmp / dzeta / r0rn;
		r0rn = ctmp;
//
// Compute w
		for(ii=jj=0; ii<ndim; ++ii, jj+=nlar)
		{
			wk[jj + 8] = wk[jj + 6] + beta * wk[jj + 5];
		}

		if (itern <= maxit) 
			goto l100;

		status = 1;
		steperr = 0;
		for(ii=jj=0; ii<ndim; ++ii, jj+=nlar)
		{
			xi[ii] = wk[jj + 10];
		}
		break;
//
	default:
		status = 0;
		steperr = -1;
//
// Compute norm and Ax0 (x0 initial guess)
		norm = (real)0.;
		for(ii=jj=0; ii<ndim; ++ii, jj+=nlar)
		{
			wk[jj + 10] = xi[ii];
			norm += b[ii].modSquared();
		}
		norm = Sqrt(norm);
		nou = 1;
		break;
	}
}
