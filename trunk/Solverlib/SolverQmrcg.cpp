#include "StdAfx.h"

#include "SolverQmrcg.h"

REGISTER_SOLVER(Qmrcg,QMRCCG)
SolverQmrcg::SolverQmrcg(void)
{
	itern = ncompte = 0;
	ndim = nlar = 0;
	xs = xr = NULL;
	wkFileSize = 10;
}

SolverQmrcg::~SolverQmrcg(void)
{
	CleanDelete2(xs);
	CleanDelete2(xr);
}

void SolverQmrcg::SetParameters(int nndim, int mmaxit)
{
	ndim = nndim;
	nlar = wkFileSize;
	maxit = mmaxit;
	xr = new Complex[ndim];
	xs = new Complex[ndim];
	wk = new Complex[ndim * wkFileSize];
}

void SolverQmrcg::Pimqmrcg(Complex *xi, Complex *b, real tol, real &tole)
{
/* **
Interface to QMR solver
ndim - dimension
matvec - external Ax routine
b - right hand side
lda - leading dimension
nlar - number of needed scratch vectors (9+1)
xi - on input initial guess, on output the result
xr - scratch vector
wrk - work array
maxit - maximum number of iterations
itno - number of iterations
tol - needed relative error
tole - achived  relative error
ncompte - cumber of Ax multiplications

History
(PJF) = Piotr J. Flatau
(PCC) = P. C. Chaumet
(AR) = A. Rahmani
February 4, 2010  P. C. Chaumet and A. Rahmani
May 6, 2010 (PJF) converted to Fortran90, introduce  DDPRECISION to
                   handle single/double precision easily, 
                   introduced pointer/target to split work array
** */
//
// FLATAU split wrk array to 2 arrays which are needed by pimzqmr
// xs => wrk(1:lda,1)
// wrk2 => wrk(1:lda,2:nlar)
//
	char cmsgnm[128];
	nou = 0;
	ncompte = 0;
	int itlast = 1;
	do
	{
		Pimzqmr(xs, xi, b, tole, tol);
		if(itern > itlast)
		{
			sprintf(cmsgnm, "IT=%8d f.err=%10.3lf", itern-2, tole);
			Wrimsg("Gmrccg", cmsgnm);
			itlast = itern;
		}
//
		if (status < 0)
		{
			fprintf(stderr, "stop nstat %d %d\n", status, steperr);
			return;
		}
		++ncompte;

		if (nt == 1)				// A \times x
			Matvec(xi, xr, NULL);
		if (nt == 2)				// A^{T} \times x
			Matvec(xi, xr, NULL);
	} while(status != 1);

	if(steperr == 0)
	{
		sprintf(cmsgnm, "IT=%6d has reached MAXIT=%6d", itern, maxit);
		Errmsg("Fatal", "Pimqmrcg", cmsgnm);
	}
//
// FLATAU after all is done assign xs (solution) to xi
	memcpy(xi, xs, ndim*sizeof(Complex));
}

void SolverQmrcg::Pimzqmr(Complex *xs, Complex *xi, Complex *b, real &tole, real tol)
{
/* ** 
     Iterative solver QMR

     Authors: P. C. Chaumet and A. Rahmani
     Date: 04/02/2010

     Purpose: iterative solver for linear system Ax=b. There is no
     condition on the matrix. Notice that the products A x and At x are provided by the user.

     Reference: if you use this routine in your research, please
     reference, as appropriate: P. C. Chaumet and A. Rahmani, Efficient
     discrete dipole approximation for magneto-dielectric scatterers
     Opt. Lett. 34, 917 (2009). R. D. Da Cunha and T. Hopkins,
     Appl. Numer. Math. 19, 33 (1995).

History
 (PJF) = Piotr J. Flatau
 (PCC) = P. C. Chaumet
 (AR) = A. Rahmani
 Originally written by  R. D. Da Cunha and T. Hopkins
 February 4, 2010 Modified by P. C. Chaumet and A. Rahmani
 May 6, 2010  (PJF) converted to Fortran90, introduce  DDPRECISION to
                    handle sing/double precision easily, introduced pointer/target
                    to split work array
** */
//

	Complex dd, den, ksi0, gamma0;
	int i, ii;

	switch(nou)
	{
	case 0:
//
// 1. lambda=1, kappa=-1, theta=-1
		lambda = coner;
		kappa = -coner;
		theta = -coner;
		norm = (real)0.;
		for(i=0; i<ndim; ++i)
			norm += b[i].modSquared();
		norm = Sqrt(norm);
//
// Loop
		status = 0;
		steperr = -1;
		itern = 0;
//
// 2. wtilde=vtilde=r=b-Ax
// r=b-Ax
// A*x=wrk(3)
		nou = 1;
		nt = 1;
		memcpy(xs, xi, ndim * sizeof(Complex));
//
// Compute A*xi
		break;

	case 1:
		for(i=ii=0; i<ndim; ++i, ii+=nlar)
		{
			wk[ii + 2] = xr[i];
			wk[ii    ] = b[i] - wk[ii + 2];
			wk[ii + 6] = wk[ii];
			wk[ii + 7] = wk[ii];
		}
//
// 3. p=q=d=s=0
		for(i=ii=0; i<ndim; ++i, ii+=nlar)
		{
			wk[ii + 1] = wk[ii + 3] = wk[ii + 4] = wk[ii + 5] = czero;
		}
// 
// 4. gamma=||vtilde||_{2}, ksi=||wtilde||_{2},
// rho=wtilde^{T}vtilde, epsilon=(A^{T}wtilde)^{T}vtilde, mu=0
		dots[0] = dots[1] = dots[2] = czero;
		dd = czero;
		for(i=ii=0; i<ndim; ++i, ii+=nlar)
		{
			dots[0] += (wk[ii + 6] * wk[ii + 6].conjg());
			dots[1] += (wk[ii + 7] * wk[ii + 7].conjg());
			dots[2] += (wk[ii + 6] * wk[ii + 7]);
			dd += (wk[ii + 2] * wk[ii + 2].conjg());
		}
//
// Compute A^{T}wtilde
// CALL TMATVEC(WRK(IWTILDE),WRK(IATWTILDE),IPAR)
		for(i=ii=0; i<ndim; ++i, ii+=nlar)
		{
			xi[i] = wk[ii + 7];
		}
		nt = 2;
		nou = 2;
		break;

	case 2:
		dots[3] = czero;
		for(i=ii=0; i<ndim; ++i, ii+=nlar)
		{
			wk[ii + 8] = xr[i];
			dots[3] += (wk[ii + 6] * wk[ii + 8]);
		}
// 
// Accumulate simultaneously partial inner-products
		gamma = dots[0].sqrt();
		ksi = dots[1].sqrt();
		rho = dots[2];
		epsilon = dots[3];
		mu = czero;
//
// 5. tau=epsilon/rho
		if (rho == czero)
		{
			itern = 0;
			status = -3;
			steperr = 5;
			return;
		}

		tau = epsilon / rho;
		if (!FunAt100(xi))
			return;
		break;

	case 3:
		for(i=ii=0; i<ndim; ++i, ii+=nlar)
		{
			wk[ii + 2] = xr[i];
			wk[ii + 6] = wk[ii + 2] - tau / gamma * wk[ii + 6];
		}
//
// 9. wtilde=q-tau/ksi*wtilde
		for(i=ii=0; i<ndim; ++i, ii+=nlar)
		{
			wk[ii + 7] = wk[ii + 3] - tau / ksi * wk[ii + 7];
		}
//
// 11. gamma=||vtilde||_{2}, ksi=||wtilde||_{2},
// rho=wtilde^{T}vtilde, epsilon=(A^{T}wtilde)^{T}vtilde
		dots[0] = dots[1] = dots[2] = czero;
		for(i=ii=0; i<ndim; ++i, ii+=nlar)
		{
			dots[0] += (wk[ii + 6] * wk[ii + 6].conjg());
			dots[1] += (wk[ii + 7] * wk[ii + 7].conjg());
			dots[2] += (wk[ii + 6] * wk[ii + 7]);
		}
//
// Compute A^{T}wtilde
		for(i=ii=0; i<ndim; ++i, ii+=nlar)
		{
			xi[i] = wk[ii + 7];
		}
        nt = 2;
		nou = 4;
		break;

	case 4:
		dots[3] = czero;
		for(i=ii=0; i<ndim; ++i, ii+=nlar)
		{
			wk[ii + 8] = xr[i];
			dots[3] += (wk[ii + 6] * wk[ii + 8]);
		}
//
// Accumulate simultaneously partial inner-products
		gamma0 = gamma;
		gamma = dots[0].sqrt();
		ksi0 = ksi;
		ksi = dots[1].sqrt();
		{
			Complex rho0 = rho;
			rho = dots[2];
			epsilon = dots[3];
//
// 12. mu=(gamma0*ksi0*rho)/(gamma*tau*rho0)
			den = gamma * tau * rho0;
		}
		if (den == czero)
		{
			status = -3;
			steperr = 12;
			return;
		}
		mu = (gamma0 * ksi0 * rho) / den;
//
// 13. tau=epsilon/rho-gamma*mu
		if (rho == czero)
		{
			status = -3;
			steperr = 13;
            return;
		}
		{
			Complex tau0 = tau;
			tau = epsilon/rho - gamma*mu;
//
// 14. theta=(|tau0|^2*(1-lambda))/(lambda*|tau|^2+|gamma|^2)
			real abstau02 = tau0.modSquared();
			den = lambda*abstau02 + gamma.modSquared();
			if (den == czero)
			{
				status = -3;
				steperr = 14;
				return;
			}
			theta = (coner - lambda) * abstau02 / den;
//
// 15. kappa=(-gamma0*CONJG(tau0)*kappa0)/(gamma0*|tau|^2+|gamma|^2)
			kappa = -(gamma0 * tau0.conjg() * kappa) / den;
//
// 16. lambda=(lambda0*|tau0|^2)/(gamma0*|tau|^2+|gamma|^2)
			lambda = lambda * abstau02 / den;
		}
//
// 17. d=theta*d+kappa*p
		for(i=ii=0; i<ndim; ++i, ii+=nlar)
		{
			wk[ii + 4] = theta * wk[ii + 4] + kappa * wk[ii + 1];
		}
//
// 18. s=theta*s+kappa*A*p
		for(i=ii=0; i<ndim; ++i, ii+=nlar)
		{
			wk[ii + 5] = theta * wk[ii + 5] + kappa * wk[ii + 2];
		}
//
// 19. x=x+d
		for(i=ii=0; i<ndim; ++i, ii+=nlar)
		{
			xs[i] += wk[ii + 4];
		}
//
// 20. r=r-s
		for(i=ii=0; i<ndim; ++i, ii+=nlar)
		{
			wk[ii] -= wk[ii + 5];
		}
//
// Criterion to stop
		{
			real tmp1 = (real)0.;
			for(i=ii=0; i<ndim; ++i, ii+=nlar)
			{
				tmp1 += wk[ii].modSquared();
			}
			tole = Sqrt(abs(tmp1)) / norm;
		}
		if (tole <= tol)
		{
			nou = 5;
			status = 1;
			return;
		}
		if (itern > maxit)
		{
			status = 1;
			steperr = 0;
			return;
		}
		FunAt100(xi);
		break;

	case 5:
		FunAt100(xi);
		break;
	}
}

bool SolverQmrcg::FunAt100(Complex *xi)
{
	int i, ii;
	++itern;
//
// 6. p=1/gamma*vtilde-mu*p
	if (gamma == czero)
	{
		status = -3;
		steperr = 6;
		return false;
	}
	for(i=ii=0; i<ndim; ++i, ii+=nlar)
	{
		wk[ii + 1] = wk[ii + 6] / gamma - mu * wk[ii + 1];
	}
//
// 7. q=1/ksi*A^{T}wtilde-(gamma*mu)/ksi*q
	if (ksi == czero)
	{
		status = -3;
		steperr = 7;
		return false;
	}
	for(i=ii=0; i<ndim; ++i, ii+=nlar)
	{
		wk[ii + 3] = (wk[ii + 8] - gamma * mu * wk[ii + 3]) / ksi;
	}
//
// 8. vtilde=Ap-tau/gamma*vtilde
	for(i=ii=0; i<ndim; ++i, ii+=nlar)
	{
		xi[i] = wk[ii + 1];
	}
	nt = 1;
	nou = 3;

	return true;
}
