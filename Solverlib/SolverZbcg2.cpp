#include "StdAfx.h"

#include "SolverZbcg2.h"

REGISTER_SOLVER(Zbcg2,PBCGS2)
SolverZbcg2::SolverZbcg2(void)
{
	info = l = n = itern = itermx = 0;
	matrix_z = NULL;
	y0 = yl = zy0 = zyl = NULL;
	wkFileSize = 10;
}

SolverZbcg2::~SolverZbcg2(void)
{
	CleanDelete2(y0);
	CleanDelete2(yl);
	CleanDelete2(zy0);
	CleanDelete2(zyl);
	for(int i=0; i<=l; ++i)
	{
		CleanDelete2(matrix_z[i]);
	}
	CleanDelete2(matrix_z);
}

bool SolverZbcg2::SetParameters(int ll, int nn, int maxiter)
{
	l = ll;
	n = nn;
	itermx = maxiter;
	if ((l < 1) || (l > 2)) 
	{
		info = -2;
		return false;
	}
	matrix_z = new Complex *[l+1];
	for(int i=0; i<=l; ++i)
	{
		matrix_z[i] = new Complex[l+1];
	}
	y0 = new Complex[l+1];
	yl = new Complex[l+1];
	zy0 = new Complex[l+1];
	zyl = new Complex[l+1];
	wk = new Complex[wkFileSize * n];
	return true;
}

void SolverZbcg2::Zbcg2(Complex *x, Complex *rhs, bool print_resid, bool nonzero_x, real &toler, int &mxmatvec)
{
/* **
Improved "vanilla" BiCGStab(2) iterative method

Copyright (c) 2001 by M.A.Botchev (http://www.math.utwente.nl/~botchev/),
                      University of Twente
Permission to copy all or part of this work is granted, provided that the copies are not made or distributed
for resale, and that the copyright notice and this notice are retained.

This is the "vanilla" version of BiCGstab(\ell) as described
in PhD thesis of D.R.Fokkema, Chapter 3 (also available as
Preprint 976, Dept. of Mathematics, Utrecht University, URL
http://www.math.uu.nl/publications/).  It includes two enhancements 
to BiCGstab(\ell) proposed by G.Sleijpen and H.van der Vorst in
1) G.Sleijpen and H.van der Vorst "Maintaining convergence 
   properties of BiCGstab methods in finite precision arithmetic",
   Numerical Algorithms, 10, 1995, pp.203-223
2) G.Sleijpen and H.van der Vorst "Reliable updated residuals in
   hybrid BiCG methods", Computing, 56, 1996, pp.141-163

{{ This code based on original work of D.R.Fokkema:
Subroutine zbistbl v1.1 1998    
Copyright (c) 1995-1998 by D.R. Fokkema.
Permission to copy all or part of this work is granted, provided that the copies are not made or distributed 
for resale, and that the copyright notice and this notice are retained.
}}

Your bug reports, comments, etc. are welcome:  m.a.botchev@math.utwente.nl

------------------------------
Description of the parameters:
------------------------------

x          (input/output) COMPLEX*16 array dimension n initial guess on input, solution on output
rhs        (input) COMPLEX*16 array dimension n the right-hand side (r.h.s.) vector
print_resid (input) LOGICAL. If print_resid=.true. the number of 
           matrix-vector multiplications done so far and residual norm will
           be printed to the standard output each iteration
l          (input) INTEGER the dimension \ell of BiCGstab(\ell)
           in this simple version it is required that l <= 2
           l=2 is often useful for systems with nonsymmetric matrices
n          (input) INTEGER size of the linear system to solve 
matvec     (input) EXTERNAL name of matrix vector subroutine to deliver y:=A*x by CALL matvec(n,x,y)
nonzero_x  (input) LOGICAL tells BiCGstab(\ell) if the initial guess x is zero or not. 
           If nonzero_x is .FALSE., initial residual equals r.h.s. vector and one MATVEC call is avoided
toler      (input/output) DOUBLE PRECISION tolerance: the iterations are 
           stopped as soon as || residual ||/|| initial residual|| <= toler,
           the norm is Euclidean.  On output, if info>=0, the value of 
           toler is set to the actually achieved residual reduction
mxmatvec   (input/output) INTEGER.  On input: maximum number of matrix 
           vector multiplications allowed to be done.  On output: 
           if info>=0, mxmatvec is set to the actual number of matrix 
           vector multiplications done
work       (workspace) COMPLEX*16 array of dimension (n,2*l+5)
info       (output) INTEGER.  info = 0 in case of succesful computations and 
           info = -m (<0) - means paramater number m has an illegal value
           info = 1 - means no convergence achieved (stopping criterion is not fulfilled)
           info = 2 - means breakdown of the algorithm (taking a larger value of parameter l usually helps)

WARNING: If the iterations are ended normally (info=0 or info=1),
the true residual norm is computed and returned as an output value 
of the parameter toler.  The true residual norm can be slightly larger
than the projected residual norm used by the algorithm to stop the
iterations.  It may thus happen that on output info=0 but the value
of toler is (slightly) larger than tolerance prescribed on input.

History:
Fortran versions history removed.
** */
	const real delta = (real)1.e-2;

	info = 0;

	if ((l < 1) || (l > 2)) info = -2;
	if (n < 1) info = -3;
	if (toler <= (real)0.) info = -9;
	if (mxmatvec < 0) info = -10;

	if (info != 0) return;

	int rr = 0;						int n_rr = n * rr;
	int r = rr + 1;					int n_r = n * r;
	int u = r + l + 1;				int n_u = n * u;
	int xp = u + l + 1;				int n_xp = n * xp;
	int bp = xp + 1;				int n_bp = n * bp;

//
// Initialize first residual
	int ih, nmatvec;
	if (nonzero_x) 
	{
		Matvec(x, wk + n_r, NULL);
		for(ih=0; ih<n; ++ih)
			wk[ih + n_r] = rhs[ih] - wk[ih + n_r];
		nmatvec = 1;
	}
	else
	{
		memcpy(wk+n_r, rhs, n*sizeof(Complex));
//		for(ih=0; ih<n; ++ih)
//			wk[ih + n_r] = rhs[ih];
		nmatvec = 0;
	}
//
// Initialize iteration loop
	memcpy(wk+n_rr, wk+n_r, n*sizeof(Complex));
	memcpy(wk+n_bp, wk+n_r, n*sizeof(Complex));
	memcpy(wk+n_xp, x, n*sizeof(Complex));
	memset(x, 0, n*sizeof(Complex));

//	for(ih=0; ih<n; ++ih)
//	{
//		wk[ih + n_rr] = wk[ih + n_r];
//		wk[ih + n_bp] = wk[ih + n_r];
//		wk[ih + n_xp] = x[ih];
//		x[ih] = czero;
//	}

	real rnrm0 = Dnorm2_bcg(n, wk + n_r);
	real rnrm = rnrm0;
//fprintf(stderr, "rnmr0 = %lf\n", rnrm0);

	real mxnrmx = rnrm0;
	real mxnrmr = rnrm0;
	bool rcmp = false;
	bool xpdt = false;

	Complex alpha = czero;
	Complex omega = coner;
	Complex sigma = coner;
	Complex rho0 = coner;

// !BTD 080716:
	itern = 0;
//
// Iterate
	int i, j, k;
	while ((rnrm > toler*rnrm0) && (itern <= itermx))
	{
// !BTD 080716:
		itern++;
		if (itern > itermx)
			Errmsg("Fatal", "Zbcg2wp", " itern > itermx");
//          
// The BiCG part ---
		rho0 = -omega*rho0;
		for(k=0; k<l; ++k)
		{
			Complex rho1 = Zdot_bcg(n, wk + n_rr, wk + (r+k)*n);
// fprintf(stderr, "rho1 = %lf %lf\n", rho1.re, rho1.im);
            if (rho0 == czero)
			{
				info = 2;
				toler = rnrm / rnrm0;
				mxmatvec = nmatvec;
				return;
			}
			Complex beta = alpha * (rho1/rho0);
			rho0 = rho1;
			for(j=0; j<=k; ++j)
			{
				int n_uj = n * (u+j);
				int n_rj = n * (r+j);
				for(ih=0; ih<n; ++ih)
				{
					wk[ih + n_uj] = wk[ih + n_rj] - beta * wk[ih + n_uj];
				}
			}
			Matvec(wk + (u+k)*n, wk + (u+k+1)*n, NULL);
			++nmatvec;

			sigma = Zdot_bcg(n, wk + n_rr, wk + (u+k+1)*n);
// fprintf(stderr, "sigma = %lf %lf\n", sigma.re, sigma.im);
			if (sigma == czero)
			{
				info = 2;
				toler = rnrm / rnrm0;
				mxmatvec = nmatvec;
				return;
			}
			alpha = rho1 / sigma;
			for(ih=0; ih<n; ++ih)
				x[ih] += (alpha * wk[ih + n_u]);
			for(j=0; j<=k; ++j)
			{
				int n_uj1 = n * (u+j+1);
				int n_rj = n * (r+j);
				for(ih=0; ih<n; ++ih)
				{
					wk[ih + n_rj] -= (alpha * wk[ih + n_uj1]);
				}
			}
			Matvec(wk + (r+k)*n, wk + (r+k+1)*n, NULL);
			++nmatvec;
			rnrm = Dnorm2_bcg(n, wk + n_r);
// fprintf(stderr, "rnrm = %lf\n", rnrm);
			mxnrmx = max(mxnrmx, rnrm);
			mxnrmr = max(mxnrmr, rnrm);
		}
//
// The convex polynomial part
// --- Z = R'R
		for(i=0; i<=l; ++i)
			for(j=0; j<=i; ++j)
				matrix_z[i][j] = (Zdot_bcg(n, wk + (r+j)*n, wk + (r+i)*n)).conjg();
//
// Lower triangular part of Z is computed; compute the rest knowing that Z^H=Z
		for(j=1; j<=l; ++j)
			for(ih=0; ih<j; ++ih)
				matrix_z[ih][j] = matrix_z[j][ih].conjg();
//
// Small vectors y0 and yl
		y0[0] = -coner;
		y0[1] = (matrix_z[1][0] / matrix_z[1][1]);								// Works only for l=2
		y0[l] = czero;

        yl[0] =  czero;
        yl[1] = (matrix_z[1][2] / matrix_z[1][1]);								// Works only for l=2
        yl[l] = -coner;
//
// --- Convex combination
//
// Compute Z*y0 and Z*yl
        for(ih=0; ih<=l; ++ih)
		{
			zy0[ih] = czero;
			zyl[ih] = czero;
			for(j=0; j<=l; ++j)
			{
				zy0[ih] += matrix_z[ih][j] * y0[j];
				zyl[ih] += matrix_z[ih][j] * yl[j];
			}
		}

		real kappa0 = Sqrt(Zdot_bcg(l+1, y0, zy0).abs());
		real kappal = Sqrt(Zdot_bcg(l+1, yl, zyl).abs());

        Complex varrho = Zdot_bcg(l+1, yl, zy0) / (kappa0*kappal);

		Complex hatgamma = varrho / varrho.abs() * max_(varrho.mod(), (real)7e-1);
		Complex tmp = hatgamma * kappa0 / kappal;

		for(ih=0; ih<=l; ++ih)
			y0[ih] -= (tmp * yl[ih]);
//
// Update
		omega = y0[l];
		for(j=0; j<l; ++j)
		{
			int n_j = n * j;
			for(ih=0; ih<n; ++ih)
			{
				wk[ih + n_u] -= (y0[j+1] * wk[ih + n_u + n_j + n]);
				x[ih]        += (y0[j+1] * wk[ih + n_r + n_j]);
				wk[ih + n_r] -= (y0[j+1] * wk[ih + n_r + n_j + n]);
			}
		}
//
// y0 has changed; compute Z*y0 once more
		for(ih=0; ih<=l; ++ih)
		{
			zy0[ih] = czero;
			for(j=0; j<=l; ++j)
			{
				zy0[ih] += matrix_z[ih][j] * y0[j];
			}
		}
		rnrm = Sqrt(Zdot_bcg(l+1, y0, zy0).abs());
//
// The reliable update part
		mxnrmx = max(mxnrmx, rnrm);
		mxnrmr = max(mxnrmr, rnrm);
		xpdt =  (rnrm < delta*rnrm0 && rnrm0 < mxnrmx);
		rcmp = ((rnrm < delta*mxnrmr && rnrm0 < mxnrmr) || xpdt);
		if (rcmp) 
		{
			Matvec(x, wk + n_r, NULL);
			++nmatvec;
			for(ih=0; ih<n; ++ih)
				wk[ih + n_r] = wk[ih + n_bp] - wk[ih + n_r];
			mxnrmr = rnrm;
			if (xpdt)
			{
				for(ih=0; ih<n; ++ih)
				{
					wk[ih + n_xp] += x[ih];
//					x[ih] = czero;
//					wk[ih + n_bp] = wk[ih + n_r];
				}
				memset(x, 0, n*sizeof(Complex));
				memcpy(wk+n_bp, wk+n_r, n*sizeof(Complex));
				mxnrmx = rnrm;
			}
		}
//
// IF(print_resid)PRINT *,nmatvec,' ',rnrm
		if (print_resid)
		{
			char cmsgnm[72];
			sprintf(cmsgnm, "IT=%8d  f.err= %10.3e", itern, rnrm/rnrm0);
			Wrimsg("Zbcg2 ", cmsgnm);
		}
	}
//
// End of iterations
	for(ih=0; ih<n; ++ih)
		x[ih] += wk[ih + n_xp];

	if (rnrm > toler*rnrm0)
		info = 1;
//
// Compute the true residual:
// --------------------- One matvec can be saved by commenting out this:
//      CALL MATVEC(X,WORK(1:N,R),N)
// --------------------- One matvec can be saved by commenting out this^
	Matvec(x, wk + n_r, NULL);
	for(ih=0; ih<n; ++ih)
		wk[ih + n_r] = rhs[ih] - wk[ih + n_r];
	rnrm = Dnorm2_bcg(n, wk + n_r);
	++nmatvec;

	toler = rnrm / rnrm0;
	mxmatvec = nmatvec;
}

real SolverZbcg2::Dnorm2_bcg(int n, Complex *zx)											// L2 norm function
{
	return Sqrt(Zdot_bcg(n, zx, zx).mod());
}

Complex SolverZbcg2::Zdot_bcg(int n, Complex *zx, Complex *zy)								// Complex inner product function
{
	Complex result;
	for(int i=0; i<n; ++i)
		result += (zx[i].conjg() * zy[i]);
//	fprintf(stdout, "%d res = %lf %lf\n", n, result.re, result.im);
	return result;
}