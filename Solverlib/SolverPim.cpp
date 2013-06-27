#include "StdAfx.h"
#include "SolverPim.h"

REGISTER_SOLVER(Pim,PETRKP)
REGISTER_SOLVER(Pim,PBCGST)

SolverPim::SolverPim(void)
{
	Init();
}

SolverPim::~SolverPim(void)
{
	if (cashedN > 0)
	{
		free(qi);
		free(gi);
		free(pi);
		free(cr);
		free(axi);
		free(r);
		free(ace);
		cashedN = -1;
	}
}

void SolverPim::Init(void)
{
	memset(spar, 0, PimSparNameEnd * sizeof(real));
	memset(ipar, 0, PimIparNameEnd * sizeof(int));
	Matvec = NULL;
	Cmatvec = NULL;
	Tmatvec = NULL;
	Preconl = NULL;
	Preconr = NULL;
	Pcsum = NULL;
	Pscnrm = NULL;
	Progress = NULL;
	wkFileSize = 10;
	qi = gi = pi = cr = axi = r = ace = NULL;
	cashedN = -1;
}

void SolverPim::SetParameters(int n, int blksz, int loclen, int basisdim, int nprocs, int procid, int precontype, int stoptype, int maxit, real epsilon)
{
	ipar[PimIparNameLda]  = ipar[PimIparNameN]  = n;
	ipar[PimIparNameBlksz]  = blksz;
	ipar[PimIparNameLoclen]  = loclen;
	ipar[PimIparNameBasisdim]  = basisdim;
	ipar[PimIparNameNprocs]  = nprocs;
	ipar[PimIparNameProcid]  = procid;
	ipar[PimIparNamePrecontype]  = precontype;
	ipar[PimIparNameStoptype]  = stoptype;
	ipar[PimIparNameMaxit]  = maxit;
	ipar[PimIparNameItno] = -1;
	ipar[PimIparNameStatus] = -1;
	ipar[PimIparNameSteperr] = -1;

	spar[PimSparNameEpsilon] = epsilon;
	spar[PimSparNameExitnorm] = -(real)1.;
}

void SolverPim::GetParameters(int &n, int &blksz, int &loclen, int &basisdim, int &nprocs, int &procid, int &precontype, int &stoptype, 
	int &maxit, int &itno, int &status, int &steperr, real &epsilon, real &exitnorm)
{
	n = ipar[PimIparNameN];
	blksz = ipar[PimIparNameBlksz];
	loclen = ipar[PimIparNameLoclen];
	basisdim = ipar[PimIparNameBasisdim];
	nprocs = ipar[PimIparNameNprocs];
	procid = ipar[PimIparNameProcid];
	precontype = ipar[PimIparNamePrecontype];
	stoptype = ipar[PimIparNameStoptype];
	maxit = ipar[PimIparNameMaxit];
	itno = ipar[PimIparNameItno];
	status = ipar[PimIparNameStatus];
	steperr = ipar[PimIparNameSteperr];
	epsilon = spar[PimSparNameEpsilon];
	exitnorm = spar[PimSparNameExitnorm];

	wk = new Complex[wkFileSize * n];
}

void SolverPim::GetCbicgstabParameters(int &loclen, int &precontype, int &stoptype, int &maxit, int &itno, int &status, int &steperr, real &epsilon, real &exitnorm)
{									   
	loclen = ipar[PimIparNameLoclen];
	precontype = ipar[PimIparNamePrecontype];
	stoptype = ipar[PimIparNameStoptype];
	maxit = ipar[PimIparNameMaxit];
	itno = ipar[PimIparNameItno];
	status = ipar[PimIparNameStatus];
	steperr = ipar[PimIparNameSteperr];
	epsilon = spar[PimSparNameEpsilon];
	exitnorm = spar[PimSparNameExitnorm];
}

void SolverPim::SetMatvecFunctions(void (*Matvec)(Complex *, Complex *, int *), void (*Cmatvec)(Complex *, Complex *, int *), void (*Tmatvec)(Complex *, Complex *, int *))
{
	this->Matvec = Matvec;
	this->Cmatvec = Cmatvec;
	this->Tmatvec = Tmatvec;
}

void SolverPim::SetPreFunctions(void (*Preconl)(Complex *, Complex *, int *), void (*Preconr)(Complex *, Complex *, int *))
{
	this->Preconl = Preconl;
	this->Preconr = Preconr;
}

void SolverPim::SetNormFunctions(void (*Pcsum)(int, Complex *), real (*Pscnrm)(int, Complex *))
{
	this->Pcsum = Pcsum;
	this->Pscnrm = Pscnrm;
}

void SolverPim::SetProgressFunction(void (*Progress)(int, int, real, Complex *, Complex *, Complex *))
{
	this->Progress = Progress;
}

bool SolverPim::Petr(Complex *x, Complex *b)
{
/* **
Solves Ax=b

Input:
	B(LOCLEN) - RHS (doesn't change on output)
	WRK(LDA, LOCLEN) - scratch array
	IPAR - defines int paramteres (see documentation of PIM)
	SPAR - define real parameters
	LDA - leading dimension of arrays
	MATVEC -- name of procedure for computing y=Ax
	CMATVEC -- name of procedure for computing y= conj(A') x

Output:
	X(LOCLEN)

Reference:
   Petravic,M and Kuo-Petravic,G. 1979, JCP, 32,263

History:
Fortran versions history removed.
end history
** */
//
	if (!Matvec || !Cmatvec)					// these functions are necessary for Petr
		return false;
	if (!Pscnrm)
		return false;
	if (!Progress)
		return false;

	int loclen = ipar[PimIparNameLoclen];
	int maxit = ipar[PimIparNameMaxit];
	real epsilon = spar[PimSparNameEpsilon];
//
	const int iace = 0;
	const int igi  = iace + loclen;
	const int ipi  = igi + loclen;
	const int iqi  = ipi + loclen;
	const int iaxi = iqi + loclen;
	const int ir   = iaxi + loclen;
//
// Compute RHSSTOP=EPSILON*|B|
	real rhsstop = Scsetrhsstop(b, NULL, epsilon);
//
// Compute conjg(A')*B
	int l;
	Cmatvec(b, wk + iace, NULL);
	real bnorm = (real)0.;
	for(l=0; l<loclen; ++l)
	{
		bnorm += b[l].modSquared();
		wk[igi + l] = wk[iace + l];
		wk[ipi + l] = wk[igi + l];
	}
//
// Compute |QI>=A|PI>
	Matvec(wk + ipi, wk + iqi, NULL);
	real qiqi = (real)0.;
	real gigi = (real)0.;
	for(l=0; l<loclen; ++l)
	{
		gigi += wk[igi + l].modSquared();
		qiqi += wk[iqi + l].modSquared();
	}
	real alphai = gigi / qiqi;
	for(l=0; l<loclen; ++l)
		x[l] += wk[ipi + l] * alphai;
//
// Compute |AX1>:
	Matvec(x, wk + iaxi, NULL);

	for(int itno=1; itno<=maxit; ++itno)
	{
		ipar[PimIparNameItno] = itno;
//
// Transfer <GI|GI> -> <GI-1|GI-1>:
		real gi1gi1 = gigi;
//
// Compute |GI>=AC|E>-AC|AXI>:
		Cmatvec(wk + iaxi, wk + igi, NULL);
//
// Compute GIGI=<GI|GI>:
		gigi = (real)0.;
		for(l=0; l<loclen; ++l)
			wk[igi + l] = wk[iace + l] - wk[igi + l];
		for(l=0; l<loclen; ++l)
			gigi += wk[igi + l].modSquared();
//
// Compute BETAI-1=<GI|GI>/<GI-1|GI-1>:
		real betai1 = gigi / gi1gi1;
//
// Compute |PI>=|GI>+BETAI-1|PI-1>:
		for(l=0; l<loclen; ++l)
			wk[ipi + l] = wk[igi + l] + wk[ipi + l] * betai1;
//
// Compute |QI>=A|PI>:
		Matvec(wk + ipi, wk + iqi, NULL);
//
// Compute <QI|QI>:
		qiqi = (real)0.;
		for(l=0; l<loclen; ++l)
			qiqi += wk[iqi + l].modSquared();
//
// Compute ALPHAI=<GI|GI>/<QI|QI>:
		alphai = gigi / qiqi;
//
// Compute |XI+1>=|XI>+ALPHAI*|PI>:
		for(l=0; l<loclen; ++l)
			x[l] += (wk[ipi + l] * alphai);
// 
// Except every 10TH iteration, compute |AXI+1>=|AXI>+ALPHAI*|QI> Draine: warning this part perhaps not checked
		if (itno % 10 == 0)
		{
			for(l=0; l<loclen; ++l)
				wk[iaxi + l] += (wk[iqi + l] * alphai);
		}
		else
            Matvec(x, wk + iaxi, NULL);
//
// Compute residual vector |RI>=A|XI>-|B>
		for(l=0; l<loclen; ++l)
			wk[ir + l] = wk[iaxi + l] - b[l];
//
// Call STOPCRIT to check whether to terminate iteration.
// Criterion: terminate if <RI|RI> < TOL**2 * <B|B>
// Return with STATUS=0 if this condition is satisfied.
		real exitnorm;
		int status = Stopcrit(b, wk + ir, NULL, NULL, NULL, wk, rhsstop, 0, exitnorm);
//
// Call PROGRESS to report "error" = sqrt(<RI|RI>)/sqrt(<B|B>)
		Progress(loclen, itno, exitnorm, NULL, NULL, NULL);

		if (status == 0)
			break;
		fprintf(stderr, "Petr %d %lf %lf\n", itno, rhsstop, exitnorm);
	}
	ipar[PimIparNameStatus] = 0;

	return true;
}

/* **
                             CCGPAK 2.0
 Conjugate gradient package for solving complex matrix equations (Fortran90)
                          Piotr J. Flatau

                       last change August 9, 2012
If you use this library in publication please reference:
Flatau, P. J., 2012, CCGPAC 2.0 - Conjugate gradient package for solving complex matrix equations, google code
http://code.google.com/p/conjugate-gradient-lib/

Copyright (C) 2012 P.J. Flatau
Copyright (c) 2013, C++ version, V.Choliy
This code is covered by the GNU General Public License.
 ----------------------------------------------------------------
library of CG implementations
 CGSQR     - Complex conjugate gradient-squared algorithm of Sonneveld,
 CGSTAB    - Complex conjugate gradient squared stabilized
 PETR      - Petravic,M and Kuo-Petravic,G. 1979, JCP, 32,263
 CORS      - The BiCOR and CORS Iterative Algorithms 
** */

bool SolverPim::Petr90(Complex *x, Complex *b, int n, Complex *cxsc)
{
// solves Ax=b
//
// Input:
// B(n) - RHS (doesn't change on output)
// MATVEC -- name of procedure for computing y=Ax
// CMATVEC -- name of procedure for computing y= conj(A') x
//
// Output:
// X(n)
//
// Reference:
//    Petravic,M and Kuo-Petravic,G. 1979, JCP, 32,263
// History records of Fortran version removed.
// =======================================================================
	
	const int arraySize = n * sizeof(Complex);
	int i;
	if (cxsc)
	{
		qi  = cxsc;
		gi  = qi + n;
		pi  = gi + n;
		cr  = pi + n;
		axi = cr + n;
		r   = axi + n;
		ace = r + n;
		cashedN = 0;
	}
	else
	{
		if (n != cashedN)
		{
			qi  = (Complex *)realloc(qi,  arraySize);
			gi  = (Complex *)realloc(gi,  arraySize);
			pi  = (Complex *)realloc(pi,  arraySize);
			cr  = (Complex *)realloc(cr,  arraySize);
			axi = (Complex *)realloc(axi, arraySize);
			r   = (Complex *)realloc(r,   arraySize);
			ace = (Complex *)realloc(ace, arraySize);
			if (!qi || !gi || !pi || !cr || !axi || !r || !ace)
			{
				fprintf(stderr, "Allocation Error Detected in conjugate gradient petr90");
				return false;
			}
			cashedN = n;
		}
	}
// Compute conjg(A')*B
	Cmatvec(b, ace, &n);
	Complex bnorm = DotProduct(b, b, n);
	memcpy(gi, ace, arraySize);
	memcpy(pi, gi, arraySize);
// Compute |QI>=A|PI>
	Matvec(pi, qi, &n);
	Complex gigi = DotProduct(gi, gi, n);
	Complex qiqi = DotProduct(qi, qi, n);
	Complex alphai = gigi / qiqi;
	for(i=0; i<n; ++i)
	{
		x[i] += +alphai * pi[i];
	}
// compute |AX1>:
	Matvec(x, axi, &n);
	int itnoMax = 0;
	for(int itno=0; itno<cgstruct.maxit; ++itno)
	{
		itnoMax = itno;
		Complex gi1gi1 = gigi;								// Transfer <GI|GI> -> <GI-1|GI-1>:
		Cmatvec(axi, gi, &n);								// Compute |GI>=AC|E>-AC|AXI>:
		for(i=0; i<n; ++i)
		{
			gi[i] = ace[i] - gi[i];							// Compute GIGI=<GI|GI>:
			gigi = DotProduct(gi, gi, n);
		}
		Complex betai1 = gigi / gi1gi1;								// Compute BETAI-1=<GI|GI>/<GI-1|GI-1>:
		for(i=0; i<n; ++i)
		{
			pi[i] = gi[i] + betai1 * pi[i];					// Compute |PI>=|GI>+BETAI-1|PI-1>:
		}
		Matvec(pi, qi, &n);									// Compute |QI>=A|PI>:
		qiqi = DotProduct(qi, qi, n);						// Compute <QI|QI>:
		alphai = gigi / qiqi;								// Compute ALPHAI=<GI|GI>/<QI|QI>:
		for(i=0; i<n; ++i)
		{
			x[i] += alphai * pi[i];							// Compute |XI+1>=|XI>+ALPHAI*|PI>:
		}
// Except every 10TH iteration, compute |AXI+1>=|AXI>+ALPHAI*|QI>
// Draine: warning this part perhaps not checked
		if (itno % 10)
		{
			for(i=0; i<n; ++i)
			{
				axi[i] += alphai * qi[i];
			}
		}
		else
		{
			Matvec(x, axi, &n);
		}
// Compute residual vector |RI>=A|XI>-|B>
		for(i=0; i<n; ++i)
		{
			r[i] = axi[i] - b[i];
		}
		Complex rnorm = DotProduct(r, r, n);
		Complex tmp = (rnorm / bnorm).sqrt();
		if (cgstruct.print)
		{
			if(cgstruct.print)
			{
				fprintf(cgstruct.ioerr, "itno = %d, sqrt(rnorm/bnorm)= %lf %lf\n", itno, tmp.re, tmp.im);
			}
		}
		if (tmp.mod() < cgstruct.epsilon_err)
			break;
	}
//
	cgstruct.itno = itnoMax;
	return true;
}

int SolverPim::Stopcrit(Complex *b, Complex *r, Complex *rtrue, Complex *x, Complex *xold, Complex *wrk, real rhsstop, int cnvrtx, real &exitnorm)
{
	int loclen = ipar[PimIparNameLoclen];
	int precontype = ipar[PimIparNamePrecontype];
	int stoptype = ipar[PimIparNameStoptype];

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
			return ((exitnorm < rhsstop) ? 0 : -99);
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
			}
			return ((exitnorm < rhsstop) ? 0 : -99);
		}
		break;

	case 4:
	case 5:
	case 6:										// ||Q1r|| < epsilon or ||Q1r|| < epsilon||b|| or ||Q1r|| < epsilon||Q1b
		exitnorm = Pscnrm(loclen, r);
		return ((exitnorm < rhsstop) ? 0 : -99);
		break;

	case 7:										// ||x-x0||<epsilon
		Ccopy(loclen, x, 1, wrk, 1);
		Caxpy(loclen, -coner, xold, 1, wrk, 1);
        exitnorm = Pscnrm(loclen, wrk);
		return ((exitnorm < rhsstop) ? 0 : -99);
		break;

	default:
		return 0;
		break;
	}
}

real SolverPim::Scsetrhsstop(Complex *b, Complex *wrk, real epsilon)
{	
	int loclen = ipar[PimIparNameLoclen];
	int stoptype = ipar[PimIparNameStoptype];
	
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

bool SolverPim::PimCbicgstab(Complex *x, Complex *b)
{
	if (!Matvec)
		return false;
	if (!Preconl)
		return false;
	if (!Pcsum || !Pscnrm)
		return false;

	real exitnorm;
	int precontype, loclen, iz, itno, status, steperr;
	do
	{
		Complex dots[2];
		real epsilon;
		int stoptype, maxit;
		real macheps = Smachcons('m');
		GetCbicgstabParameters(loclen, precontype, stoptype, maxit, itno, status, steperr, epsilon, exitnorm);
//
// Check consistency of preconditioning and stop types
		if (((precontype == 0) || (precontype == 2)) && (stoptype == 6))
		{
			itno = 0;
			status = -4;
			steperr = 0;
			break;
		}
//
// Does not need conversion Y=Q2X for residual
		int cnvrtx = 0;
//
// Set indices for mapping local vectors into wrk
		int ir = 0;
		int irtilde = ir + loclen;
		int ip = irtilde + loclen;
		int iq = ip + loclen;
		int is = iq + loclen;
		int it = is + loclen;
		int iv = it + loclen;
		int iw = iv + loclen;
		iz = iw + loclen;
		int ixold = iz + loclen;
//
// Set rhs of stopping criteria
		real rhsstop = Scsetrhsstop(b, wk+ir, epsilon);
//
// 1. r=Q1(b-AQ2x)
		if (stoptype != 6)
		{
			switch(precontype)
			{
			case 0:													// r=b-Ax
				Ccopy(loclen, b, 1, wk+ir, 1);
				Matvec(x, wk+iw, NULL);
				Caxpy(loclen, -coner, wk+iw, 1, wk+ir, 1);
				break;

			case 1:													// r=Q1(b-Ax)
				Ccopy(loclen, b, 1, wk+iz, 1);
				Matvec(x, wk+iw, NULL);
				Caxpy(loclen, -coner, wk+iw, 1, wk+iz, 1);
				Preconl(wk+iz, wk+ir, NULL);
				break;

			case 2:													// r=b-AQ2x
				Ccopy(loclen, b, 1, wk+ir, 1);
				Preconr(x, wk+iw, NULL);
				Matvec(wk+iw, wk+iz, NULL);
				Caxpy(loclen, -coner, wk+iz, 1, wk+ir, 1);
				break;

			case 3:													// r=Q1(b-AQ2x)
				Ccopy(loclen, b, 1, wk+ip, 1);
				Preconr(x, wk+iw, NULL);
				Matvec(wk+iw, wk+iz, NULL);
				Caxpy(loclen, -coner, wk+iz, 1, wk+ip, 1);
				Preconl(wk+ip, wk+ir, NULL);
				break;
			}
		}
		else														// r has been set to Qb in the call to dsetrhsstop
		{
			switch(precontype)
			{
			case 1:													// r=Q1(b-Ax)
				Matvec(x, wk+iw, NULL);
				Preconl(wk+iw, wk+iz, NULL);
				Caxpy(loclen, -coner, wk+iz, 1, wk+ir, 1);
				break;

			case 3:													// r=Q1(b-AQ2x)
				Preconr(x, wk+iz, NULL);
				Matvec(wk+iz, wk+iw, NULL);
				Preconl(wk+iw, wk+iz, NULL);
				Caxpy(loclen, -coner, wk+iz, 1, wk+ir, 1);
				break;
			}
		}
//
// 2. rtilde=r
		Ccopy(loclen, wk+ir, 1, wk+irtilde, 1);
//
// 3. p=v=0
		Cinit(loclen, czero, wk+ip, 1);
		Cinit(loclen, czero, wk+iv, 1);
//
// 4. rho=alpha=omega=1
		Complex rho = coner;
		Complex alpha = coner;
		Complex omega = coner;
//
// Loop
		status = 0;
		exitnorm = -(real)1.;
		steperr = -1;
		bool bOk = true;
		for(itno=0; itno<maxit; ++itno)
		{
//
// 5. rho=dot(rtilde,r)
			Complex rho0 = rho;
			dots[0] = Cdotc(loclen, wk+irtilde, 1, wk+ir, 1);
			Pcsum(1, dots);
			rho = dots[0];
//
// 6. beta=rho*alpha/(rho0*omega)
			Complex kappa = rho0 * omega;
			if (kappa == czero)
			{
				status = -3;
				steperr = 6;
				bOk = false;
				break;
			}
			Complex beta = rho * alpha / kappa;
//
// 7. p=r+beta*(p-omega*v)
			Caxpy(loclen, -omega, wk+iv, 1, wk+ip, 1);
			Ccopy(loclen, wk+ip, 1, wk+iw, 1);
			Ccopy(loclen, wk+ir, 1, wk+ip, 1);
			Caxpy(loclen, beta, wk+iw, 1, wk+ip, 1);
//
// 8. v=Q1AQ2p
			switch(precontype)
			{
			case 0:
				Matvec(wk+ip, wk+iv, NULL);
				break;

			case 1:
				Matvec(wk+ip, wk+iz, NULL);
				Preconl(wk+iz, wk+iv, NULL);
				break;

			case 2:
				Preconr(wk+ip, wk+iz, NULL);
				Matvec(wk+iz, wk+iv, NULL);
				break;

			case 3:
				Preconr(wk+ip, wk+iv, NULL);
				Matvec(wk+iv, wk+iz, NULL);
				Preconl(wk+iz, wk+iv, NULL);
				break;
			}
//
// 9. xi=dot(rtilde,v)
			dots[0] = Cdotc(loclen, wk+irtilde, 1, wk+iv, 1);
			Pcsum(1, dots);
			Complex xi = dots[0];
//
// 10. alpha=rho/xi
			if (xi == czero)
			{
				status = -3;
				steperr = 10;
				bOk = false;
				break;
			}
			alpha = rho/xi;
//
// 11. s=r-alpha*v
			Ccopy(loclen, wk+ir, 1, wk+is, 1);
			Caxpy(loclen, -alpha, wk+iv, 1, wk+is, 1);
//
// 12. if ||s||<breaktol then soft-breakdown has occurred
			kappa = Complex(Pscnrm(loclen, wk+is), (real)0.);
			if (kappa.mod() < macheps)
			{
				status = -2;
				steperr = 12;
				Caxpy(loclen, alpha, wk+ip, 1, x, 1);
				bOk = false;
				break;
			}
//
// 13. t=Q1AQ2s
			switch(precontype)
			{
			case 0:
                Matvec(wk+is, wk+it, NULL);
				break;

			case 1:
				Matvec(wk+is, wk+iz, NULL);
				Preconl(wk+iz, wk+it, NULL);
				break;

			case 2:
				Preconr(wk+is, wk+iz, NULL);
				Matvec(wk+iz, wk+it, NULL);
				break;

			case 3:
				Preconr(wk+is, wk+it, NULL);
				Matvec(wk+it, wk+iz, NULL);
				Preconl(wk+iz, wk+it, NULL);
				break;
			}
//
// 14. omega=dot(t,s)/dot(t,t)
			dots[0] = Cdotc(loclen, wk+it, 1, wk+it, 1);
			dots[1] = Cdotc(loclen, wk+it, 1, wk+is, 1);
//
// Accumulate simultaneously partial values
			Pcsum(2, dots);
			if (dots[0] == czero)
			{
				status = -3;
				steperr = 14;
				bOk = false;
				break;
			}
			omega = dots[1] / dots[0];
//
// 15. x=x+alpha*p+omega*s
			Ccopy(loclen, x, 1, wk+ixold, 1);
			Caxpy(loclen, alpha, wk+ip, 1, x, 1);
			Caxpy(loclen, omega, wk+is, 1, x, 1);
//
// 16. r=s-omega*t
			Ccopy(loclen, wk+is, 1, wk+ir, 1);
			Caxpy(loclen, -omega, wk+it, 1, wk+ir, 1);
//
// 17. check stopping criterion
			status = Stopcrit(b, wk+ir, wk+iz, x, wk+ixold, wk+iw, rhsstop, cnvrtx, exitnorm);
//
// Call monitoring routine
			Progress(loclen, itno, exitnorm, x, wk+ir, wk+iz);
			if (status == 0)
			{
				bOk = false;
				break;
			}
		}
		if (itno > maxit)
		{
			status = -1;
			itno = maxit;
		}
		if (bOk == false)
			break;
	} while(0);

	if ((precontype == 2) || (precontype == 3))
	{
		Ccopy(loclen, x, 1, wk+iz, 1);
		Preconr(wk+iz, x, NULL);
	}
//
// Set output parameters
	ipar[PimIparNameItno] = itno;
	ipar[PimIparNameStatus] = status;
	ipar[PimIparNameSteperr] = steperr;
	spar[PimSparNameExitnorm]= exitnorm;

	return true;
}

real SolverPim::Smachcons(char what)
{
	const real macheps   = (real)1.19209e-07;				// These values are for IEEE-754 arithmetic
	const real underflow = (real)0.11754945e-37;
	const real overflow  = (real)1.7014118e38;

	real result = (real)0.;
	switch(what)
	{
	case 'M':
	case 'm':
		result = macheps;
		break;

	case 'U':
	case 'u':
		result = underflow;
		break;

	case 'O':
	case 'o':
		result = overflow;
		break;

	default:
		result = -(real)1.;
		break;
	}

	return result;
}