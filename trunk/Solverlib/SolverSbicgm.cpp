#include "StdAfx.h"

#include "SolverSbicgm.h"

REGISTER_SOLVER(Sbicgm,SBICGM)
SolverSbicgm::SolverSbicgm(void)
{
	cap = cp = cr = NULL;
	nlar = 0;
	wkFileSize = 3;
}

SolverSbicgm::~SolverSbicgm(void)
{
	if (cashedN > 0)
	{
		free(cap);
		free(cp);
		free(cr);
		cashedN = -1;
	}
}

//
// Purpose: SBICG --- simplified biconjugate gradient method specialized to symmetric matrices
bool SolverSbicgm::Sbicg90(Complex *cx, Complex *cy, int n, Complex *cxsc)
{
	const int arraySize = n * sizeof(Complex);
	int i;
	if (cxsc)
	{
		cap = cxsc;
		cp  = cap + n;
		cr  = cp + n;
		cashedN = 0;
	}
	else
	{
		if (n != cashedN)
		{
			cap = (Complex *)realloc(cap, arraySize);
			cp  = (Complex *)realloc(cp,  arraySize);
			cr  = (Complex *)realloc(cr,  arraySize);
			if (!cap || !cp || !cr)
			{
				fprintf(stderr, "Allocation Error Detected in conjugate gradient sbicg90");
				return false;
			}
			cashedN = n;
		}
	}
//
	Complex ay = DotProduct(cy, cy, n);
	Matvec(cx, cr, &n);
	for(i=0; i<n; ++i)
	{
		cr[i] = cy[i] - cr[i];
		cp[i] = cr[i];
	}
	Complex csk;
	for(i=0; i<n; ++i)
	{
		csk += cr[i] * cr[i];
	}
//
	for(int iter=0; iter<cgstruct.maxit; ++iter)
	{
		Matvec(cp, cap, &n);
		Complex cak;
		for(i=0; i<n; ++i)
		{
			cak += cap[i] * cp[i];
		}
		cak = csk / cak;
		for(i=0; i<n; ++i)
		{
			cx[i] += cak * cp[i];
			cr[i] -= cak * cap[i];
		}
		Complex csk2;
		for(i=0; i<n; ++i)
		{
			csk2 += cr[i] * cr[i];
		}
		Complex ek = DotProduct(cr, cr, n);
		Complex tmp = (ek/ay).sqrt();
		if (cgstruct.print)
		{
			fprintf(cgstruct.ioerr, "iter %d sqrt(ek/ay)= %lf %lf\n", iter, tmp.re, tmp.im);
		}
		cgstruct.itno = iter;
		if (tmp.mod() < cgstruct.epsilon_err) 
			break;
		Complex cbk = csk2 / csk;
		for(i=0; i<n; ++i)
		{
			cp[i] = cr[i] + cbk * cp[i];
		}
		csk = csk2;
	}
//
	return true;
}

void SolverSbicgm::SetMatvecFunctions(void (*Matvec)(Complex *, Complex *, int *), void (*Cmatvec)(Complex *, Complex *, int *), void (*Tmatvec)(Complex *, Complex *, int *))
{
	this->Matvec = Matvec;
	this->Cmatvec = Cmatvec;
	this->Tmatvec = Tmatvec;
}