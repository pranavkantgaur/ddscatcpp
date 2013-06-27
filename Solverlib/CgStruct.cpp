#include "StdAfx.h"

#include "CgStruct.h"

CgStruct::CgStruct(void)
{
	ioerr = stderr;
	n = blksz = loclen = basisdim = nprocs = procid = maxit = itno = steperr = 0;
	precontype = Precontype_End;
	stoptype = Stoptype_End;
	status = StopStatus_End;
	epsilon_err = exitnorm = mineig_real = maxeig_real = mineig_imag = maxeig_imag = (real)0.;
	print = false;
}

CgStruct::~CgStruct(void)
{
}

void CgStruct::QuickParam(int n, Stoptype stoptype, Precontype precontype, real tol, bool bPrint, int maxiter, int basisdim)
{
	this->n = n;
	this->ioerr = stderr;
	this->stoptype = stoptype;
	this->precontype = precontype;	
	this->epsilon_err = tol;
	this->print = bPrint;
	this->maxit = maxiter;
	this->basisdim = basisdim;
}