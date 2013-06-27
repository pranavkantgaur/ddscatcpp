#ifndef __SOLVERQMRCG_H__
#define __SOLVERQMRCG_H__

#include "AbstractSolver.h"

class SOLVERLIB_API SolverQmrcg : public AbstractSolver
{
protected:
	int ncompte, itern, status, steperr;
	int ndim, nlar, maxit, nou, nt;
	real norm;
	Complex *xs, *xr, dots[4];
	Complex lambda, kappa, theta, gamma, ksi, rho, epsilon, mu, tau;

public:
	SolverQmrcg(void);
	virtual ~SolverQmrcg(void);

public:
	inline int GetNumberOfMultiplications() { return ncompte; }
	inline int GetNumberOfIterations() { return itern; }

public:
	void Pimqmrcg(Complex *xi, Complex *b, real tol, real &tole);
	void SetParameters(int nndim, int mmaxit);

protected:
	void Pimzqmr(Complex *xs, Complex *xi, Complex *b, real &tole, real tol);
	bool FunAt100(Complex *xi);
};

#endif
