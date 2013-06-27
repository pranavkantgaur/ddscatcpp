#ifndef __SOLVERTANGCG_H__
#define __SOLVERTANGCG_H__

#include "AbstractSolver.h"

class SOLVERLIB_API SolverTangcg : public AbstractSolver
{
protected:
	int ncompte, itern, ndim, maxit, nou, steperr, status;
	Complex *xr;
	Complex alpha, beta, dzeta, eta, r0rn;
	real norm, residu;

public:
	SolverTangcg(void);
	virtual ~SolverTangcg(void);

public:
	void Tangcg(Complex *xi, Complex *b, real tol, real &tole);
	void SetParameters(int nndim, int mmaxit);

public:
	inline int GetNumberOfMultiplications() { return ncompte; }
	inline int GetNumberOfIterations() { return itern; }

protected:
	void Gpbicg(Complex *xi, Complex *b, real tol);
};

#endif
