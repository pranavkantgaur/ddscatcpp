#ifndef __SOLVERZBCG2_H__
#define __SOLVERZBCG2_H__

#include "AbstractSolver.h"

class SOLVERLIB_API SolverZbcg2 : public AbstractSolver
{
protected:
	int info;
	int l, n, itern, itermx;
	Complex **matrix_z;
	Complex *y0, *yl, *zy0, *zyl;

public:
	SolverZbcg2(void);
	virtual ~SolverZbcg2(void);

public:
	void Zbcg2(Complex *x, Complex *rhs, bool print_resid, bool nonzero_x, real &toler, int &mxmatvec);
	real Dnorm2_bcg(int n, Complex *zx);
	Complex Zdot_bcg(int n, Complex *zx, Complex *zy);
	bool SetParameters(int ll, int nn, int maxiter);

public:
	inline int GetInfo() { return info; }
	inline int GetNumberOfIterations() { return itern; }
};

#endif
