#ifndef __ABSTRACTSOLVER_H__
#define __ABSTRACTSOLVER_H__

#include "Solverlib.h"
#include "Enumerator.h"
#include "Functions.h"
#include "Complex.h"
#include "Blas.h"

typedef enum _PimIparName
{
	PimIparNameLda, PimIparNameN, PimIparNameBlksz, 
	PimIparNameLoclen, PimIparNameBasisdim, PimIparNameNprocs, 
	PimIparNameProcid, PimIparNamePrecontype, PimIparNameStoptype, 
	PimIparNameMaxit, PimIparNameItno, PimIparNameStatus, 
	PimIparNameSteperr, PimIparNameEnd
} PimIparName;

typedef enum _PimSparName
{
	PimSparNameEpsilon, PimSparNameExitnorm, PimSparNameSmallEigen, PimSparNameLargeEigen, PimSparNameSmallImaginary, 
	PimSparNameLargeImaginary, PimSparNameEnd
} PimSparName;

#define REGISTER_SOLVER(x,y) \
AbstractSolver *Create##y() \
{ \
	AbstractSolver *result = new Solver##x; \
	return result; \
} \
const bool reg##y = AbstractSolver::RegisterSolver(SolMethod_##y, Create##y);

extern const Complex coner;
extern const Complex czero;

class AbstractSolver;

typedef AbstractSolver *(*SolverCreator)();

class SOLVERLIB_API AbstractSolver
{
protected:
	static SolMethod method;
	static AbstractSolver *solver;
	AbstractSolver(void);
	virtual ~AbstractSolver(void);
	void(*Matvec)(Complex *, Complex *, int *);

protected:
	int wkFileSize;
	Complex *wk;
	static real errscal;
	static map<SolMethod, SolverCreator> *factory;

public:
	static AbstractSolver *GetInstance(SolMethod method);
	static void Kill(void);

public:
	static void DeleteSolver();
	static bool RegisterSolver(SolMethod method, SolverCreator creator);
	static real &Errscal() { return errscal; }
	void SetMatvecFunction(void (*Matvec)(Complex *, Complex *, int *));
	Complex DotProduct(Complex *a, Complex *b, int n);
};

#endif
