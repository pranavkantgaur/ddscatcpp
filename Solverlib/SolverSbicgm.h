#ifndef __SOLVERSBICGM_H__
#define __SOLVERSBICGM_H__

#include "AbstractSolver.h"
#include "CgStruct.h"

class SOLVERLIB_API SolverSbicgm : public AbstractSolver
{
protected:
	Complex *cap, *cp, *cr;
	int cashedN;
	CgStruct cgstruct;
	int nlar;

protected:
	void (*Cmatvec)(Complex *, Complex *, int *);
	void (*Tmatvec)(Complex *, Complex *, int *); 

public:
	SolverSbicgm(void);
	virtual ~SolverSbicgm(void);

public:
	inline CgStruct &GetCgStruct() { return cgstruct; }
	inline int GetNumberOfIterations() { return cgstruct.itno; }

public:
	bool Sbicg90(Complex *xi, Complex *b, int n, Complex *cxsc = NULL);
	void SetMatvecFunctions(void (*Matvec)(Complex *, Complex *, int *), void (*Cmatvec)(Complex *, Complex *, int *), void (*Tmatvec)(Complex *, Complex *, int *));
};

#endif
