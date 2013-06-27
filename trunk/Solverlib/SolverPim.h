#ifndef __SOLVERPIM_H__
#define __SOLVERPIM_H__

#include "AbstractSolver.h"
#include "CgStruct.h"

class SOLVERLIB_API SolverPim : public AbstractSolver
{
protected:
	real spar[PimSparNameEnd];
	int ipar[PimIparNameEnd];
	CgStruct cgstruct;
	Complex *qi, *gi, *pi, *cr, *axi, *r, *ace;				// Petr90 need them
	int cashedN;											// Petr90 need it

protected:
	void (*Cmatvec)(Complex *, Complex *, int *);
	void (*Tmatvec)(Complex *, Complex *, int *); 
	void (*Preconl)(Complex *, Complex *, int *);
	void (*Preconr)(Complex *, Complex *, int *); 
	void (*Pcsum)(int, Complex *);
	real (*Pscnrm)(int, Complex *);
	void (*Progress)(int, int, real, Complex *, Complex *, Complex *);

public:
	SolverPim(void);
	virtual ~SolverPim(void);

public:
	bool Petr(Complex *x, Complex *b);
	bool Petr90(Complex *x, Complex *b, int n, Complex *cxsc = NULL);
	bool PimCbicgstab(Complex *x, Complex *b);

public:
	inline int GetIntParameter(PimIparName index) { return ipar[(int)index]; }
	inline void SetIntParameter(PimIparName index, int value) { ipar[(int)index] = value; }
	inline real GetRealParameter(PimSparName index) { return spar[(int)index]; }
	inline int GetNumberOfIterations() { return ipar[PimIparNameItno]; }
	inline CgStruct &GetCgStruct() { return cgstruct; }

public:
	void SetParameters(int n, int blksz, int loclen, int basisdim, int nprocs, int procid, int precontype, int stoptype, int maxit, real epsilon);
	void GetParameters(int &n, int &blksz, int &loclen, int &basisdim, int &nprocs, int &procid, int &precontype, int &stoptype, int &maxit, 
		int &itno, int &status, int &steperr, real &epsilon, real &exitnorm);
	void GetCbicgstabParameters(int &loclen, int &precontype, int &stoptype, int &maxit, int &itno, int &status, int &steperr, real &epsilon, real &exitnorm);
	void SetMatvecFunctions(void (*Matvec)(Complex *, Complex *, int *), void (*Cmatvec)(Complex *, Complex *, int *), void (*Tmatvec)(Complex *, Complex *, int *));
	void SetPreFunctions(void (*Preconl)(Complex *, Complex *, int *), void (*Preconr)(Complex *, Complex *, int *)); 
	void SetNormFunctions(void (*Pcsum)(int, Complex *), real (*Pscnrm)(int, Complex *));
	void SetProgressFunction(void (*Progress)(int, int, real, Complex *, Complex *, Complex *));
	int Stopcrit(Complex *b, Complex *r, Complex *rtrue, Complex *x, Complex *xold, Complex *wrk, real rhsstop, int cnvrtx, real &exitnorm);
	real Scsetrhsstop(Complex *b, Complex *wrk, real epsilon);

protected:
	void Init(void);
	real Smachcons(char what);
};

#endif

/* **
           PIM -- The Parallel Iterative Methods package
           ---------------------------------------------

                      Rudnei Dias da Cunha
               Centro de Processamento de Dados,
         Universidade Federal do Rio Grande do Sul, Brasil
                              and
     Computing Laboratory, University of Kent at Canterbury, U.K.

                          Tim Hopkins
     Computing Laboratory, University of Kent at Canterbury, U.K.

  Description of parameter arrays
   IPAR (INPUT)  : int
     ipar( 1): lda    (Leading dimension of a)
           2 : n      (Number of rows/columns of a)
           3 : blksz  (Size of block of data; used when data is partitioned using cyclic mode)
           4 : loclen (Number of elements stored locally;
                       *PARALLEL: Equal to at least m/nprocs or n/procs depending if row or column partitioning is used or,
                                  in the case of cyclic partitioning, it is a multiple of either m/(blksz*nprocs) or n/(blksz*nprocs).
                       *SEQUENTIAL: equal to n)
           5 : basisdim (Dimension of orthogonal basis, used in GMRES)
           6 : nprocs (Number of processors)
           7 : procid (Processor identification)
           8 : precontype (Type of preconditioning; one of
                           0 : No preconditioning,
                           1 : Left preconditioning,
                           2 : Right preconditioning,
                           3 : Symmetric preconditioning)
           9 : stoptype (Type of stopping criteria used)
          10 : maxit  (Maximum number of iterations allowed)

   IPAR (OUTPUT) : int
     ipar(11): itno   (Number of iterations executed)
          12 : status (On exit of iterative method, one of
                        0: converged to solution
                       -1: no convergence has been achieved
                       -2: "soft"-breakdown, solution may have been found
                       -3: "hard"-breakdown, no solution)
                       -4: conflict in preconditioner and stopping criterion selected
                       -5: error in stopping criterion 3, r^{T}z<0)
          13 : steperr (If status is either -2 or -3, it gives the step number of the respective algorithm
                         where a breakdown has occurred. Refer to the User's Guide for further information)

   RPAR/DPAR (INPUT)  : real/real
     spar( 1): epsilon (The value of epsilon for use in the stopping criterion)

   RPAR/DPAR (OUTPUT) : real/real
     spar( 2): exitnorm (The value of a norm of either the residual vector or the difference between two successive solution estimates according to the value of stoptype)
           3, 4 : smallest and largest eigenvalues of Q1AQ2 (in the symmetric case) OR smallest and largest real parts (in the nonsymmetric case)
           5, 6 : smallest and largest imaginary parts (only in the nonsymmetric case)
** */
