#ifndef __CGSTRUCT_H__
#define __CGSTRUCT_H__

#include "Solverlib.h"
#include "Definitions.h"

//
// none = 0			// No preconditioning'
// left = 1			// Left preconditioning'
// right = 2		// Right preconditioning'
// symmetric = 3	// Symmetric preconditioning
typedef enum _Precontype
{
	Precontype_None, Precontype_Left, Precontype_Right, Precontype_Symmetric, Precontype_End	
} Precontype;

//
// r_eps = 1			// Stopping criterion: ||r||<epsilon'
// r_eps_b = 2			// Stopping criterion: ||r||<epsilon||b||'
// sqrt_rz_eps_b = 3	// Stopping criterion : sqrt(r''z)<epsilon ||b||'
// z_eps = 4			// Stopping criterion : ||z||<epsilon'
// z_eps_b = 5			// Stopping criterion : ||z||<epsilon ||b||'
// z_eps_qb = 6			// Stopping criterion : ||z||<epsilon ||qb||'
// diff_iter_eps = 7	// Stopping criterion : ||x(k)-x(k-1)||<epsilon'
typedef enum _Stoptype
{
	Stoptype_None, 
	Stoptype_RltEpsilon, Stoptype_RltEpsilonB, Stoptype_SqrtRltEpsilonB, 
	Stoptype_ZltEpsilon, Stoptype_ZltEpsilonB, Stoptype_ZltEpsilonBQ, 
	Stoptype_TwoXltEpsilon, Stoptype_End
} Stoptype;

//
// converged           =  0		// Exit status: converged to solution'
// no_convergence      = -1		// Exit status: no convergence achieved'
// soft_breakdown      = -2		// Exit status: soft breakdown solution may have been found'
// hard_breakdown      = -3		// Exit status: hard breakdown no solution found'
// precon_stop_error   = -4		// Conflict between preconditioner and selecte stopping criterion'
// stop_error          = -5		// Error in stopping criterion, r''z<0'
typedef enum _StopStatus
{
	StopStatus_Converged, StopStatus_NoConvergence, StopStatus_SoftBreakdown, 
	StopStatus_HardBreakdown, StopStatus_PreconStopError, StopStatus_StopError, 
	StopStatus_End
} StopStatus;

class SOLVERLIB_API CgStruct
{
public:
	int n, blksz, loclen, basisdim, nprocs, procid, maxit, itno, steperr;
	real epsilon_err, exitnorm, mineig_real, maxeig_real, mineig_imag, maxeig_imag;
	Precontype precontype;
	Stoptype stoptype;
	StopStatus status;
	FILE *ioerr;
	bool print;

public:
	CgStruct(void);
	~CgStruct(void);

public:
	void QuickParam(int n, Stoptype stoptype, Precontype precontype, real tol, bool bPrint, int maxiter, int basisdim);
};

#endif // __CGSTRUCT_H__
