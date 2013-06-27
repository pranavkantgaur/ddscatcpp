#include "StdAfx.h"

#include "Definitions.h"
#include "DDscatParameters.h"
#include "DDscatCommons.h"
#include "DDscatMain.h"
#include "AbstractTarget.h"
#include "AbstractSolver.h"
#include "SolverPim.h"
#include "SolverTangcg.h"
#include "SolverQmrcg.h"
#include "SolverZbcg2.h"
#include "SolverSbicgm.h"
#include "Cgcommon2.h"
#include "TimerManager.h"
#include "ScatterManager.h"
#include "Cprod.h"
#include "DielectricManager.h"

template <typename T>
void Copyit(T *src, T *dst, int len)
{
	memcpy(dst, src, len*sizeof(T));
}

void Getfml(AbstractTarget *currentTarget, DDscatParameters *param, Vect3<real> &ak_tf, real ak3, Vect3<real> *aks_tf, 
	Matrix *theMatrix, Matrix *theTensor, DipoleData *theDipoleData, Vect3<Complex> &cxe01_tf, Vect3<Complex> &cxe02_tf, 
	DielectricManager *dielec, FourArray *cxfData, Complex *cxsc, Complex *cxscr1, Vect3<real> *em1_tf, Vect3<real> *em2_tf, 
	int ibeth, int ibeth1, int iphi, int iphi1, int &itask, int *itnum, PeriodicBoundaryFlag jpbc, int lace, int myid, int &navg, SumPackage &sumPackage)
{
/* **
Function of GETFML is to 
 (1) obtain a solution to the scattering problem for given target orientation, and 
 (2) return the scattering function fml in preselected scattering directions.

 Given:
       AK_TF(1-3)=(k_x,k_y,k_z)*d for incident wave [in Target Frame]
       AK3     =(kd)**3
       AKS_TF(1-3,1-NSCAT)=(k_x,k_y,k_z)*d for NSCAT scattering direction
       GAMMA = parameter controlling summation over replica dipoles (smaller alpha -> longer summation)
       BETADF(1-NAT)=orientation angle beta (radians) describing orientation of "Dielectric Frame" relative to Target Frame for dipoles 1-NAT0
       PHIDF(1-NAT)=orientation angle phi (radians) describing orientation of Dielectric Frame relative to Target Frame for dipole 1-NAT0
       THETADF(1-NAT)=orientation angle theta (radians) describing orientation of Dielectric Frame relative to Target Frame for dipoles 1-NAT0
       DX(1-3) =(d_x/d,d_y/d,d_z/d) where d_x,d_y,d_z=x,y,z lattice spacing and d=(d_x*d_y*d_z)**(1/3)=effective lattice spacing
       CALPHA  =descriptor of method used for assigning polarizabilitie = LATTDR or DRAI88 or GOBR88
       CMDSOL  =descriptor of method used for iterative solution
       CMDFFT  =descriptor of method used for FFTs
       CMDTRQ  =descriptor of whether or not to compute torques
       CSHAPE  =descriptor of target shape, needed by subroutine ALPHA
       CXE01_TF(1-3)=incident polarization state 1 at origin [in TF]
       CXE02_TF(1-3)=incident polarization state 2 at origin [in TF]
       CXEPS(1-NCOMP)=dielectric constant for compositions 1-NCOMP
       EM1_TF(1-3,1-NSCAT)=unit scat. polarization vectors 1 in TF
       EM2_TF(1-3,1-NSCAT)=unit scat. polarization vectors 2 in TF
       ETASCA  =parameter controlling number of scattering angles used for calculation of radiation force, <cos>, <cos^2>, and radiation torque
       IBETH   =MYID+IBETH1 if first time through combined BETA/THETA orientation loop
       IBETH1  =starting value of IBETH (see above)
       ICOMP(1-NAT3)=x,y,z "composition" for sites 1-NAT
       IDVOUT  =device for running output
       INIT    =0,1,2 for choice of |x0> for CCG method
       IOCC(1-NAT)=0,1 if site in extended target is vacant,occupied
       IORTH   =1,2 to do 1,2 orthogonal pol states
       IPHI    =which phi value for target orientation
       IPHI1   =starting value of IPHI (see above)
       IPBC    = 0 if PBC not used
               = 1 if PBC used
       JPBC    = 0 if PBC not used
               = 1 if PBC used in y direction only
               = 2 if PBC used in z direction only
               = 3 if PBC used in both y and z directions
               **NOTE** JPBC is not used by GETFML When JPBC > 0, GETFML computes f_ml for the Target Unit Cell.
       ITASK    not used
       IXYZ0(1-MXNAT,1-3)=lattice coordinates for sites 1-NAT (in TF) with sites 1-NAT0 in list being the occupied sites.
       DX(1-3) = location/DX(1-3) in TF corresponding to lattice site (IX,IY,IZ)=(0,0,0)
       LACE,LAXI,LCLM,LGI,LPI,LQI,LSC0=integers used to assign location in scratch array
       MXCOMP  =dimensioning info: maximum number of allowed compositions
       MXCXSC  =dimensioning info for Complex scratch space
       MXN3    =dimensioning info (3*max number of sites)
       MXNAT   =dimensioning info (max number of sites)
       MXNX    =dimensioning info (max x extent of target)
       MXNY    =dimensioning info (max y extent of target)
       MXNZ    =dimensioning info (max z extent of target)
       MXPHI   =dimensioning info (max number of phi vals for targ. ori
       MXSCA   =dimensioning info (max number of scat. dirs for f_ml)
       MYID    =parallel process identifier (=0 if only 1 process)
       NAT     =number of sites in extended target
       NAT0    =number of sites in original target
       NAT3    =3*NAT
       NCOMP   =number of different dielectric tensor elements in target
       NSCAT   =number of scat. dirs for f_ml
       NX      =x extent of extended target
       NY      =y extent of extended target
       NZ      =z extent of extended target
       PHI(1-NPHI)=phi values for target orientation
       PIA2    =\pi*(3*NAT/4\pi)^{2/3} (NAT=number of sites in original
       SHPAR(1-10)=target shape parameters, needed by subroutine ALPHA
       TOL     =tolerance for terminating conjugate gradient iteration

 and scratch space:
       CXRLOC(1-MXCOMP+1,3,3)
       CXSC(1-MXCXSC)=scratch space
       CXSCR1(1-MXN3)=scratch space
       CXZC(1->MXNX+1,1->MXNY+1,1->MXNZ+1,1-6)=scratch space
       CXZW(1-MXNAT,1-24)=scratch space
       SCRRS1(1-MXNAT,3)=scratch space
       SCRRS2(1-MXNAT)=scratch space

 Returns:
       CXADIA(1-3,1-NAT)=diagonal elements of "A matrix"
       CXAOFF(1-3,1-NAT)=off-diagonal elements of 3x3 blocks on diagonal of "A matrix"
       CXALPH(1-3,1-NAT)=diagonal elements of 3x3 polarizability tensor for dipoles 1-NAT (in TF)
       CXALOF(1-3,1-NAT)=off-diagonal elements of 3x3 polarizability tensor for dipoles 1-NAT
       CXE_TF(1-NAT3)=incident x,y,z E field at dipoles 1-NAT (in TF)
       CXF11(1-NSCAT)=scattering matrix element f_11 for NSCAT dirs
       CXF12(1-NSCAT)=                          f_12
       CXF21(1-NSCAT)=                          f_21
       CXF22(1-NSCAT)=                          f_22
       CXPOL_TF(1-NAT3)=x,y,z polarization of dipoles 1-NAT in TF
       QABS(1-2)=Q_abs=C_abs/(PIA2*d**2) for incident pols 1,2
       QBKSCA(1-2)=diff.scatt.cross section/(PIA2*d^2) for backscat for inc.pols.1,2
       QEXT(1-2)=Q_ext=C_ext/(PIA2*d**2) for incident pols 1,2
       QPHA(1-2)=phase shift cross section for incident pols 1,2
       QSCA(1-2)=Q_sca=C_sca/(PIA2*d**2) for incident pols 1,2
       QSCAG(1-3,1-2)=<cos(theta)>*Q_sca for incident pols 1,2
                      <sin(theta)*cos(phi)>*Q_sca for incident pols 1,2
                      <sin(theta)*sin(phi)>*Q_sca for incident pols 1,2
                      [QSCAG(1-3,1-2) is in Lab Frame]
       QSCAG2(1-2)   =<cos^2(theta)>*Q_sca for inciden tpols 1,2
       QTRQAB(1-3,1-2)=vector torque cross section/(PIA2*d^2) for
                      torque due to absorption, along axes 1-3, for
                      incident pols 1,2
                      [QTRQAB(1-3,1-2) is in Lab Frame]
       QTRQSC(1-3,1-2)=vector torque cross section/(PIA2*d^2) for
                      torque due to scattering, along axes 1-3, for
                      incident pols 1,2
                      [QTRQSC(1-3,1-2) is in Lab Frame]
        TIMERS(1-12)=timing information
       ITNUM(1-2)= number of iterations taken for polarizations 1 and 2
       MXITER    = maximum number of iterations allowed.
 if JO=1:
       CXSCR1(1->3*NAT0)=reduced polarization vector for incident pol. 1
       CXSC(LACE->LACE-1+3*NAT0)=reduced polarization vector for incident pol. 2
 Requires:
	MATVEC = external routine to compute A*x
	CMATVEC = external routine to compute conjg(A')*x where A is matrix such that A*pol=E and x is arbitrary vector

B.T.Draine, Princeton Univ. Observatory, 1990.

History records of Fortran vectors removed.

Copyright (C) 1993,1994,1995,1996,1997,1998,2003,2004,2006,2007,2008,
              2010,2011 B.T. Draine and P.J. Flatau
Copyright (C) 2012,2013, C++ versions, Choliy V.

This code is covered by the GNU General Public License.
** */

	const real pia2 = Pi * Pow((real)0.75 * currentTarget->Nat0() / Pi, (real)(2./3.));

	static real e02;
	int nat03 = 3 * currentTarget->Nat0();
	int nat3 = 3 * currentTarget->Nat();

	ScatterManager *scatter = new ScatterManager(currentTarget);
	int i, mxmatvec;
	real dtime, cbksca, csca, cscag2;
	Vect3<real> cscag, ctrqab, ctrqsc;
	bool nonzero_x;
	char cmsgnm[256];

//fprintf(stderr, " --- AAA ---\n");
//theDipoleData->Debug(stderr);

	Vect3<Complex> cxe0_tf, cxe01_tf1, cxe02_tf1;
	Complex cxa, cxb;
	fprintf(stderr, " >Getfml x0 = %lf %lf %lf\n", currentTarget->X0().data[0], currentTarget->X0().data[1], currentTarget->X0().data[2]);
	currentTarget->Unreduce(theDipoleData->Cxxi());

//fprintf(stderr, " --- UUU ---\n");
//theDipoleData->Debug(stderr);

// 
// Compute fml directly only if IPHI=1 Otherwise use previously computed fml values to obtain new fml
	if (iphi == 0)
	{
// store incident polarization states for use when IPHI > 1
		cxe01_tf1 = cxe01_tf;
		cxe02_tf1 = cxe02_tf;
		for(int jo=0; jo<param->Iorth(); ++jo)
		{
			cxe0_tf = (jo == 0) ? cxe01_tf : cxe02_tf;
			e02 = (real)0.;
			for(i=0; i<3; ++i)
			{
				e02 += (cxe0_tf.data[i] * cxe0_tf.data[i].conjg()).re;
			}
//
// Compute quantity used for "normalizing" error:
			AbstractSolver::Errscal() = Sqrt(currentTarget->Nat0() * e02);
//
// Call EVALE for NAT sites
			currentTarget->Evale(cxe0_tf, ak_tf, theDipoleData->Cxe_tf());
//			fprintf(stderr, " --- BBB ---\n");
//			theDipoleData->Debug(stderr);
//
// Call ALPHA to determine polarizabilities at NAT sites
//     (note: this has to be within the loop over directions and
//      polarizations because Lattice Dispersion Relation polarizabilities
//      depend on direction of propagation and on polarization state)
			Alphadiag(currentTarget, ak_tf, theTensor, cxe0_tf, dielec, myid);
//
// Call EVALA to prepare A matrix elements
			theMatrix->Evala(theTensor, currentTarget->Nat());

// PIMSETPAR sets parameters used by PIM package
// ipar(1) LDA         LDA (Leading dimension of a)
// ipar(2) N           N   (Number of rows/columns of a)
// ipar(3) BLKSZ       N   (Size of block of data; used when data is partitioned using cyclic mode)
// ipar(4) LOCLEN      N   (Number of elements stored locally; for sequential=n)
// ipar(5) BASISDIM    C=basis=10 (Dimension of orthogonal basis, used in GMRES)
// ipar(6) NPROCS      -1 (Number of processors)
// ipar(7) PROCID      -1 (Processor identification)
// ipar(8) PRECONTYPE PRET (=1 Type of preconditioning 0-no preconditioning, 1-left, 2-right, 3-symmetric
// ipar(9) STOPTYPE   STOPT (Type of stopping criteria used)
// ipar(10) MAXIT     MAXIT (=int(n/2) Maximum number of iterations allowed
//
// on return from PIM, following are defined:
//      IPAR(11) = itno   (Number of iterations executed)
//      IPAR(12) = exit status 
//                 0: converged
//                -1: no convergence has been achieved
//                -2: soft breakdown, solution may have been found
//                -3: hard breakdown, no solution
//                -4 to -11 : other conditions
//      IPAR(13) = if IPAR(12) = -2 or -3, gives the step number in the
//                 algorithm where a breakdown has occurred.

// Set upper limit on iterations to smaller of 1e4 and NAT0

//			mxiter = min_(10000, nat03);

//  pass information on max.no. of iterations allowed to DDCOMMON_9

// BTD 080716: changed
//             ITERMX=IPAR(10)
//			AbstractSolver::Itermx() = mxiter;
//
// Clean CXE_TF:
			if (currentTarget->Nat0() < currentTarget->Nat())
				Nuller(theDipoleData->Cxe_tf(), currentTarget->Iocc(), currentTarget->Nat());
//fprintf(stderr, " --- CCC ---\n");
//theDipoleData->Debug(stderr);
/* **
{
	fprintf(stderr, "Nat = %d nat0 = %d\n", currentTarget->Nat(), currentTarget->Nat0());
	for(int iu=0; iu<3*currentTarget->Nat(); ++iu)
	{
		fprintf(stderr, "%3d %20.8e %20.8e %20.8e %20.8e %20.8e %20.8e %20.8e %20.8e %20.8e %20.8e %20.8e %20.8e %3d\n", 
			iu, theDipoleData->Cxxi(iu).re, theDipoleData->Cxxi(iu).im, 
			theDipoleData->Cxe_tf(iu).re, theDipoleData->Cxe_tf(iu).im, 
			theMatrix->Diagonal(iu).re, theMatrix->Diagonal(iu).im, 
			theMatrix->OffDiagonal(iu).re, theMatrix->OffDiagonal(iu).im,
			theTensor->Diagonal(iu).re, theTensor->Diagonal(iu).im, 
			theTensor->OffDiagonal(iu).re, theTensor->OffDiagonal(iu).im,
			currentTarget->Icomp().Value(iu/3, iu%3));
	}
}
** */
//
// Iterate to improve CXPOL_TF
			real tolr = (real)1.;
			switch(param->Cmdsol())
			{
			case SolMethod_PETRKP:
				{
					SolverPim *solver = (SolverPim *)(AbstractSolver::GetInstance(SolMethod_PETRKP));
// Petr					solver->SetParameters(nat03, nat03, nat03, 10, -1, -1, 1, 5, param->Mxiter(), param->Tol());
					solver->GetCgStruct().QuickParam(nat03, Stoptype_RltEpsilonB, Precontype_None, param->Tol(), true, param->Mxiter(), 2);
					Cinit(nat03, Complex(), theDipoleData->Cxxi(), 1);
					dtime = Timeit("Petrkp");
					solver->SetMatvecFunctions(Matvec, Cmatvec, NULL);
					solver->SetNormFunctions(NULL, Pscnrm2);
					solver->SetProgressFunction(Progress);
// Petr 					solver->Petr(theDipoleData->Cxxi(), theDipoleData->Cxe_tf());					// we need to call ->SetParameters before Petr
					solver->Petr90(theDipoleData->Cxxi(), theDipoleData->Cxe_tf(), nat3, cxsc);					// we need to call ->QuickParam before Petr90
					dtime = Timeit("Petrkp");
					TimerManager::GetInstance()->SetValue(0, dtime);
					TimerManager::GetInstance()->SetValue(1, (real)solver->GetNumberOfIterations());
					itnum[jo] = solver->GetNumberOfIterations();
				}
				break;

			case SolMethod_SBICGM:
				{
					SolverSbicgm *solver = (SolverSbicgm *)(AbstractSolver::GetInstance(SolMethod_SBICGM));
					solver->GetCgStruct().QuickParam(nat03, Stoptype_RltEpsilonB, Precontype_None, param->Tol(), true, param->Mxiter(), 2);
					Cinit(nat03, Complex(), theDipoleData->Cxxi(), 1);
					dtime = Timeit("Sbicgm");
					solver->SetMatvecFunctions(Matvec, Cmatvec, NULL);
					solver->Sbicg90(theDipoleData->Cxxi(), theDipoleData->Cxe_tf(), nat03, cxsc);
					dtime = Timeit("Sbicgm");
					TimerManager::GetInstance()->SetValue(0, dtime);
					TimerManager::GetInstance()->SetValue(1, (real)solver->GetNumberOfIterations());
				}
				break;

			case SolMethod_PBCGS2:			// Tested Ok
				{
					dtime = Timeit("Pbcgs2");
// !BTD 07.08.15 Use TOLR to input TOL with ZBCG2 (ZBCG2 returns actual achieved tolerance in TOLR)
					tolr = param->Tol();
// correspondence with DDSCAT
// xi - initial guess of cxpol on input; cxpol  on output
// b  - cxe (right hand side, not changed)
// xr - A xi (use matvec)  work vector
// lda -  nat3
// ndim - nat3 ?
// nlar  work array dimension >= 12
// wrk - cxsc NOTE!!! THAT MXCXSC has to be set =12 in main DDSCAT(large!)
// maxit - itermx ? mxiter?
// nloop - itern
// tol -   tolr (input)
// tole -  achieved relative error
// ipar(12) - convergence;  ipar(12)=0 converged
					mxmatvec = 4 * param->Mxiter();
					nonzero_x = false;
					SolverZbcg2 *solver = (SolverZbcg2 *)AbstractSolver::GetInstance(SolMethod_PBCGS2);
					solver->SetMatvecFunction(Matvec);
					solver->SetParameters(2, nat3, param->Mxiter());
					Cinit(nat3, Complex(), theDipoleData->Cxxi(), 1);
					solver->Zbcg2(theDipoleData->Cxxi(), theDipoleData->Cxe_tf(), true, nonzero_x, tolr, mxmatvec);
					itnum[jo] = solver->GetNumberOfIterations();
					dtime = Timeit("Pbcgs2");
				}
				break;

			case SolMethod_GPBICG:
				{
					Cinit(nat03, Complex(), theDipoleData->Cxxi(), 1);
					dtime = Timeit("Pbcgs2");
					tolr = param->Tol();
//
// correspondence with DDSCAT
// xi - initial guess of cxpol on input; cxpol  on output
// b  - cxe (right hand side, not changed)
// xr - A xi (use matvec)  work vector
// lda -  nat3
// ndim - nat3 ?
// nlar  work array dimension >= 12
// wrk - cxsc NOTE!!! THAT MXCXSC has to be set =12 in main DDSCAT(large!)
// maxit - itermx ? mxiter?
// nloop - itern
// tol -   tolr (input)
// tole -  achieved relative error
// ipar(12) - convergence;  ipar(12)=0 converged
// TODO:
// change precision
// change mat to matrix multiplication
					SolverTangcg *solver = (SolverTangcg *)(AbstractSolver::GetInstance(SolMethod_GPBICG));
					solver->SetMatvecFunction(Matvec);
					solver->SetParameters(nat03, param->Mxiter());
					solver->Tangcg(theDipoleData->Cxxi(), theDipoleData->Cxe_tf(), param->Tol(), tolr);
					itnum[jo] = solver->GetNumberOfIterations();
					dtime = Timeit("Gpbicg");
				}
				break;

			case SolMethod_QMRCCG:
				{
					Cinit(nat03, Complex(), theDipoleData->Cxxi(), 1);
	   				dtime = Timeit("Qmrccg");
					SolverQmrcg *solver = (SolverQmrcg *)(AbstractSolver::GetInstance(SolMethod_QMRCCG));
					solver->SetParameters(nat03, param->Mxiter());
					solver->SetMatvecFunction(Matvec);
					solver->Pimqmrcg(theDipoleData->Cxxi(), theDipoleData->Cxe_tf(), param->Tol(), tolr);
					itnum[jo] = solver->GetNumberOfIterations();
					dtime = Timeit("Qmrccg");
				}
				break;

			case SolMethod_PBCGST:
				{
// CALL PIMSGETPAR(IPAR,SPAR,LDA,N,BLKSZ,LOCLEN,BASISDIM,NPROCS,PROCID,PRECONTYPE,STOPTYPE,MAXIT,ITNO,STATUS,STEPERR,EPSILON,EXITNORM)
// flatau warning something is wring with timers because they also time eself - check this we need to write better timming routine
// flatau I am just hardwiring number of iterations (mxiter is not used anymore)
// BTD 07.08.12 minor revision
					SolverPim *solver = (SolverPim *)(AbstractSolver::GetInstance(SolMethod_PBCGST));
					int no_cg_restart = 5;
					int no_cg_iter = 0;
					solver->SetParameters(nat03, nat03, nat03, 10, -1, -1, 1, 2, no_cg_restart, param->Tol());
					Cinit(nat03, Complex(), theDipoleData->Cxxi(), 1);
					dtime = Timeit("Pbcgst");
					solver->SetMatvecFunctions(Matvec, NULL, NULL);
					solver->SetPreFunctions(Diagl, NULL);
					solver->SetNormFunctions(Pcsum, Pscnrm2);
					solver->SetProgressFunction(Progress);
					for(int iter=0; iter<=(param->Mxiter() / no_cg_restart); ++iter)
					{
						if (solver->GetIntParameter(PimIparNameStatus) != 0)
						{
							if(solver->GetNumberOfIterations() >= param->Mxiter())
							{
                                sprintf(cmsgnm, "itern=%6d > mxiter=%6d\n", solver->GetNumberOfIterations(), param->Mxiter());
								Errmsg("Fatal", "Getfml", cmsgnm);
							}
							if (iter > 0)
								fprintf(stderr, "restart PimCbicgstab: %d\n", iter);
							solver->PimCbicgstab(theDipoleData->Cxxi(), theDipoleData->Cxe_tf());
						}
					}
					itnum[jo] = solver->GetNumberOfIterations();
					dtime = Timeit("Pbcgst");
					TimerManager::GetInstance()->SetValue(0, dtime);
					TimerManager::GetInstance()->SetValue(1, (real)no_cg_iter);
				}
				break;

			default:
				fprintf(stderr, "Error -- invalid CMDSOL in getfml");
				return;
				break;
			}
//
// Reduce the vectors CXADIA, CXAOFF, CXPOL_TF and CXE_TF to retain only occupied sites, to eliminate unnecessary calculations for
// unoccupied sites when calling EVALQ and SCAT

//fprintf(stderr, " --- QQQ ---\n");
//theDipoleData->Debug(stderr);

			currentTarget->Reduce(theMatrix->Diagonal());
			currentTarget->Reduce(theMatrix->OffDiagonal());
			currentTarget->Reduce(theDipoleData->Cxe_tf());
			currentTarget->Reduce(theDipoleData->Cxxi());

//fprintf(stderr, " --- XXX ---\n");
//theDipoleData->Debug(stderr);

//
// Store reduced polarization vector P for use in subsequent computations for other phi values
//     CXSCR1     for JO=1
//     CXSC(LACE) for JO=2 (this scratch space no longer needed by solver)
			if (jo == 0)
			{
				Copyit(theDipoleData->Cxxi(), cxscr1, nat03);
			}
			if (jo == 1)
			{
				Copyit(theDipoleData->Cxxi(), cxsc + lace, nat03);
			}
			theMatrix->Evalq(theDipoleData, ak_tf, nat03, e02, sumPackage.Qabs().Q()[jo], sumPackage.Qext().Q()[jo], sumPackage.Qpha().Q()[jo], 1);
			sumPackage.Qabs().Q()[jo] /= pia2;
			sumPackage.Qext().Q()[jo] /= pia2;
			sumPackage.Qpha().Q()[jo] /= pia2;
	  		sprintf(cmsgnm, " Q_abs = %11.4e Q_ext= %11.4e Q_pha= %11.4e", sumPackage.Qabs().Q()[jo], sumPackage.Qext().Q()[jo], sumPackage.Qpha().Q()[jo]);
            Wrimsg("Getfml", cmsgnm);
//
// First call to timeit:
			dtime = Timeit(" Scat ");
			if (jo == 0)
			{
// Call SCAT to compute CXF11,CXF21
//fprintf(stderr, "000 111 000\n");
				scatter->Scat(ak_tf, aks_tf, em1_tf, em2_tf, e02, param->Etasca(), cbksca, csca, cscag, cscag2, ctrqab, ctrqsc, 
					theDipoleData, cxe01_tf, cxfData->Cx11(), cxfData->Cx21(), myid, jpbc, navg, param->Nscat());
			}
			if (jo == 1)
			{
// Call SCAT to compute CXF12,CXF22
//fprintf(stderr, "000 222 000\n");
				scatter->Scat(ak_tf, aks_tf, em1_tf, em2_tf, e02, param->Etasca(), cbksca, csca, cscag, cscag2, ctrqab, ctrqsc, 
					theDipoleData, cxe01_tf, cxfData->Cx12(), cxfData->Cx22(), myid, jpbc, navg, param->Nscat());
			}
// Second call to timeit
			dtime = Timeit(" Scat ");
			TimerManager::GetInstance()->SetValue(3, dtime);
//
			if (jpbc == PeriodicNo)
			{
// Note: at present, SCAT only calculates CBKSCA and CSCA when JPBC=0 therefore only calculate QBKSCA(JO), etc, when JPBC=0
				sumPackage.Qbksca().Q()[jo] = cbksca / pia2;
				sumPackage.Qsca().Q()[jo] = csca / pia2;
//
// Since Q_ext and Q_abs are calculated "exactly", when particles are "large" we obtain Q_sca from (Q_ext-Q_abs) rather than
// scattering cross section computed by SCAT.  This is particularly advantageous when the particles are large with Complex scattering
// patterns, so that angular integration of the differential scatterin cross section would be either very time-consuming or inaccurate.
// However, when particles are "small" the scattering cross section is small compared to Q_ext and Q_abs, so that the difference
// (Q_ext-Q_abs) may be inaccurate.  Accordingly, when particles are "small" we instead accept the scattering cross section
// computed by SCAT (probably quite accurate when albedo is small beca the scattering pattern tends to be dipolar and therefore well-resol
// by a reasonable number of scattering directions ICTHM and IPHI). 
// We use an ad-hoc procedure to smoothly make the transitions between these regimes, based on the albedo.
// The "dividing" line of albedo=.03 is chosen on basis of numerical experiments.  If Q_ext and Q_abs are each calculated to
// accuracy of 1 times 10^{-4} (say), then the difference between Q_ext and Q_abs would be accurate to
// 0.0001*Q_ext*sqrt(2) .  When albedo=.03 (Q_ext-Q_abs)/Q_sca would then be accurate to .0001*sqrt(2)/.03=.005
				real falb = sumPackage.Qsca().Q()[jo] / ((real)0.03 * sumPackage.Qext().Q()[jo]);
				falb = falb * falb;
				sumPackage.Qsca().Q()[jo] = (sumPackage.Qsca().Q()[jo] + (sumPackage.Qext().Q()[jo] - sumPackage.Qabs().Q()[jo]) * falb) / ((real)1. + falb);
//
// Compute QSCAG(1-3,JO), QTRQSC(1-3,JO), QSCAG2
				sumPackage.Qscag().Q()[jo] = cscag * sumPackage.Qsca().Q()[jo] / csca;
				sumPackage.Qscag2().Q()[jo] = cscag2 * sumPackage.Qsca().Q()[jo] / csca;
				if (param->Cmdtrq() == TorqMethod_DOTORQ)
				{
					sumPackage.Qtrqab().Q()[jo] = ctrqab / pia2;
					sumPackage.Qtrqsc().Q()[jo] = ctrqsc / pia2;
				}
			}												// end if(jpbc==0)
		}													// end do j=1,iorth for iphi=1
	}
	else												// iphi >= 2	
	{
// Have previously computed dipole polarization vectors for first phi value (and reduced them to NAT03 elements).
// Use these two solutions to obtain Qabs,Qext, Qpha, Qsca, g*Qsca, and fml for any additional values of target rotation phi.
//     CXSCR1     contains polarization vector for IPHI=1 and JO=1
//     CXSC(LACE) contains polarization vector for IPHI=1 and JO=2
//
// Incident polarization state 1 = CXE01_TF 
// Call EVALE to obtain E at NAT0 occupied lattice sites first call to timing routine
		dtime = Timeit(" Evale");
//
// Call EVALE for NAT0 occupied sites
		currentTarget->Evale(cxe01_tf, ak_tf, theDipoleData->Cxe_tf());
//
// Second call to timing routine
		dtime = Timeit(" Evale");
		TimerManager::GetInstance()->SetValue(4, dtime);
//
// Compute polarizabilities at NAT sites, then reduce to NAT0 occupied sites:
		dtime = Timeit(" Alpha");
		Alphadiag(currentTarget, ak_tf, theTensor, cxe01_tf, dielec, myid);
		currentTarget->Reduce(theTensor->Diagonal());
		currentTarget->Reduce(theTensor->OffDiagonal());
		dtime = Timeit(" Alpha");
		TimerManager::GetInstance()->SetValue(5, dtime);
//
// Construct polarization for incident pol. state 1
		cxa.clear();
		cxb.clear();
		for(i=0; i<3; ++i)
		{
			cxa += cxe01_tf1.data[i].conjg() * cxe01_tf.data[i];
			cxb += cxe02_tf1.data[i].conjg() * cxe01_tf.data[i];
		}
		for(i=0; i<nat03; ++i)
		{
			theDipoleData->Cxxi()[i] = cxscr1[i] * cxa + (cxsc + lace)[i] * cxb;
		}
//
		dtime = Timeit(" Evalq");											// First call to timing routine
		theMatrix->Evalq(theDipoleData, ak_tf, nat03, e02, sumPackage.Qabs().Q()[0], sumPackage.Qext().Q()[0], sumPackage.Qpha().Q()[0], 1);
		dtime = Timeit(" Evalq");											// Second call to timing routine
		TimerManager::GetInstance()->SetValue(6, dtime);
 		sumPackage.Qabs().Q()[0] /= pia2;
		sumPackage.Qext().Q()[0] /= pia2;
		sumPackage.Qpha().Q()[0] /= pia2;
		dtime = Timeit(" Scat");											// First call to timing routine
//
// Now call SCAT to compute CXF11,CXF21 CXPOL_TF has been reduced to NAT03 elements first NAT03 elements of IXYZ0 correspond to physical sites
//fprintf(stderr, "000 333 000\n");
		scatter->Scat(ak_tf, aks_tf, em1_tf, em2_tf, e02, param->Etasca(), cbksca, csca, cscag, cscag2, ctrqab, ctrqsc, 
			theDipoleData, cxe01_tf, cxfData->Cx11(), cxfData->Cx21(), myid, jpbc, navg, param->Nscat());
		dtime = Timeit(" Scat");											// Second call to timing routine
		TimerManager::GetInstance()->SetValue(7, dtime);

		sumPackage.Qbksca().Q()[0] = cbksca / pia2;
		sumPackage.Qsca().Q()[0] = csca / pia2;
//
// As above, compute Q_sca from either SCAT or (Q_ext-Q_abs) depending on whether albedo << .03 or albedo >> .03
		real falb = sumPackage.Qsca().Q()[0] / ((real)0.03 * sumPackage.Qext().Q()[0]);
		falb = falb * falb;
		sumPackage.Qsca().Q()[0] = (sumPackage.Qsca().Q()[0] + (sumPackage.Qext().Q()[0] - sumPackage.Qabs().Q()[0]) * falb) / ((real)1. + falb);
		sumPackage.Qscag().Q()[0] = cscag * sumPackage.Qsca().Q()[0] / csca;
		sumPackage.Qtrqab().Q()[0] = ctrqab / pia2;
		sumPackage.Qtrqsc().Q()[0] = ctrqsc / pia2;
		sumPackage.Qscag2().Q()[0] = sumPackage.Qsca().Q()[0] * cscag2 / csca;
//
// Polarization state 2:
		dtime = Timeit(" Evale");											// First call to timing routine
//
// Call EVALE for NAT0 occupied sites to obtain appropriate E vector for pol. state 2
		currentTarget->Evale(cxe02_tf, ak_tf, theDipoleData->Cxe_tf());
		dtime = Timeit(" Evale");											// Second call to timing routine
		TimerManager::GetInstance()->SetValue(8, dtime);
//
// Compute polarizabilities at NAT sites, then reduce to NAT0 occupied sites:
		dtime = Timeit(" Alpha");
		Alphadiag(currentTarget, ak_tf, theTensor, cxe02_tf, dielec, myid);
		currentTarget->Reduce(theTensor->Diagonal());
		currentTarget->Reduce(theTensor->OffDiagonal());
		dtime = Timeit(" Alpha");
		TimerManager::GetInstance()->SetValue(9, dtime);
//
// Construct incident polarization vector for pol. state 2
		cxa.clear();
		cxb.clear();
		for(i=0; i<3; ++i)
		{
			cxa += cxe01_tf1.data[i].conjg() * cxe02_tf.data[i];
			cxb += cxe02_tf1.data[i].conjg() * cxe02_tf.data[i];
		}
		for(i=0; i<nat03; ++i)
		{
			theDipoleData->Cxxi()[i] = cxscr1[i] * cxa + (cxsc + lace)[i] * cxb;
		}
//
		dtime = Timeit(" Evalq ");											// First call to timing routine
		theMatrix->Evalq(theDipoleData, ak_tf, nat03, e02, sumPackage.Qabs().Q()[1], sumPackage.Qext().Q()[1], sumPackage.Qpha().Q()[1], 1);
		dtime = Timeit(" Evalq ");											// Second call to timing routine
		TimerManager::GetInstance()->SetValue(10, dtime);
		sumPackage.Qabs().Q()[1] /= pia2;
		sumPackage.Qext().Q()[1] /= pia2;
		sumPackage.Qpha().Q()[1] /= pia2;
//
// Call SCAT to compute CXF12,CXF22 CXPOL_TF has NAT03 elements IXYZ0 was initially given NAT03 elements
		dtime = Timeit(" Scat");											// First call to timing routine
//fprintf(stderr, "000 444 000\n");
		scatter->Scat(ak_tf, aks_tf, em1_tf, em2_tf, e02, param->Etasca(), cbksca, csca, cscag, cscag2, ctrqab, ctrqsc, 
			theDipoleData, cxe01_tf, cxfData->Cx12(), cxfData->Cx22(), myid, jpbc, navg, param->Nscat());
		dtime = Timeit(" Scat");											// Second call to timing routine
		TimerManager::GetInstance()->SetValue(11, dtime);
//
		sumPackage.Qbksca().Q()[1] = cbksca / pia2;
		sumPackage.Qsca().Q()[1] = csca / pia2;
//
// As above, compute Q_sca from either SCAT or (Q_ext-Q_abs) depending on whether albedo << .03 or albedo >> .03
		falb = sumPackage.Qsca().Q()[1] / ((real)0.03 * sumPackage.Qext().Q()[1]);
		falb = falb * falb;
		sumPackage.Qsca().Q()[1] = (sumPackage.Qsca().Q()[1] + (sumPackage.Qext().Q()[1] - sumPackage.Qabs().Q()[1]) * falb) / ((real)1. + falb);
		sumPackage.Qscag().Q()[1] = cscag * sumPackage.Qsca().Q()[1] / csca;
		sumPackage.Qtrqab().Q()[1] = ctrqab / pia2;
		sumPackage.Qtrqsc().Q()[1] = ctrqsc / pia2;
		sumPackage.Qscag2().Q()[1] = sumPackage.Qsca().Q()[1] * cscag2 / csca;
	}
	sprintf(cmsgnm, " returning from getfml");
	Wrimsg("Getfml", cmsgnm);

	delete scatter;
	AbstractSolver::DeleteSolver();
}
