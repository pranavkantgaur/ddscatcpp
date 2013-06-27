#include "StdAfx.h"

#include "Definitions.h"
#include "DDscatMain.h"
#include "DDscatCommons.h"
#include "DDscatParameters.h"
#include "TargetManager.h"
#include "Cprod.h"
#include "GreenFunctionManager.h"

void Diagl(Complex *u, Complex *v, int *)
{
/* **
Cg package

Left preconditioning subroutine - division by diagonal
History:
Fortran versions history removed.
end history
Copyright (C) 1993,2000,2004,2005,2006,2007,2008, B.T. Draine and P.J. Flatau
This code is covered by the GNU General Public License
** */
	Matrix *theMatrix = DDscatParameters::GetInstance()->GetMatrix();
	Cdiv(u, v, theMatrix->Diagonal(), theMatrix->GetSize());
}

void Cdiv(Complex *u, Complex *v, Complex *cxa, int n)
{
/* **
Part of conjugate gradient package needed by left preconditioning subroutine
Copyright (C) 2003 B.T. Draine and P.J. Flatau
This code is covered by the GNU General Public License
** */

	for(int i=0; i<n; ++i)
	{
		v[i] = u[i] / cxa[i];
	}
}

void Matvec(Complex *cxx, Complex *cxy, int *)
{
/* **
matvec  --> cxy = A cxx          'N'
tmatvec --> cxy = A' cxx         'T'
cmatvec --> cxy = conjg(A') cxx  'C'

Notice that DDSCAT has cases 'N' and 'C' implemented in CPROD
Notice, however, that in DDSCAT 'T'='N'.
Some CG methods use cases N, C (like Petravic).
Some CG methods use cases N, T (like PIM)
Some CG use all three products N, C, T  (parts of CCGPAK)

Information about matrix A is transfered via commons /M1/-/M7/ with
the main program. Arrays cxadia, cxzc, cxzwm, iocc are properly
dimensioned inside "cprod". You may get warning during the compilation
about "misaligned common". Ignore these. Also notice mxnatf,
mxn3f,mxnxf, mxnyf, mxnzf. These are defined by "parameter
statement" in main code and this is the reason why they have "f" suffix
ipar - information array, not used, defined for compatibility

History:
Fortran versions history removed.
Copyright (C) 1993,1997,2000,2004,2005,2006,2007,2008,2011, B.T. Draine and P.J. Flatau
This code is covered by the GNU General Public License.
** */
	AbstractTarget *currentTarget = TargetManager::GetInstance()->GetCurrentTarget();
	DDscatParameters *param = DDscatParameters::GetInstance();
	Matrix *theMatrix = param->GetMatrix();
	Cprod(Common1::GetInstance()->Ak_tf(), param->Gamma(), currentTarget->Dx(), param->Cmdfft(), 'n', theMatrix->Diagonal(), theMatrix->OffDiagonal(), 
		cxx, cxy, Common4::GetInstance()->Cxzw(), currentTarget->Iocc(), currentTarget->Nat(), currentTarget->Nat0(), 3 * currentTarget->Nat(), 
		currentTarget->Pyd(), currentTarget->Pzd());
}

void Cmatvec(Complex *cxx, Complex *cxy, int *)
{
/* **
Subroutine CMATVEC

matvec  --> cxy = A cxx
tmatvec --> cxy = A' cxx
cmatvec --> cxy = conjg(A') cxx

History:
Fortran versions history removed.
Copyright (C) 1993,1997,2000,2004,2005,2006,2007,2008,2011, B.T. Draine and P.J. Flatau
This code is covered by the GNU General Public License.
** */
	AbstractTarget *currentTarget = TargetManager::GetInstance()->GetCurrentTarget();
	DDscatParameters *param = DDscatParameters::GetInstance();
	Matrix *theMatrix = param->GetMatrix();
	Cprod(Common1::GetInstance()->Ak_tf(), param->Gamma(), currentTarget->Dx(), param->Cmdfft(), 'c', theMatrix->Diagonal(), theMatrix->OffDiagonal(), 
		cxx, cxy, Common4::GetInstance()->Cxzw(), currentTarget->Iocc(), currentTarget->Nat(), currentTarget->Nat0(), 3 * currentTarget->Nat(), 
		currentTarget->Pyd(), currentTarget->Pzd());
}

void Cprod(Vect3<real> &akr, real gamma, Vect3<real> &dx, FftMethod cmethd, char cwhat, Complex *cxadia, Complex *cxaoff, 
	Complex *cxx, Complex *cxy, Array4Stacked<Complex> *cxzw, bool *iocc, int nat, int nat0, int nat3, real pyd, real pzd)
{
/* **
Given:
   AKR(1-3)= k(1-3)*d , where k=k vector in vacuo d=effective lattice spacing=(dx*dy*dz)**(1/3)
   GAMMA  = parameter controlling PBC summations over replica dipoles interaction is suppressed by factor exp(-(gamma*k*r)^4) and summations are carried out to r=2/(gamma*k)
   DX(1-3)=(dx/d, dy/d, dz/d) where dx,dy,dz=lattice spacing in x,y,z directions. 
		By definition of effective lattice spacing d, must have DX(1)*DX(2)*DX(3)=1
   CMETHD         = character variable specifying method used for FFT evaluation (see ESELF)
   CWHAT          = 'N' or 'C' determines what is to be done (see below
   CXADIA(1-NAT3) = diagonal elements of A(j,k) (diagonal elements of inverse of polarizability tensor for dipoles 1-NAT)
   CXAOFF(1-NAT3) = off-diagonal elements of 3x3 diagonal blocks of A(j,k): CXAOFF(J,1-3)=a_{23},a_{31},a_{12} from 3x3 matrix a_{ij} for dipole J.
		Recall that a_{ij} is inverse of 3x3 polarizability tensor alpha_{ij} for dipole J.
   CXX(1-NAT3)    = Complex vector
   CXZC           = array of Fourier transformed Green function coefficients used internally by ESELF, and generate by ESELF each time called with new k vector
		If IPBC=0: size = (NX+1)*(NY+1)*(NZ+1)*6 If IPBC=1: size = (2*NX)*(2*NY)*(2*NZ)*6
   *****Not to be overwritten between calls to ESELF**
   CXZW           = Complex workspace required by ESELF
   IDVOUT      = device number for output
   IOCC(1-NAT) = 0 or 1 if site is unoccupied or occupied
   MXN3,MXNX,MXNY,MXNZ = dimensioning information (see DDSCAT)
   NAT         = number of sites in extended target
   NAT0        = number of occupied sites in extended target
   NAT3        = 3*NAT
   NX,NY,NZ    = extended target size is (NX*DX(1)) by (NY*DX(2)) by (NZ*DX(3))
   PYD         = 0. to do isolated target = (periodicity in y direction)/d(2) for periodic b.c.
   PZD         = 0. to do isolated target = (periodicity in z direction)/d(3) for periodic b.c.

 Returns:
   CXY(1-NAT3)
       If CWHAT = 'N' : CXY(j)=A(j,k)*CXX(k) (sum over k)
                = 'C' : CXY(j)=conjg(A(k,j))*CXX(k) (sum over k)

It is assumed that matrix A(j,k) is symmetric Diagonal elements of A(j,k) are in vector CXADIA(j)
Off-diagonal terms coupling dipole component x with direction y, etc are stored in vector CXAOFF(j), CXAOFF(j+NAT), CXAOFF(j+2*NAT)
Other off-diagonal elements of A(j,k) are produced internally by subroutine ESELF, which uses 3d-FFT to accomplish convolution
of CXX with A-matrix terms coupling dipole j with other dipoles (i.e., E field at j due to other dipoles)

Original version of CPROD created by P.J.Flatau, Colorado State Univ.
Subsequently modified by B.T.Draine, Princeton Univ. Obs.
History:
Fortran versions history removed.
Copyright (C) 1993,1997,1998,2000,2004,2005,2006,2007,2008,2011 B.T. Draine and P.J. Flatau
This code is covered by the GNU General Public License.
** */
//
	real akd = akr.Mod();
//
// We assume that input vector CXX has zeroes for elements corresponding to vacuum sites.
	int j1;
	switch(cwhat)
	{
	case 'n':
	case 'N':						// CXX(1-NAT3) = polarizations
	{
//int iii;
//fprintf(stderr, " - BEFORE - %d\n", nat3);
//for(iii=0; iii<nat3; ++iii)
//{
//	fprintf(stderr, "%d %lf %lf %lf %lf %lf %lf %lf %lf\n", iii, cxx[iii].re, cxx[iii].im, cxy[iii].re, cxy[iii].im, cxadia[iii].re, cxadia[iii].im, cxaoff[iii].re, cxaoff[iii].im);
//}
		GreenFunctionManager::GetInstance()->GetElectric()->Self(cmethd, cxx, gamma, pyd, pzd, akr, akd, dx, cxzw, cxy);

//fprintf(stderr, " - AFTER - %d\n", nat3);
//for(iii=0; iii<nat3; ++iii)
//{
//	fprintf(stderr, "%d %lf %lf %lf %lf %lf %lf %lf %lf\n", iii, cxx[iii].re, cxx[iii].im, cxy[iii].re, cxy[iii].im, cxadia[iii].re, cxadia[iii].im, cxaoff[iii].re, cxaoff[iii].im);
//}

// CXY(1-NAT3) is now Efield produced by other dipoles (including replica dipoles if PYZ or PZD are nonzero)
// ESELF omitted contribution from diagonal terms: add these. Also remember that ESELF computes -A_jk*x_k
		for(j1=0; j1<nat3; j1++)
		{
			cxy[j1] = cxadia[j1] * cxx[j1] - cxy[j1];
		}
// Add contribution from off-diagonal elements of 3x3 diagonal blocks
// Version below assumes that off-diagonal elements a_ij are stored as
// CXAOFF(J,1-3)=(a_23,a_31,a_12) for dipole J.
// 97.12.28 (BTD):
		for(j1=0; j1<nat; ++j1)
		{
			int index = 3*j1;
			cxy[index  ] += (cxaoff[index+1] * cxx[index+2] + cxaoff[index+2] * cxx[index+1]);
			cxy[index+1] += (cxaoff[index+2] * cxx[index  ] + cxaoff[index  ] * cxx[index+2]);
			cxy[index+2] += (cxaoff[index  ] * cxx[index+1] + cxaoff[index+1] * cxx[index  ]);
		}
	}
		break;
//
	case 'c':
	case 'C':
	{
//
// Need to temporarily replace CXX by conjg(CXX)
		for(j1=0; j1<nat3; ++j1)
		{
			cxx[j1] = cxx[j1].conjg();
		}
		GreenFunctionManager::GetInstance()->GetElectric()->Self(cmethd, cxx, gamma, pyd, pzd, akr, akd, dx, cxzw, cxy);
//
// ESELF omitted contribution from diagonal terms: add these. Also remember that ESELF computed -A_jk*conjg(x_k), while
// we seek to compute conjg(A_kj)* x_k.  Since A_kj=A_jk, we have conjg(A_kj)*x_k=conjg(A_jk*conjg(x_k))
		for(j1=0; j1<nat; j1++)
		{
			int index = 3*j1;
			cxy[index  ] = cxadia[index  ] * cxx[index  ] - cxy[index  ];
			cxy[index+1] = cxadia[index+1] * cxx[index+1] - cxy[index+1];
			cxy[index+2] = cxadia[index+2] * cxx[index+2] - cxy[index+2];
		}
//
// Add contribution from off-diagonal elements of 3x3 diagonal blocks
// It is assumed that we have stored off-diagonal elements a_ij as
// CXAOFF(J,1-3)=(a_23,a_31,a_12) for element J
// Note that we could *probably* speed up the code by
// (1) "reducing" product vector CXY immediately after calling ESELF
// (2) restricting CXADIA to occupied sites
// (3) restricting CXAOFF to occupied sites
// but this would require modifications outside this routine, and therefore requires further study...
// 97.12.28 (BTD)
		for(j1=0; j1<nat; ++j1)
		{
			int index = 3*j1;
			cxy[index  ] += (cxaoff[index+1] * cxx[index+2] + cxaoff[index+2] * cxx[index+1]);
			cxy[index+1] += (cxaoff[index+2] * cxx[index  ] + cxaoff[index  ] * cxx[index+2]);
			cxy[index+2] += (cxaoff[index  ] * cxx[index+1] + cxaoff[index+1] * cxx[index  ]);
		}
//
// Now compute conjugate of CXY, and restore CXX:
		for(j1=0; j1<nat3; ++j1)
		{
			cxy[j1] = cxy[j1].conjg();
			cxx[j1] = cxx[j1].conjg();
		}
	}
		break;
//
	default:
		fprintf(stdout, " Error in CPROD: CWHAT= %c\n", cwhat);
		break;
	}
//
// Now must "clean" product vector CXY: (zero out elements corresponding to vacuum sites)
	if (nat0 < nat)
		Nuller(cxy, iocc, nat);
}
