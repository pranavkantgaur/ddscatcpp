#include "StdAfx.h"

#include "Definitions.h"
#include "DDscatMain.h"
#include "Matrix.h"

void Evalq(Matrix *theMatrix, Vect3<real> &ak, int nat3, real e02, DipoleData *theDipoleData, real &cabs, real &cext, real &cpha, int imethd)
{
/* **
Given: 
	CXADIA(J,1-3)=(a_11,a_22,a_33) for dipoles J=1-NAT, where symmetric 3x3 matrix a_ij is inverse of Complex polarizability tensor alpha_ij for dipole J
	CXAOFF(J,1-3)=(a_23,a_31,a_12) for dipoles J=1-NAT
	AK(1-3) = (k_x,k_y,k_z)*d d = (d_x*d_y*d_z)**(1/3) = effective lattice spacing
	NAT3 = 3*number of dipoles
	E02 = |E_0|^2 , where E_0 = incident Complex E field.
	CXE(1-NAT3) = components of E field at each dipole, in order E_1x,E_2x,...,E_NATx,E_1y,E_2y,...,E_NATy,E_1z,E_2z,...,E_NATz
	CXP(1-NAT3) = components of polarization vector at each dipole, in order P_1x,P_2x,...,P_NATx,P_1y,P_2y,...,P_NATy,P_1z,P_2z,...,P_NATz
	IMETHD = 0 or 1

Finds:
        CEXT = extinction cross section /d**2 
  and, if IMETHD=1, also computes
        CPHA = phase-lag cross section /d**2
        CABS = absorption cross section /d**2

Note: present definition of CPHA is 1/2 of Martin's definition

B.T.Draine, Princeton Univ. Obs., 87/1/4

History:
Fortran version history removed.

Copyright (C) 1993,1997,1998,2008 B.T. Draine and P.J. Flatau
This code is covered by the GNU General Public License.
** */

//
// Compute magnitude of kd

	Complex *cxe = theDipoleData->Cxe_tf();
	Complex *cxp = theDipoleData->Cxxi();

	int j1;
	int nat = nat3 / 3;
	cext = cabs = cpha = (real)0.;
	real ak2 = ak.ModSquared();
	real ak1 = Sqrt(ak2);
	real ak3 = ak1 * ak2;
//
// Initialization complete.
	if (imethd == 0)
	{
		for(j1=0; j1<nat3; ++j1)
		{
			cext = cext + cxp[j1].im * cxe[j1].re - cxp[j1].re * cxe[j1].im;
		}
		cext = FourPi * ak1 * cext / e02;
	}
	else							// imethd == 1
	{
// C_abs=  (4*pi*k/|E_0|^2)*sum_J { Im[P_J*conjg(a_J*P_J)] - (2/3)*k^3*|P_J|^2 } = (4*pi*k/|E_0|^2)*sum_J {-Im[conjg(P_J)*a_J*P_J] - Im[i*(2/3)*k^3*P_J*conjg(P_J)]}
//      = -(4*pi*k/|E_0|^2)*Im{ sum_J [ conjg(P_J)*(a_J*P_J + i*(2/3)*k^3*P_J) ] }
		Complex *cxadia = theMatrix->Diagonal();
		Complex *cxaoff = theMatrix->OffDiagonal();
		Complex cxa;
		Complex rabs((real)0., ak3 / 1.5);
		for(j1=0; j1<nat; ++j1)
		{
			int index = 3*j1;
			Complex dcxa = cxp[index  ].conjg() * ((cxadia[index  ] + rabs) * cxp[index  ] + cxaoff[index+1] * cxp[index+2] + cxaoff[index+2] * cxp[index+1]) +
				           cxp[index+1].conjg() * ((cxadia[index+1] + rabs) * cxp[index+1] + cxaoff[index+2] * cxp[index  ] + cxaoff[index  ] * cxp[index+2]) + 
				           cxp[index+2].conjg() * ((cxadia[index+2] + rabs) * cxp[index+2] + cxaoff[index  ] * cxp[index+1] + cxaoff[index+1] * cxp[index  ]);
			cxa += dcxa;
		}
		cabs = -FourPi * ak1 * cxa.im / e02;
//
		cxa.clear();
		for(j1=0; j1<nat3; ++j1)
		{
			cxa += cxp[j1] * cxe[j1].conjg();
		}
		cext = FourPi * ak1 * cxa.im / e02;
		cpha =  TwoPi * ak1 * cxa.re / e02;
	}
}
