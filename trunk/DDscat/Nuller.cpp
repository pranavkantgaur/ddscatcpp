#include "StdAfx.h"

#include "Definitions.h"
#include "Complex.h"

void Nuller(Complex *cxvec, bool *iocc, int nat)
{
/* **
!***********************************************************************
! This routine takes as input a Complex vector with NAT3 elements
! and sets to zero those elements corresponding to "vacuum" sites.
! This is accomplished by multiplying by vector IOCC(J), whose
! elements are either 1 or 0 depending on whether site is occupied
! or "vacuum".
! Note: CXVEC should be given dimension CXVEC(NAT,3) here
!       so that first 3*NAT elements of CXVEC are employed.

! B.T.Draine, Princeton Univ. Obs., 90/11/1
! History:
! 90/12/15 (BTD): Corrected error in dimensioning of CXVEC.

! Copyright (C) 1993, B.T. Draine and P.J. Flatau
! This code is covered by the GNU General Public License.
!***********************************************************************
** */
	for(int j=0; j<nat; ++j)
	{
		if (iocc[j] == false)
		{
			cxvec[3*j  ].clear();
			cxvec[3*j+1].clear();
			cxvec[3*j+2].clear();
		}
	}
}