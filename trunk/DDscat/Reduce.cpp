#include "StdAfx.h"

#include "Definitions.h"
#include "Vect3.h"
#include "Complex.h"

void Reduce(Complex *cxv, bool *iocc, int nat, int nat0)
{
/* **
!***********************************************************************
! Given:
!      CXV(1-3*NAT) defined for NAT lattice sites with ordering
!          (v_x1,...,v_xj,...,v_xNAT,v_y1,...,v_yj,...,v_yNAT,
!           v_z1,...,v_zj,...,v_zNAT)
!          where index j runs over occupied and unoccupied lattice sites
!      IOCC(1-NAT) = 0 or 1 depending on whether lattice site
!                    is vacant or occupied
!      MXN3,MXNAT = dimensioning information
!      NAT = number of lattice sites
!      NAT0 = number of occupied lattice sites

! Returns
!      CXV(1-3*NAT0) defined for NAT0 occupied lattice sites
!         with ordering
!         (v_x1,...,v_xj,...,v_xNAT0,v_y1,...,v_yj,...,v_yNAT0,
!          v_z1,...,v_zj,...,v_zNAT0)
!         where index j now runs over only occupied lattice sites

! B.T.Draine, Princeton Univ. Observatory, 90.11.29
! History:
! 90.12.03 (BTD): Modified for new ordering of vectors.
! 90.12.05 (BTD): Corrected errors
! 06.04.11 (BTD): Cosmetic changes
! end history

! Copyright (C) 1993, 2006 B.T. Draine and P.J. Flatau
! This code is covered by the GNU General Public License.
!***********************************************************************
** */

	if (nat == nat0)
		return;

	int joc = 0;
	for(int j=0; j<nat; ++j)
	{
		if (iocc[j] == true)
		{
			cxv[joc  ] = cxv[3*j];
			cxv[joc+1] = cxv[3*j+1];
			cxv[joc+2] = cxv[3*j+2];
			joc += 3;
		}
	}
}
