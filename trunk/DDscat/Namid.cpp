#include "StdAfx.h"

#include "Definitions.h"

void Namid(int myid, char *cfllog)
{
/* **
!***********************************************************************
! Purpose: to generate unique file names for log file written
!          each MPI process

! Present version allows up to 1000 MPI processes (000-999)

! Given:
!         MYID=int (0-999) MPI ID number

! Returns:
!          CFLLOG=name for output file containing running output
!                 from DDSCAT
!                 = ddscat.log_nnn
!                   where nnn = MYID

! B.T.Draine, Princeton Univ. Observatory, 2003.04.12
! History:
! 03.04.12 (BTD): Created using NAMER as an example
! end history

! Copyright (C) 2003, B.T. Draine and P.J. Flatau
! This code is covered by the GNU General Public License.
!***********************************************************************
** */

	sprintf(cfllog, "ddscat.log_%03d", myid);
}
