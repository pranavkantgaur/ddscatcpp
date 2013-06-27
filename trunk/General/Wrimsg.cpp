#include "StdAfx.h"

void Wrimsg(const char *csubrt, const char *cmsgnm)
{
/* **
! Standard procedure for writing messages
! History:
! 96.11.14 (PJF) Remove "getset" and hardwire" ioout"
! 96.11.20 (BTD) change IOOUT to IDVOUT

! Copyright (C) 1993,1996 B.T. Draine and P.J. Flatau
! This code is covered by the GNU General Public License.
!***********************************************************************
** */

// ! Note: on some systems (e.g., Solaris) IDVOUT=0 generates unbuffered
// ! output to "standard output".  On other systems IDVOUT=0 may not be
// ! valid; then set IDVOUT=6 to get "standard output", probably buffered.

	const char *Format9000 = " >%s %s\n";
	fprintf(stdout, Format9000, csubrt, cmsgnm);
}