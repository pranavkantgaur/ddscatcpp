#include "StdAfx.h"

#include "Complex.h"
#include "DDscatMain.h"

void Cxfft3_mkl(Array3F<Complex> &cx, int mx, int my, int mz, int isign)
{
/* **
!=======================================================================
! Purpose: This is a dummy routine to substitute for CXFFT3_MKL to
!          allow DDSCAT to be used on systems where the Intel MKL
!          library is not available.  If called, it will generate
!          a fatal error with an explanatory error message.
! History
! 08.06.05 (BTD): created
! end history
! Copyright (C) 2008
!               B.T. Draine and P.J. Flatau
! This code is covered by the GNU General Public License.
!=======================================================================
** */

	char cmsgnm[72];

	sprintf(cmsgnm, "FATAL ERROR: DDSCAT compiled with cxfft3_mkl_fake");
	Wrimsg("Cxfft3_mkl", cmsgnm);
	sprintf(cmsgnm, " *** option FFTMKL cannot be used with dummy routine");
	Wrimsg("Cxfft3_mkl", cmsgnm);
	sprintf(cmsgnm, " *** to enable option FFTMKL it is necessary to:");
	Wrimsg("Cxfft3_mkl", cmsgnm);
	sprintf(cmsgnm, "     have Intel Math Kernel Library (MKL) installed on system");
	Wrimsg("Cxfft3_mkl", cmsgnm);
	sprintf(cmsgnm, "     and edit Makefile to use cxfft3_mkl.f90 and mkl_dfti.f90");
	Wrimsg("Cxfft3_mkl", cmsgnm);
	exit(0);
}
