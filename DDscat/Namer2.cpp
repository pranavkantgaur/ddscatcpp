#include "StdAfx.h"

#include "Definitions.h"

void Namer2(int iwav, int irad, int idir, char *cflfml)
{
/* **
!***********************************************************************
! Purpose: to generate file names for specific (wavelength,
!                                               size,
!                                               incident direction)
! Present version allows up to 1000 wavelengths, (000-999)
!                              1000 sizes        (000-999)
!                              1000 directions   (000-999)
! Given:
!         IWAV=int (0-999) identifying wavelength
!         IRAD=int (0-999) identifying size
!         IDIR=int (0-999) identifying direction
! Returns:
!         CFLFML =name for output file containing Complex scattering
!                 amplitudes f_ml for selected directions
!
! B.T.Draine, Princeton Univ. Observatory, 2007
! History:
! 07.08.31 (BTD): Adapted from NAMER
! Copyright (C) 2007 B.T. Draine and P.J. Flatau
! This code is covered by the GNU General Public License.
!***********************************************************************
** */

	sprintf(cflfml, "w%03dr%03dk%03d.fml", iwav, irad, idir);
}