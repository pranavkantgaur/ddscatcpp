#include "StdAfx.h"

#include "Definitions.h"

void Namer(int iwav, int irad, int idir, char *cflpol1, char *cflpol2, char *cflsca, char *cflavg, 
	char *cfle1, char *cfle2, char *cfleb1, char *cfleb2)
{
/* **
Purpose: to generate file names for specific (wavelength, size, incident direction)
Present version allows up to 1000 wavelengths, (000-999)
                             1000 sizes        (000-999)
                             1000 directions   (000-999)
Given:
        IWAV=int (0-999) identifying wavelength
        IRAD=int (0-999) identifying size
        IDIR=int (0-...) identifying orientation
Returns:
        CFLSCA =name for output file containing Qext,Qabs,Qpha,Qsca,g
                and fml values for selected directions
        CFLAVG =name for output file for given wavelength and size,
                containing scattering etc. averaged over inc.dir.
        CFLPOL1=name for output file for given wavelength and size,
                containing polarization vector for input pol 1
        CFLPOL2=name for output file for given wavelength and size,
                containing polarization vector for input pol 2
        CFLE1  =name for "nearfield" output file with E in rectangular
                volume for input pol 1
        CFLE2  =name for "nearfield" output file with E in rectangular
                volume for input pol 2
        CFLEB1  =name for "nearfield" output file with E and B in rectangular
                volume for input pol 1
        CFLEB2  =name for "nearfield" output file with E and B in rectangular
                volume for input pol 2

B.T.Draine, Princeton Univ. Observatory, 1988

ChB: History records of Fortran versions removed.

Copyright (C) 1993,2005,2006,2007,2011 B.T. Draine and P.J. Flatau
Copyright (C) 2012,2013, C++ versions, Choliy V.
This code is covered by the GNU General Public License.
** */

// *** Filename CFLSCA will be of the form w001r001k001
//      (for IWAV=2,IRAD=2,IDIR=2)

// *** Set IWAV0=IWAV-1 to run over range 0-999
//     Set IRAD0=IRAD-1 to run over range 0-999
//     Set IDIR0=IDIR-1 to run over range 0-999

	sprintf(cflavg,  "w%03dr%03d.avg",       iwav, irad);
	sprintf(cflsca,  "w%03dr%03dk%03d.sca",  iwav, irad, idir);
	sprintf(cflpol1, "w%03dr%03dk%03d.pol1", iwav, irad, idir);
	sprintf(cflpol2, "w%03dr%03dk%03d.pol2", iwav, irad, idir);
	sprintf(cfle1,   "w%03dr%03dk%03d.E1",   iwav, irad, idir);
	sprintf(cfle2,   "w%03dr%03dk%03d.E2",   iwav, irad, idir);
	sprintf(cfleb1,  "w%03dr%03dk%03d.EB1",  iwav, irad, idir);
	sprintf(cfleb2,  "w%03dr%03dk%03d.EB2",  iwav, irad, idir);
}