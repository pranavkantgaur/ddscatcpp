#include "StdAfx.h"

#include "FileNamer.h"
#include "Definitions.h"

FileNamer *FileNamer::item = NULL;
FileNamer::FileNamer(void)
{
	norichar = 0;
	cashedIwav = cashedIrad = cashedIdir = -1;
}

FileNamer::~FileNamer(void)
{

}

FileNamer *FileNamer::GetInstance()
{
	if (!item)
		item = new FileNamer;

	return item;
}

void FileNamer::Kill(void)
{
	CleanDelete(item);
}

void FileNamer::Init(int nr)
{
	norichar = 1 + (int)log10((real)nr);
}

void FileNamer::Namer(int iwav, int irad, int idir)				// and Namer2
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
        CFLSCA =name for output file containing Qext,Qabs,Qpha,Qsca,g and fml values for selected directions
        CFLAVG =name for output file for given wavelength and size, containing scattering etc. averaged over inc.dir.
        CFLPOL1=name for output file for given wavelength and size, containing polarization vector for input pol 1
        CFLPOL2=name for output file for given wavelength and size, containing polarization vector for input pol 2
        CFLE1  =name for "nearfield" output file with E in rectangular volume for input pol 1
        CFLE2  =name for "nearfield" output file with E in rectangular volume for input pol 2
        CFLEB1  =name for "nearfield" output file with E and B in rectangular volume for input pol 1
        CFLEB2  =name for "nearfield" output file with E and B in rectangular volume for input pol 2
        CFLFML =name for output file containing Complex scattering amplitudes f_ml for selected directions

B.T.Draine, Princeton Univ. Observatory, 1988

History records of Fortran versions removed.

Copyright (C) 1993,2005,2006,2007,2011,2012,2013 B.T. Draine and P.J. Flatau
Copyright (C) 2012,2013, C++ versions, Choliy V.
This code is covered by the GNU General Public License.
** */

// *** Filename CFLSCA will be of the form w001r001k001
//      (for IWAV=2,IRAD=2,IDIR=2)

// *** Set IWAV0=IWAV-1 to run over range 0-999
//     Set IRAD0=IRAD-1 to run over range 0-999
//     Set IDIR0=IDIR-1 to run over range 0-999

	if ((iwav == cashedIwav) && (irad == cashedIrad) && (idir == cashedIdir))
		return;

	char Format[32], Buffer[128];
	sprintf(Format, "w%%03dr%%03dk%%0%dd.%%s", max(norichar, 3));

    sprintf(Buffer, "w%03dr%03d.avg", iwav, irad);		cflavg = string(Buffer);
	sprintf(Buffer, Format, iwav, irad, idir, "sca");	cflsca = string(Buffer);
	sprintf(Buffer, Format, iwav, irad, idir, "pol1");	cflpol1 = string(Buffer);
	sprintf(Buffer, Format, iwav, irad, idir, "pol2");	cflpol2 = string(Buffer);
	sprintf(Buffer, Format, iwav, irad, idir, "E1");	cfle1 = string(Buffer);
	sprintf(Buffer, Format, iwav, irad, idir, "E2");	cfle2 = string(Buffer);
	sprintf(Buffer, Format, iwav, irad, idir, "EB1");	cfleb1 = string(Buffer);
	sprintf(Buffer, Format, iwav, irad, idir, "EB2");	cfleb2 = string(Buffer);
	sprintf(Buffer, Format, iwav, irad, idir, "fml");	cflfml = string(Buffer);	

	cashedIwav = iwav;
	cashedIrad = irad;
	cashedIdir = idir;
}
