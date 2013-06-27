// Readnf.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include "Vtr.h"
#include "Vect3.h"
#include "Complex.h"
#include "Functions.h"
#include "ArrayF.h"
#include "ProcessHelper.h"
#include "Linia.h"

int main(int argc, const char *argv[])
{
/* **
Program READNF
Purpose: to read binary output files created by nearfield calculation

given (from parameter file readnf.par)
   CFLENAME = name of file with precomputed P, E, and possibly B
   IVTR     = 0 to skip creation of VTR file
            = 1 to create VTR file with |E|
            = 2 to create VTR file with |E|^2
   ILINE    = 1 to evaluate E field at points along 1 or more lines
   XA,YA,ZA = (x,y,z)_TF (physical units) for starting point of line
   XB,YB,ZB = (x,y,z)_TF (physical units) for endpoint of line
   NAB      = number of points along line, including A and B

extracts from file:
   AEFF             = effective radius of target (phys. units)
   NAMBIENT         = (real) refractive index of ambient medium
   WAVE             = wavelength in vacuo of incident wave (phys. units)
   DPHYS            = interdipole separation (phys. units)
   NAT0             = number of dipoles in physical target
   NCOMP            = number of distinct compositions present
   NX,NY,NZ         = dimensions/d of computational volume
                      (computational volume has NXYZ=NX*NY*NZ points)
   X0(1-3)          = (x/d,y/d,z/d) in Target frame for index I,J,K=0,0,0
                      thus (x,y,z)=[ x0(1,2,3) + (I,J,K) ]*d
   AK_TF(1-3)         = (k_x,k_y,k_z)*d in the Target Frame
   CXE0_TF(1-3)       = E_inc (complex) in the Target Frame at (x_TF,y_TF,z_TF)=(0,0,0)
   ICOMP(1-3*NXYZ)  = composition identifier for all points and directions
   CXEINC(1-3*NXYZ) = complex incident macroscopic E field at all points
   CXESCA(1-3*NXYZ) = complex radiated macroscopic E field at all points
   CXPOL(1-3*NXYZ)  = complex polarization/d^3 at all points
   CXADIA(1-3*NXYZ) = diagonal element of polarizability/d^3 at all pts
   CXBINC(1-3*NXYZ) = complex incident B field at all points
   CXBSCA(1-3*NXYZ) = complex radiated B field at all points
   ICOMP(1-3*NXYZ)  = composition identifier at all points (= 0 for vacuum)

using extracted information, readnf uses interpolation to evaluate,
for each of NAB points on line:

   CXE_INC(1-3)= (complex) incident E field at location (x_TF,y_TF,z_TF)
   CXE_SCA(1-3)= (complex) radiated E field at location (x_TF,y_TF,z_TF)
   CXP(1-3)    = (complex) polarization/d^3 at location (x_TF,y_TF,z_TF)
 
and, if NRFLDB=1:

   CXB_INC(1-3) = (complex) incident B field at location (x_TF,y_TF,z_TF)
   CXB_SCA(1-3) = (complex) radiated B field at location (x_TF,y_TF,z_TF)

current version writes out textfile readnf_E.out with one line per point along track:

   x_TF,y_TF,z_TF,CXE(1-3)           [if NRFLDB=0: only E field is available]
   x_TF,y_TF,z_TF,CXE(1-3),CXB(1-3)  [if NRFLDB=1: both E and B are available]

where CXE=CXE_INC+CXE_SCA = total macroscopic E field at point
      CXB=CXB_INC+CXB_SCA = total B field at point

*** Note: Existing code writes out E and B at points along a straight line
          as a simple example with output of limited volume.
          Users who wish to write out other information for purposes
          of display or analysis should go to the end of the existing
          program and modify it to write out whatever is desired
          (e.g., you may want E or B at points on a 2-D plane, or 3-D volume).
          
          
If IVTR > 0, then generate VTR output files

B.T. Draine, Princeton Univ. Observatory, 2011.08.30

History records of Fortran versions removed.

Copyright (C) 2011,2012,2013 B.T. Draine and P.J. Flatau
Copyright (C) 2012,2013 C++ versions, Choliy V.

This code is covered by the GNU General Public License.
** */
	const int nrword = sizeof(real);						// NRWORD = length (bytes) of REAL word used to check compatibility with word length in stored data
	const char *cflpar_default = "Readnf2.par";
//
	char cflpar[64];
	strcpy(cflpar, cflpar_default);
	if(argc == 2)
	{
		strcpy(cflpar, argv[1]);
	}
	printf(">Readnf2 using parameter file = %s\n", cflpar);
//
	char Buffer[256], cflename[64];
	int iline;
	FILE *paramFile = fopen(cflpar, "r");
	if (paramFile == NULL)
	{
		printf("Error opening Readnf2 parameter file.\n");
		return -1;
	}
	fgets(Buffer, 255, paramFile);
	ExtractFirstWord(Buffer, '\'');
	strcpy(cflename, Buffer);
	printf(">Readnf2 input data from %s\n", cflename);
//
// Read iline
	fgets(Buffer, 255, paramFile);
	sscanf(Buffer, "%d", &iline);
	printf(">Readnf2 now will calculate %d lines\n", iline);
//
//              >>>>> Important Note! <<<<<
// The structure of the READ statements below *must* conform to the
// structure of the corresponding WRITE statements in nearfield.f90 
// Any changes must be made in both modules.
//
	ProcessHelper *helper = new ProcessHelper;
	bool bOk = helper->Load(cflename, nrword);
	if (bOk == false)
	{
		fprintf(stderr, "Error opening binaryFile in Readnf1::Load, nrwords are: %d & %d\n", nrword, helper->GetNrwordNf());
		delete helper;
		return -1;
	}
	printf(">Readnf2 text Version of input file %s\n", helper->GetStringVersion());
	printf(">Readnf2 numeric Version of input file %d\n", helper->GetIntVersion());
//
	helper->PrepareNormsAndPointing();
	printf("Readnf2 8*pi*|<s>| for incident wave = %11.3lf\n", helper->GetSnorm());
//
	if (iline > 0)
	{
		int nline = 0;
		while(1)
		{
			Linia linia;
			char *ia = fgets(Buffer, 255, paramFile);
			if (!ia && feof(paramFile))
				break;
			linia.Scanf(Buffer);
			if (helper->IsLiniaOk(linia))
			{
				sprintf(Buffer, "Readnf2_%03d.out",  ++nline);
				helper->WriteLiniaData(Buffer, linia);
			}
			else
			{
				fprintf(stderr, "Fatal error: requested track extends beyond computational volume\n");
				helper->OutMinMax();
				linia.Debug();
				fprintf(stderr, "This line is skipped.\n");		
			}
		}
	}
	delete helper;

	return 0;
}
