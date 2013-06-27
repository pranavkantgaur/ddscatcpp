// Readnf1.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include "Vect3.h"
#include "Complex.h"
#include "Functions.h"
#include "ArrayF.h"
#include "ProcessHelper.h"

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
	const char *cflpar_default = "Readnf1.par";
//
	char cflpar[64];
	strcpy(cflpar, cflpar_default);
	if(argc == 2)
	{
		strcpy(cflpar, argv[1]);
	}
	printf(">Readnf1 using parameter file = %s\n", cflpar);
//
	char Buffer[256], cflename[64], cflvtr[64];
	int ivtr;
	FILE *paramFile = fopen(cflpar, "r");
	if (paramFile == NULL)
	{
		printf("Error opening Readnf1 parameter file.\n");
		return -1;
	}
	fgets(Buffer, 255, paramFile);
	ExtractFirstWord(Buffer, '\'');
	strcpy(cflename, Buffer);
	printf(">Readnf1 input data from %s\n", cflename);
//
	fgets(Buffer, 255, paramFile);
	ExtractFirstWord(Buffer, '\'');
	strcpy(cflvtr, Buffer);
	printf(">Readnf1 Vtk output name %s\n", cflvtr);
//
	fgets(Buffer, 255, paramFile);
	sscanf(Buffer, "%d", &ivtr);
	if ((ivtr < 1) || (ivtr > 3))
	{
		printf("Readnf1 wrong ivtr = %d option.\n", ivtr);
		return -1;
	}
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
		fprintf(stderr, "Error opening binaryFile in Readnf1::Load\n");
		delete helper;
		return -1;
	}
	printf(">Readnf1 text Version of input file %s\n", helper->GetStringVersion());
	printf(">Readnf1 numeric Version of input file %d\n", helper->GetIntVersion());
//
	helper->PrepareNormsAndPointing();
	printf(">Readnf1 8*pi*|<s>| for incident wave = %11.3lf\n", helper->GetSnorm());
//
	int jj1 = helper->SanityCheck();
	if (jj1 != helper->GetNat0())
	{
		fprintf(stderr, "Sanity failure: inconsistent j1=%d and nat0=%d\n", jj1, helper->GetNat0());
		delete helper;
		return -1;
	}
//
// Write information to unit IDVOUT (usually a log file)
	printf(">Readnf1 %11.4e = normalized error |P/alpha-E|^2/|E_inc|^2\n", helper->GetSumerr2());
	printf(">Readnf1 %11.4e = AEFF (vol. equiv. radius, phys. units)\n", helper->GetAeff());
	printf(">Readnf1 %11d  = NAT0 (number of physical dipoles in target)\n", helper->GetNat0());
	printf(">Readnf1 %11.4e = d = interdipole separation (phys. units)\n", helper->GetDphys());
	printf(">Readnf1 %11.4e = wavelength in vacuo (phys. units)\n", helper->GetWave());
	printf(">Readnf1 %11.4e = wavelength in ambient medium (phys. units)\n", helper->GetWave() / helper->GetNambient());
//
// have completed:
// 1. reading data from file
// 2. evaluating E_inc, E_sca, and P at points along track
// 3. writing track data to ascii output file

// if IVTR > 0: now write VTK file for graphics define mesh for graphics
	if (ivtr > 0)
	{
		helper->PrepareVTKFile(cflvtr, ivtr);
// Bruce: in the future we should output CXE itself as a vector.
// I am not doing it now, so we can start
		printf(">Readnf1 completed writing VTK file\n");
	}
	else
        return -3;

	delete helper;

	return 0;
}
