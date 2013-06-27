// Postprocess.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include "Definitions.h"
#include "Vect3.h"
#include "Complex.h"
#include "Functions.h"
#include "ProcessHelper.h"
#include "Vtr.h"

//
// Purpose: 
// to use subroutine READNF to read data from near-field files written 
// by DDSCAT and then
// postprocess as desired for visualization, etc.
//
// allocatable arrays are passed through modules READNF_ECOM and READNF_BCOM
// other information is passed through READNF argument list
//
// To the user: if desired, add additional processing code
// at bottom of program.  See comments there.
int main(int argc, const char *argv[])
{
/* **
Program DDPOSTPROCESS v2
Purpose: to read binary output files created by nearfield calculation

given (from parameter file ddpostprocess.par)
   CFLENAME = name of file with precomputed P, E, and possibly B
   IVTR     = 0 to skip creation of VTR file
            = 1 to create VTR file with |E|
            = 2 to create VTR file with |E|^2
            = 3 to create VTR file with time-average <Re(E)xRe(B)>/(4*pi)
   ILINE    = 1 to evaluate E field at points along 1 or more lines
   XA,YA,ZA = (x,y,z)_TF (physical units) for starting point of line
   XB,YB,ZB = (x,y,z)_TF (physical units) for endpoint of line
   NAB      = number of points along line, including A and B

extracts from file (by using subroutine READNF):
   AEFF             = effective radius of target (phys. units)
   NAMBIENT         = (real) refractive index of ambient medium
   WAVE             = wavelength in vacuo of incident wave (phys. units)
   DPHYS            = interdipole separation (phys. units)
   NAT0             = number of dipoles in physical target
   NCOMP            = number of distinct compositions present
   NX,NY,NZ         = dimensions/d of computational volume (computational volume has NXYZ=NX*NY*NZ points)
   X0(1-3)          = (x/d,y/d,z/d) in Target frame for index I,J,K=0,0,0, thus (x,y,z)=[ x0(1,2,3) + (I,J,K) ]*d
   AK_TF(1-3)       = (k_x,k_y,k_z)*d in the Target Frame
   CXE0_TF(1-3)     = E_inc (complex) in the Target Frame at (x_TF,y_TF,z_TF)=(0,0,0)
   ICOMP(1-3*NXYZ)  = composition identifier for all points and directions
   CXEINC(1-NXYZ,3) = complex incident macroscopic E field at all points
   CXESCA(1-NXYZ,3) = complex radiated macroscopic E field at all points
   CXPOL(1-NXYZ,3)  = complex polarization/d^3 at all points
   CXADIA(1-NXYZ,3) = diagonal element of polarizability/d^3 at all pts

if the stored file contained magnetic field information (NRFLDB=1) then also return

   CXBINC(1-NXYZ,3) = complex incident B field at all points
   CXBSCA(1-NXYZ,3) = complex radiated B field at all points
   ICOMP(1-3*NXYZ)  = composition identifier at all points (= 0 for vacuum)

using extracted information, ddpostprocess uses simple interpolation to 
evaluate, for each of NAB points on line:

   CXE_INC(1-3)= (complex) incident E field at location (x_TF,y_TF,z_TF)
   CXE_SCA(1-3)= (complex) radiated E field at location (x_TF,y_TF,z_TF)
   CXP(1-3)    = (complex) polarization/d^3 at location (x_TF,y_TF,z_TF)
 
and, if NRFLDB=1:

   CXB_INC(1-3) = (complex) incident B field at location (x_TF,y_TF,z_TF)
   CXB_SCA(1-3) = (complex) radiated B field at location (x_TF,y_TF,z_TF)

current version writes out textfile ddpostprocess_E.out with one line per point along track:

   x_TF,y_TF,z_TF,CXE(1-3)           [if NRFLDB=0: only E field is available]
   x_TF,y_TF,z_TF,CXE(1-3),CXB(1-3)  [if NRFLDB=1: both E and B are available]

where CXE=CXE_INC+CXE_SCA = total macroscopic E field at point
      CXB=CXB_INC+CXB_SCA = total B field at point

If IVTR > 0, then generate VTR output files

*** NB: Existing code simply writes out E and B at points along a 
        straight line as a simple example with limited output.

       Users who wish to write out other information for purposes
       of display or analysis should go to the end of the existing
       program and add additional code to write out whatever is desired
       (e.g., you may want E or B at points on a 2-D plane, or 3-D volume).
          
DDPOSTPROCESS is adapted from program originally first named READE, 
then renamed READNF.  READE was first written 2011.08.30
B.T. Draine, Princeton University Observatory

Fortran versions history records removed.

Copyright (C) 2011,2012,2013 B.T. Draine and P.J. Flatau
Copyright (C) 2013, C++ version, V.Choliy
This code is covered by the GNU General Public License.
** */

	const char *cflpar_default = "DDpostprocess.par";
	char cflpar[20];
	strcpy(cflpar, cflpar_default);
	if(argc == 2)
	{
		strcpy(cflpar, argv[1]);
	}
	printf(">DDpostprocess using parameter file = %s\n", cflpar);
//
// NRWORD = length (bytes) of REAL word used to check compatibility with word length in stored data
//	const int nrword = sizeof(real);
//
	char Buffer[256], cflename[64], cflvtr[64];
	int ivtr, iline;
	FILE *paramFile = fopen(cflpar, "r");
	if (paramFile == NULL)
	{
		printf("Error opening Readnf parameter file.\n");
		return -1;
	}
	fgets(Buffer, 255, paramFile);
	char *ia = strchr(Buffer+1, '\'');
	if (ia) *ia = '\0';
	strcpy(cflename, Buffer+1);
	printf(">DDpostprocess input data from %s\n", cflename);
	fgets(Buffer, 255, paramFile);
	ia = strchr(Buffer+1, '\'');
	if(ia) *ia = '\0';
	strcpy(cflvtr, Buffer+1);
	printf(">DDpostprocess Vtk output name %s\n", cflvtr);
	fgets(Buffer, 255, paramFile);
	sscanf(Buffer, "%d", &ivtr);
//
// If IVTR > 0, then create VTR files for subsequent visualization
	fgets(Buffer, 255, paramFile);
	sscanf(Buffer, "%d", &iline);
	printf(">DDpostprocess now read file = %s\n", cflename);
//
	ProcessHelper *helper = new ProcessHelper;
	bool bOk = helper->Load(cflpar, sizeof(real));
	if (bOk == false)
	{
		fprintf(stderr, "Error opening binaryFile in DDpostprocess::Load\n");
		delete helper;
		return -1;
	}
	printf(">DDpostprocess text Version of input file %s\n", helper->GetStringVersion());
	printf(">DDpostprocess numeric Version of input file %d\n", helper->GetIntVersion());

// determine E_inc, E_sca, and P at points along defined track

// if NRFLDB=1, also calculate B and time-averaged S along track, where
// S=Poynting vector normalized by Poynting vector of incident plane wave
// Let E=Re[Ec*e^(-iwt)]=E1cos(wt)+E2sin(wt)  where Ec=E1+iE2
//     B=Re[Bc*e^(-iwt)]=B1cos(wt)+B2sin(wt)  where Bc=B1+iB2
//     S=(1/4pi)*ExB
//    <S>=(1/4pi)*(1/2)(E1xB1+E2xB2)
//       =(1/8pi)*( E1xB1 + E2xB2 )
//       =(1/8pi)*Re( Ec x conjg(Bc) )
// we will omit the (1/8pi) because we will always normalize by incident S

// calculate |<S>| for incident wave
	helper->PrepareNormsAndPointing();
	printf("DDpostprocess 8*pi*|<s>| for incident wave = %11.3lf\n", helper->GetSnorm());
//
	sprintf(Buffer, "Ncomp = %d compositions\n", helper->GetNcomp());
	Wrimsg("DDpostprocess", Buffer);
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
//
// have completed:
// 1. reading data from file
// 2. evaluating E_inc, E_sca, and P at points along track
// 3. writing track data to ascii output file

// *** Additional code below here to write out additional files for display or analysis.
	if (ivtr <= 0)
	{
		delete helper;
		return -3;
	}

// if IVTR > 0: now write VTK file for graphics define mesh for graphics

// VTK supplementary arrays 
// mesh x,y,z assuming that we are on rectangular grid)
// variables are dimensioned nx,ny,nz

	int jj1 = helper->SanityCheck();
	if (jj1 != helper->GetNat0())
	{
		fprintf(stderr, "Sanity failure: inconsistent j1=%d and nat0=%d\n", jj1, helper->GetNat0());
		delete helper;
		return -1;
	}
//
// Write information to unit IDVOUT (usually a log file)
	printf(">DDpostprocess %11.4e = normalized error |P/alpha-E|^2/|E_inc|^2\n", helper->GetSumerr2());
	printf(">DDpostprocess %11.4e = AEFF (vol. equiv. radius, phys. units)\n", helper->GetAeff());
	printf(">DDpostprocess %11d  = NAT0 (number of physical dipoles in target)\n", helper->GetNat0());
	printf(">DDpostprocess %11.4e = d = interdipole separation (phys. units)\n", helper->GetDphys());
	printf(">DDpostprocess %11.4e = wavelength in vacuo (phys. units)\n", helper->GetWave());
	printf(">DDpostprocess %11.4e = wavelength in ambient medium (phys. units)\n", helper->GetWave() / helper->GetNambient());
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
		printf(">DDpostprocess completed writing VTK file\n");
	}
	else
        return -3;

	delete helper;

	return 0;
}
