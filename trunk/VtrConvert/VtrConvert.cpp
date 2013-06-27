// VtrConvert.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include "Vtr.h"
#include "Vect3.h"
#include "ArrayF.h"

int main(int argc, const char *argv[])
{
/* **
Purpose:
This code converts DDSCAT target shape file to VTK format

Calling sequence:
vtrconvert input_ddscat_shape output_vtr_file
where output_vtr_file is an output file name (without prefix) input_ddscat_shape is file created by "calltarget"

Thus typical sequance would be:
calltarget < sphere40x40x40.shp" (to create target.out file - DDSCAT shape file)
vtrtarget,"target.out","sphere40x40x40" (to do the conversion)
where the input file "sphere40x40x40.shp"  could be
ELLIPSOID
40 40 40
0 0 0

To compile the VTRCONVER:
"gfortran vtr.f90 vtrconvert.f90 -o vtrconvert"

NOTES:
You can use "paraview", "mayavi2" to visualize the data using VTK format

History:
Written by PJF 2011
(PJF) Feb 2012 Added capability to output a1, a2 vectors to separate file
Copyright (C) 2011 B.T. Draine and P.J. Flatau
This code is covered by the GNU General Public License.
** */

	Vtr vtr;					// VtrFileHandle *fd = vtr.FileHandle();
	char cddscat[81], cvtr[81];
	strcpy(cddscat, "target.out");
	strcpy(cvtr, "target");
	const int ipad = 5;
	if(argc == 3)
	{
		strcpy(cddscat, argv[1]);
		strcpy(cvtr, argv[2]);
		printf(" VTRTARGET>> we will use input DDSCAT file shape format name: %s\n", cddscat);
		printf(" VTRTARGET>> we will output VTK name: %s\n", cvtr);
	}
	else
	{
		printf(" VTRTARGET>> you have to define 2 files on the command line output\n");
		printf(" VTRTARGET>> vtrtarget  ddscat_file vtk_file (without extension)\n");
		printf(" VTRTARGET>> for example\n");
		printf(" VTRTARGET>> vtrtarget target.out  output\n");
		printf(" exiting now\n");
		return -1;
	}
//
	FILE *ioshpFile = fopen(cddscat, "r");
	if (!ioshpFile)
	{
		printf("Cannot open ioshpFile.\n");
		return -2;
	}
//
	char Buffer[256], cdescr[81];
	Vect3<real> a1, a2, dx, x0;
	int nat, jxx, icomp1, icomp2, icomp3, ixyz1, ixyz2, ixyz3;

	fgets(Buffer, 255, ioshpFile);
	fgets(Buffer, 255, ioshpFile);
	sscanf(Buffer, "%d", &nat);
	fgets(Buffer, 255, ioshpFile);
	a1.Load(Buffer, realFormat);	
	fgets(Buffer, 255, ioshpFile);
	a2.Load(Buffer, realFormat);
	fgets(Buffer, 255, ioshpFile);
	dx.Load(Buffer, realFormat);       
	fgets(Buffer, 255, ioshpFile);
	x0.Load(Buffer, realFormat);	
	fgets(cdescr, 255, ioshpFile);
	long position = ftell(ioshpFile);
	fgets(Buffer, 255, ioshpFile);
	sscanf(Buffer, "%d%d%d%d%d%d%d", &jxx, &ixyz1, &ixyz2, &ixyz3, &icomp1, &icomp2, &icomp3);

	int mi1 = ixyz1;
	int mi2 = ixyz2;
	int mi3 = ixyz3;
	int mx1 = ixyz1;
	int mx2 = ixyz2;
	int mx3 = ixyz3;
	for(int jx=1; jx<nat; ++jx)
	{
		fgets(Buffer, 255, ioshpFile);
		sscanf(Buffer, "%d%d%d%d%d%d%d", &jxx, &ixyz1, &ixyz2, &ixyz3, &icomp1, &icomp2, &icomp3);
		if (ixyz1 < mi1) mi1 = ixyz1;
		if (ixyz2 < mi2) mi2 = ixyz2;
		if (ixyz3 < mi3) mi3 = ixyz3;
		if (ixyz1 > mx1) mx1 = ixyz1;
		if (ixyz2 > mx2) mx2 = ixyz2;
		if (ixyz3 > mx3) mx3 = ixyz3;
	}
	printf("%d %d %d %d %d %d\n", mi1, mi2, mi3, mx1, mx2, mx3);
	fseek(ioshpFile, position, SEEK_SET);
//
// Allocate array and write
// I am padding to allocate a bit of space around the target so the graphics codes do not make objects "holow" (check ipad=0 to see what I am talking about)
	int mi1ext = mi1 - ipad;
	int mi2ext = mi2 - ipad;
	int mi3ext = mi3 - ipad;

	int mx1ext = mx1 + ipad;
	int mx2ext = mx2 + ipad;
	int mx3ext = mx3 + ipad;

	Array3F<real> vtr8;
	vtr8.Dimension(mx1ext - mi1ext + 1, mx2ext - mi2ext + 1, mx3ext - mi3ext + 1);
	real *x = new real[mx1ext-mi1ext+1];
	real *y = new real[mx2ext-mi2ext+1];
	real *z = new real[mx3ext-mi3ext+1];
//
	for(int jx=0; jx<nat; ++jx)
	{
		fgets(Buffer, 255, ioshpFile);
		sscanf(Buffer, "%d%d%d%d%d%d%d", &jxx, &ixyz1, &ixyz2, &ixyz3, &icomp1, &icomp2, &icomp3);
		vtr8.Value(ixyz1-1-mi1ext, ixyz2-1-mi2ext, ixyz3-1-mi3ext) = icomp1;
	}
	fclose(ioshpFile);
//
	for(int jx=mi1ext; jx<=mx1ext; ++jx)
	{
		x[jx-mi1ext] = (real)jx;
	}
	for(int jx=mi2ext; jx<=mx2ext; ++jx)
	{
		y[jx-mi2ext] = (real)jx;
	}
	for(int jx=mi3ext; jx<=mx3ext; ++jx)
	{
		z[jx-mi3ext] = (real)jx;
	}
//
// Copy to kind=8 precision for graphics
	vtr.OpenFile(cvtr);
	vtr.WriteMesh3d(x, y, z, mx1ext-mi1ext+1, mx2ext-mi2ext+1, mx3ext-mi3ext+1);
	vtr.WriteScalar3d("Composition", vtr8);
	vtr.CloseFile();
	vtr.CollectFile();
//
// Write vectors  a1, a2 to separate files copy to kind=8 precision for graphics
	vtr.OpenFile("a1a2");
	real xvect = (real)0.;
	real yvect = (real)0.;
	real zvect = (real)0.;
	vtr.WriteMesh3d(&xvect, &yvect, &zvect, 1, 1, 1);
//
	Array3F<real> u, v, w;
	u.Dimension(1, 1, 1);
	v.Dimension(1, 1, 1);
	w.Dimension(1, 1, 1);
	u.Value(0,0,0) = a1.data[0];
	v.Value(0,0,0) = a1.data[1];
	w.Value(0,0,0) = a1.data[2];
	vtr.WriteVector3d("a1", u, v, w, 1, 1, 1);
	u.Value(0,0,0) = a2.data[0];
	v.Value(0,0,0) = a2.data[1];
	w.Value(0,0,0) = a2.data[2];
	vtr.WriteVector3d("a2", u, v, w, 1, 1, 1);
	vtr.CloseFile();
	vtr.CollectFile();
	
	return 0;
}
