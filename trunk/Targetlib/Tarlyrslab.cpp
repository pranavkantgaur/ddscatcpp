#include "StdAfx.h"

#include "Tarlyrslab.h"
#include "TargetManager.h"

/* **
Routine to construct rectangular prism from "atoms" with layered structure
slab has dimension XV*d   in x direction
                   YV*d   in y direction
                   ZV*d   in z direction
current version allows for up to 4 layers
slab is layered in x direction top surface has x=0
x0(1-3) = point at center of top surface

     0       > x > - f1*a is composition 1
  -f1*a      > x > -(f2+f1)*a is composition 2
-(f2+f1)*a   > x > -(1-f4)*a is composition 3
  -f4*a      > x >  -a   is composition 4

for 1 layers, let f2=f3=f4=0
    2             f3=f4=0
    3             f4=0

It is required that f1+f2+f3+f4=1.

Input:
       XV=(x-length)/d    (d=lattice spacing)
       YV=(y-length)/d    (d=lattice spacing)
       ZV=(z-length)/d    (d=lattice spacing)
       F1,F2,F3,F4 = fraction of thickness contributed by compositions 1,2,3,4
       DX(1-3)=(dx,dy,dz)/d where dx,dy,dz=lattice spacings in x,y,z directions, and d=(dx*dy*dz)**(1/3)=effective lattice spacing
       IOSHP=device number for "target.out" file =-1 to suppress printing of "target.out"
       MXNAT=dimensioning information (max number of atoms)

Output:
       A1(1-3)=unit vector (1,0,0) defining target axis 1 in Target Fr
       A2(1-3)=unit vector (0,1,0) defining target axis 2 in Target Fr
       NAT=number of atoms in target
       IXYZ(1-NAT,1-3)=x/d,y/d,z/d for atoms of target
       ICOMP(1-NAT,1-3)=composition
       X0(3)=location/d in TF of lattice site with IXYZ=(0,0,0)

Fortran history records removed.

Copyright (C) 2007 B.T. Draine and P.J. Flatau
Copyright (c) 2012, C++ version, V.Choliy
This code is covered by the GNU General Public License.
** */

Tarlyrslab::Tarlyrslab(TargetManager *man) : AbstractTarget(man) 
{
	shortDescr = string("Tarlyrslab");
	longDescr = string("Layered slab");
}

void Tarlyrslab::Sizer(void)
{
	dx = manager->CashedDx();
	AllocateArrays(nint_(shpar[0] / dx.data[0]), nint_(shpar[1] / dx.data[1]), nint_(shpar[2] / dx.data[2]));
}

void Tarlyrslab::Descriptor(void)
{
	sprintf(freeDescr, " Layered slab; NX,NY,NZ=%4d%4d%4d\n", nx, ny, nz);
}

void Tarlyrslab::Vector(void)
{
	VectorA();
}

void Tarlyrslab::VectorX(void)
{
//
// Specify X0(3) = (x,y,z)/d in TF corresponding to IXYZ=(0,0,0)
// Top surface of slab (0.5d above topmost dipole layer) is assumed to have x=0
// Dipoles with JX=-1 are assumed to be at x=-dx/2
// Therefore X0(1)=0.5*DX(1)
// Set y=0 and z=0 to be at middle of slab
// JY runs from 1 to NY -> X0(2) = -(1+NY)/2
// JZ runs from 1 to NZ -> X0(3) = -(1+NZ)/2
	x0.Set(half_ * dx.data[0], -half_ * (real)(1 + ny), -half_ * (real)(1 + nz));
}

void Tarlyrslab::Allocator(void)
{
	if (Fabs(shpar[3] + shpar[4] + shpar[5] + shpar[6] - onex_) > (real)1.e-6)
	{
		fprintf(stderr, "f1,f2,f3,f4=%lf%lf%lf%lf\n", shpar[3], shpar[4], shpar[5], shpar[6]);
		Errmsg("Fatal", shortDescr.c_str(), " F1+F2+F3+F4 != 1");
	}

	if ((shpar[3] < zero_) || (shpar[4] < zero_) || (shpar[5] < zero_) || (shpar[6] < zero_))
		Errmsg("Fatal", shortDescr.c_str(), " one or more F are negative");

	int ncomp = 4;
	if (shpar[6] == zero_) ncomp = 3;
	if (shpar[5] == zero_) ncomp = 2;
	if (shpar[4] == zero_) ncomp = 1;

/* **	Original verison code runs wrong
	int nx1 = nint_(shpar[3] * nx) - nx;
	int nx2 = (ncomp >= 2) ? nint_((shpar[3] + shpar[4]) * nx) - nx : 0;
	int nx3 = (ncomp >= 3) ? nint_((shpar[3] + shpar[4] + shpar[5]) * nx) - nx : 0;
** */
	int nx1 = nint_(shpar[3] * nx) - nx;
	int nx2 = (ncomp >= 2) ? nint_((shpar[3] + shpar[4]) * nx) - nx : nx1;
	int nx3 = (ncomp >= 3) ? nint_((shpar[3] + shpar[4] + shpar[5]) * nx) - nx : nx2;
//
// Now populate lattice:
// Top of lattice (in TF) consists of composition 1 as we descend in x (in TF) we pass through compositions 2,3,4
// Top layer will have JX=-1 Lowest layer will have JX=-NX

	minJx = -nx + 1;
	maxJx = 0;
	minJy = 1;
	maxJy = ny;
	minJz = 1;
	maxJz = nz;

	int curSize = nx * ny;
	ixyz.Dimension(curSize, 3);

    int ic = 0;
	for(int jx=-nx+1; jx<=0; ++jx)
	{
		if (jx <= nx1) 
			ic = 1;
		if ((ncomp >= 2) && (jx > nx1) && (jx <= nx2)) 
			ic = 2;
		if ((ncomp >= 3) && (jx > nx2) && (jx <= nx3)) 
			ic = 3;
		if ((ncomp >= 4) && (jx > nx3)) 
			ic = 4;
		for(int jy=1; jy<=ny; ++jy)
		{
			for(int jz=1; jz<=nz; ++jz)
			{
				if (nat0 >= curSize)
				{
					ixyz.Extend(nz);
					curSize = ixyz.GetSize(0);
				}
				ixyz.Fill3(nat0, jx, jy, jz);
				int index = GetLinearAddress(nat0);
				Composer(index, ic);
				iocc[index] = true;
				++nat0;
			}
		}
	}
	icomp.Close(nat0);
}

const char *TargetVerboseDescriptor_Layrdslab(int num)
{

	return NULL;
}

REGISTER_TARGET(Layrdslab,7,false,-1,-6,"rect. block with 2,3, or 4 layers")
void Target_Layrdslab::SayHello(FILE *stream)
{
	fprintf(stream, "LAYRDSLAB = block with 2,4, or 4 layers\n");
	fprintf(stream, "The target %s parameters are:\n", shortDescr.c_str());
	fprintf(stream, "x, y, z length/d:\n");
	fprintf(stream, "%20.16lf %20.16lf %20.16lf\n", shpar[0], shpar[1], shpar[2]);
	fprintf(stream, "fractions f1,f2,f3 (f1+f2+f3.le.1), (for bilayer slab set f1+f2=1, f3=0):\n");
	fprintf(stream, "%20.16lf %20.16lf %20.16lf\n", shpar[3], shpar[4], shpar[5]);
}

const char *TargetVerboseDescriptor_Lyrslbpbc(int num)
{

	return NULL;
}

REGISTER_TARGET(Lyrslbpbc,9,false,7,-6,"TUC = block with max 4 layers")
void Target_Lyrslbpbc::SayHello(FILE *stream)
{
	fprintf(stream, "LYRSLBPBC TUC=block with 4 layers\n");
	fprintf(stream, "The target %s parameters are:\n", shortDescr.c_str());
	fprintf(stream, "x, y, z length/d:\n");
	fprintf(stream, "%20.16lf %20.16lf %20.16lf\n", shpar[0], shpar[1], shpar[2]);
	fprintf(stream, "fractions f1,f2,f3,f4 (f1+f2+f3+f4=1):\n");
	fprintf(stream, "%20.16lf %20.16lf %20.16lf\n", shpar[3], shpar[4], shpar[5]);
	fprintf(stream, "next two parameters s.b. zeros,  %20.16lf %20.16lf:\n", shpar[6], shpar[7]);
}
