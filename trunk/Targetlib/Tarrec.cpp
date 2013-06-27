#include "StdAfx.h"

#include "Tarrec.h"
#include "TargetManager.h"

/* **
Routine to construct rectangular prism from "atoms"
Input:
       XV=(x-length)/d    (d=lattice spacing)
       YV=(y-length)/d
       ZV=(z-length)/d
       DX(1-3)=(dx,dy,dz)/d where dx,dy,dz=lattice spacings in x,y,z
               directions, and d=(dx*dy*dz)**(1/3)=effective lattice spacing
       IOSHP=device number for "target.out" file
            =-1 to suppress writing of "target.out"
       MXNAT=dimensioning information (max number of atoms)
Output:
       A1(1-3)=unit vector (1,0,0) defining target axis 1 in Target Fr
       A2(1-3)=unit vector (0,1,0) defining target axis 2 in Target Fr
       NAT=number of atoms in target
       IXYZ(1-NAT,1-3)=x/d,y/d,z/d for atoms of target
       X0(1-3)=location/d in TF of dipole with IXYZ=0,0,0
               "upper" surface of the rectangular block is in x_TF=0 plane
               origin of the TF is set to be at the midpoint of the upper surface

B.T.Draine, Princeton Univ. Obs.
Fortren version records removed.

Copyright (C) 1993,1995,1996,1997,1998,2000,2007,2008, B.T. Draine and P.J. Flatau
Copyright (C) 2012, 2013, C++ version, Choliy V.
This code is covered by the GNU General Public License.
** */

Tarrec::Tarrec(TargetManager *man) : AbstractTarget(man)
{ 
	shortDescr = string("Tarrec");
	longDescr = string("Rectangular prism");
}

void Tarrec::Sizer(void)
{
	dx = manager->CashedDx();
	AllocateArrays((int)(shpar[0] / dx.data[0] + half_), (int)(shpar[1] / dx.data[1] + half_), (int)(shpar[2] / dx.data[2] + half_));
}

void Tarrec::Descriptor(void)
{
	sprintf(freeDescr, " %s; NX,NY,NZ=%4d%4d%4d", longDescr.c_str(), nx, ny, nz);
}

void Tarrec::Vector(void)
{
	VectorA();
}

void Tarrec::VectorX(void)
{
//
// Top dipole layer will be at x_TF=-d/2, and will have JX=NX therefore set X0(1)=-NX-0.5
// we want x_TF axis (i.e., y_TF=z_TF=0) to run through middle of target
// JY runs from 1 to NY : middle is at 0.5*(1+NY)
//                        X0(2)=-0.5*(1+NY)
//                        X0(3)=-0.5*(1+NZ)
	x0.Set(-(real)nx - half_, -half_ * (real)(ny + 1), -half_ * (real)(nz + 1));
}

void Tarrec::Allocator(void)
{
//
// Now populate lattice with homogeneous, isotropic composition.
	int curSize = nx * ny;
	ixyz.Dimension(curSize, 3);
//
	minJx = minJy = minJz = 1;
	maxJx = nx;
	maxJy = ny;
	maxJz = nz;
//
	nat0 = 0;
	for(int jx=0; jx<nx; ++jx)
	{
		for(int jy=0; jy<ny; ++jy)
		{
			for(int jz=0; jz<nz; ++jz)
			{
				if (nat0 >= curSize)
				{
					ixyz.Extend(nz);
					curSize = ixyz.GetSize(0);
				}
				ixyz.Fill3(nat0, jx+1, jy+1, jz+1);					// Positions are expressed starting from one, not zero 
				int index = GetLinearAddress(nat0);
				Composer(index, 0);
				iocc[index] = true;
				++nat0;
			}
		}
	}
	ixyz.Close(nat0);
}

void Tarrec::OutShpar(FILE *file)
{
	fprintf(file, " AX,AY,AZ=%8.4lf%8.4lf%8.4lf\n", shpar[0], shpar[1], shpar[2]); 
}

const char *TargetVerboseDescriptor_Rctglprsm(int num)
{
	static const char *descr[3] = 
	{
		"x-length", "y-length", "z-length"
	};
	if (num < 3)
		return descr[num];
	else
		return NULL;
}

REGISTER_TARGET(Rctglprsm,3,false,-1,0,"rectangular prism") 
void Target_Rctglprsm::SayHello(FILE *stream)
{
	fprintf(stream, "The target %s parameters are:\n", shortDescr.c_str());
	fprintf(stream, "x-length, y-length, z-length:\n");
	fprintf(stream, "%20.16lf %20.16lf %20.16lf\n", shpar[0], shpar[1], shpar[2]);
}

const char *TargetVerboseDescriptor_Rctglpbc(int num)
{
	static const char *descr[5] = 
	{
		"block dimensions X/d", "block dimensions Y/d", "block dimensions Z/d", "period in x/d", "period in y/d"
	};
	if (num < 5)
		return descr[num];
	else
		return NULL;
}

REGISTER_TARGET(Rctglpbc,5,false,3,0,"1 block (one TUC)")
void Target_Rctglpbc::SayHello(FILE *stream)
{
	fprintf(stream, "The target %s parameters are:\n", shortDescr.c_str());
	fprintf(stream, "block dimensions X/d, Y/d, Z/d:\n");
	fprintf(stream, "%20.16lf %20.16lf %20.16lf\n", shpar[0], shpar[1], shpar[2]);
	fprintf(stream, "next two parameters s.b. zeros, %20.16lf %20.16lf\n", shpar[3], shpar[4]);
}

void Target_Rctglpbc::PreparePyzd()
{
	pyd = shpar[3] / dx.data[1];
	pzd = shpar[4] / dx.data[2];
}
