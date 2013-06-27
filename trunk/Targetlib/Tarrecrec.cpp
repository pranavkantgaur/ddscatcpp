#include "StdAfx.h"

#include "Tarrecrec.h"
#include "TargetManager.h"

/* **
Routine to generate target consisting of rectangular block of
       composition 1 resting on top of rectangular block of composition 2
       TF origin = center of top surface of block 1 i.e. 0.5d above the top layer of dipoles

Input:
        XV1=(x-length)/d  of block 1    (d=lattice spacing)
        YV1=(y-length)/d           1
        ZV1=(z-length)/d           1
        XV2=(x-length)/d  of block 2
        YV2=(y-length)/d  of block 2
        ZV2=(z-length)/d  of block 2

        DX(1-3)=(dx,dy,dz)/d where dx,dy,dz=lattice spacings in x,y,z directions, and d=(dx*dy*dz)**(1/3)=effective lattice spacing
        IOSHP = device number for "target.out" file
              = -1 to suppress printing of "target.out"
        MXNAT=dimensioning information (max number of atoms)

Output:
        A1(1-3)=unit vector (1,0,0) defining target axis 1 in Target Fr
        A2(1-3)=unit vector (0,1,0) defining target axis 2 in Target Fr
        NAT=number of atoms in target
        IXYZ(1-NAT,1-3)=x/d,y/d,z/d for atoms of target
        X0(3)=(x,y,z)_TF corresponding to lattice site IXYZ=(0,0,0)

History records of Fortran versions removed.

Copyright (C) 1993,1995,1996,1997,1998,2000,2007,2008,2012, B.T. Draine and P.J. Flatau
Copyright (c) 2012,2013 C++ versions, Choliy V.
This code is covered by the GNU General Public License.
** */

Tarrecrec::Tarrecrec(TargetManager *man) : AbstractTarget(man) 
{ 
	shortDescr = string("Tarrecrec");
	longDescr = string("two rectangular solids");
	jx1min = jx1max = jx2min = jx2max = 0;
	jy1min = jy1max = jy2min = jy2max = 0;
	jz1min = jz1max = jz2min = jz2max = 0;
}

void Tarrecrec::Sizer(void)
{
	dx = manager->CashedDx();
	int nx1 = nint_(shpar[0] / dx.data[0]);
	int ny1 = nint_(shpar[1] / dx.data[1]);
	int nz1 = nint_(shpar[2] / dx.data[2]);
	int nx2 = nint_((shpar[0] + shpar[3]) / dx.data[0]) - nx1;
	int ny2 = nint_(shpar[4] / dx.data[1]);
	int nz2 = nint_(shpar[5] / dx.data[2]);
//
// Bottom of upper box and top of lower box will be in x_TF=0 plane JX=0 corresponds to x_TF=dx/2
//     upper box:  JX runs from 0 to NX1-1 for material 1 example: NX1=1, JX runs from 0 to 1
//     lower box:  JX runs from -NX2 to -1 for material 2 example: NX2=2, JX runs from -2 to -1
//
// If NY1 is even, JY=0 corresponds to y_TF = -dy/2
//     upper box:  JY runs from 1-(NY1/2) to (NY1/2 example: NY1=2, JY runs from 0 to +1
//     lower box:  JY runs from 1-int(NY2/2) to NY2-int(NY2/2) example1: NY2=2, JY runs from 0 to 1 example2: NY2=3, JY runs from 0 to 2
//
//    NY1 is  odd, JY=0 corresponds to y_TF=0
//     upper box:  JY runs from -(NY1-1)/2 to (NY1-1)/2 example: NY1=3, JY runs from -1 to +1
//     lower box:  JY runs from -int[(NY2-1)/2] to NY2-1-int[(NY2-1)/2]
//
// If NZ1 is even, JZ=0 corresponds to z_TF = -dz/2
//     upper box:  JZ runs from 1-(NZ1/2) to (NZ1/2)
//     lower box:  JZ runs from 1-int(NZ2/2) to NZ2-int(NZ2/2)
//
//    NZ1 is odd,  JZ=0 corresponds to z_TF=0
//     upper box:  JZ runs from -(NZ1-1)/2 to (NZ1-1)/2
//     lower box:  JZ runs from -int[(NZ2-1)/2] to NZ2-1-int[(NZ2-1)/2]

	x0.data[0] = half_;

	jx1min = 0;
	jx1max = nx1 - 1;
	jx2min = -nx2;
	jx2max = -1;
	if (ny1 % 2 == 0)											// ny1 is even
	{
		x0.data[1] = -half_;
		jy1min = 1 - ny1/2;
		jy1max = ny1/2;
		jy2min = 1 - ny2/2;
		jy2max = ny2 - ny2/2;
	}
	else														// ny1 is odd
	{
		x0.data[1] = zero_;
		jy1min = -(ny1-1)/2;
		jy1max =  (ny1-1)/2;
		jy2min = -(ny2-1)/2;
		jy2max =  ny2 - 1 - (ny2-1)/2;
	}
//
	if (nz1 % 2 == 0)											// nz1 is even
	{
		x0.data[2] = -half_;
		jz1min = 1 - nz1/2;
		jz1max = nz1/2;
		jz2min = 1 - nz2/2;
		jz2max = nz2 - nz2/2;
	}
	else														// nz2 is odd
	{
         x0.data[2] = zero_;
         jz1min = -(nz1-1) / 2;
         jz1max =  (nz1-1) / 2;
         jz2min = -(nz2-1) / 2;
         jz2max = nz2 - 1 - (nz2-1)/2;
	}
	AllocateArrays(max(jx2max, jx1max) - min(jx2min, jx1min) + 1, max(jy2max, jy1max) - min(jy2min, jy1min) + 1, max(jz2max, jz1max) - min(jz2min, jz1min) + 1); 
	minJx = min_(jx2min, jx1min);
	maxJx = max_(jx2max, jx1max);
	minJy = min_(jy2min, jy1min);
	maxJy = max_(jy2max, jy1max);
	minJz = min_(jz2min, jz1min);
	maxJz = max_(jz2max, jz1max);
}

void Tarrecrec::Descriptor(void)
{
	sprintf(freeDescr, "Two rect.slabs NX1,NY1,NZ1=%4d%4d%4d NX2,NYZ,NZ2=%4d%4d%4d", minJx, minJy, minJz, maxJx, maxJy, maxJz);
}

void Tarrecrec::Vector(void)
{
	VectorA();
}

void Tarrecrec::VectorX(void)
{

}

void Tarrecrec::Allocator(void)
{
//
// Now populate upper rectangular slab:
	int curSize = nx * ny;
	ixyz.Dimension(curSize, 3);

	int jx, jy, jz;
	nat0 = 0;
	for(jx=jx1min; jx<=jx1max; ++jx)
	{
		for(jy=jy1min; jy<=jy1max; ++jy)
		{
			for(jz=jz1min; jz<=jz1max; ++jz)
			{
				if (nat0 >= curSize)
				{
					ixyz.Extend(nz);
					curSize = ixyz.GetSize(0);
				}
				ixyz.Fill3(nat0, jx, jy, jz);			// ! homogeneous, isotropic composition:
				int index = GetLinearAddress(nat0);
				Composer(index, 0);
				iocc[index] = true;
				++nat0;
			}
		}
	}
//
// Lower rectangular slab
	for(jx=jx2min; jx<=jx2max; ++jx)
	{
		for(jy=jy2min; jy<=jy2max; ++jy)
		{
			for(jz=jz2min; jz<=jz2max; ++jz)
			{
				if (nat0 >= curSize)
				{
					ixyz.Extend(nz);
					curSize = ixyz.GetSize(0);
				}
				ixyz.Fill3(nat0, jx, jy, jz);			// ! homogeneous, isotropic composition:
				int index = GetLinearAddress(nat0);
				Composer(index, 1);
				iocc[index] = true;
				++nat0;
			}
		}
	}
	ixyz.Close(nat0);
}

void Tarrecrec::OutShpar(FILE *file)
{
	fprintf(file, " NX1,NY1,NZ1 = %4d%4d%4d NX2,NY2,NZ2 = %4d%4d%4d\n", (int)shpar[0], (int)shpar[1], (int)shpar[2], (int)shpar[3], (int)shpar[4], (int)shpar[5]);
}

const char *TargetVerboseDescriptor_Rctgrctg(int num)
{
	static const char *descr[6] = 
	{
		"x length/d for first block", "y length/d for first block", "z length/d for first block", 
		"x length/d for second block", "y length/d for second block", "z length/d for second block"
	};
	if (num < 6)
		return descr[num];
	else
		return NULL;
}

REGISTER_TARGET(Rctgrctg,6,false,-1,0,"Two rectangular solids")
void Target_Rctgrctg::SayHello(FILE *stream)
{
	fprintf(stream, "The target %s parameters are:\n", shortDescr.c_str());
	fprintf(stream, "x,y,z length/d for first block:\n");
	fprintf(stream, "%lf %lf %lf\n", shpar[0], shpar[1], shpar[2]);
	fprintf(stream, "x,y,z length/d for second block:\n");
	fprintf(stream, "%lf %lf %lf\n", shpar[3], shpar[4], shpar[5]);
}

const char *TargetVerboseDescriptor_Recrecpbc(int num)
{
	static const char *descr[8] = 
	{
		"x length/d for first block", "y length/d for first block", "z length/d for first block", 
		"x length/d for second block", "y length/d for second block", "z length/d for second block",
		"this parameter s.b. zero", "this parameter s.b. zero too"
	};
	if (num < 8)
		return descr[num];
	else
		return NULL;
}

REGISTER_TARGET(Recrecpbc,8,false,6,0,"2 stacked blocks (one TUC)")
void Target_Recrecpbc::SayHello(FILE *stream)
{
	fprintf(stream, "The target %s parameters are:\n", shortDescr.c_str());
	fprintf(stream, "x,y,z length/d for first block:\n");
	fprintf(stream, "%lf %lf %lf\n", shpar[0], shpar[1], shpar[2]);
	fprintf(stream, "x,y,z length/d for second block:\n");
	fprintf(stream, "%lf %lf %lf\n", shpar[3], shpar[4], shpar[5]);
	fprintf(stream, "next two parameters s.b. zeros, %lf %lf\n", shpar[6], shpar[7]);
}
