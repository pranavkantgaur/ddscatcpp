#include "StdAfx.h"

#include "Tarslabhole.h"
#include "TargetManager.h"

/* **
Routine to construct rectangular slab with dimensions A*d, B*d, C*d
in x_tf,y_tf,z_tf directions with circular hole of radius R*d
hole axis is in x_tf direction, through center of slab

x0(1-3) = point at center of hole, at top surface of slab

Input:
       A=(x-length)/d     (d=lattice spacing)
       BA=B/A  
       CA=C/A
       RA=R/A

       DX(1-3)=(dx,dy,dz)/d where dx,dy,dz=lattice spacings in x,y,z
               directions, and d=(dx*dy*dz)**(1/3)=effective lattice spacing
       IOSHP = device number for "target.out" file
             = -1 to suppress printing of "target.out"
       MXNAT=dimensioning information (max number of atoms)
Output:
       A1(1-3)=unit vector (1,0,0) defining target axis 1 in Target Frame
       A2(1-3)=unit vector (0,1,0) defining target axis 2 in Target Frame
       NAT=number of atoms in target
       IXYZ(1-NAT,1-3)=location indices for atoms of target
       ICOMP(1-NAT,1-3)=composition
       X0(3)=location/d in TF of lattice site with IXYZ=(0,0,0)
             X(1-3)=[IXYZ(J,1-3)+X0(1-3)]*d

B.T.Draine, Princeton Univ. Obs.
Fortren version records removed.

Copyright (C) 2008,2010, B.T. Draine and P.J. Flatau
Copyright (C) 2012, 2013, C++ version, Choliy V.
This code is covered by the GNU General Public License.
** */

Tarslabhole::Tarslabhole(TargetManager *man) : AbstractTarget(man) 
{ 
	shortDescr = string("Tarslabhole");
	longDescr = string("rectangular block with circular hole");
}

void Tarslabhole::Sizer(void)
{
	dx = manager->CashedDx();
	AllocateArrays(nint_(shpar[0] / dx.data[0]), nint_(shpar[1] * shpar[0] / dx.data[1]), nint_(shpar[2] * shpar[0] / dx.data[2]));
}

void Tarslabhole::Descriptor(void)
{
	strcpy(freeDescr, "SLAB_HOLE: rect. block with cylindrical hole");
}

void Tarslabhole::Vector(void)
{
	VectorA();
}

void Tarslabhole::VectorX(void)
{
	x0.Set(half_ * dx.data[0], -half_ * (real)(1 + ny), -half_ * (real)(1 + nz));
}

void Tarslabhole::Allocator(void)
{
	if (twox_ * shpar[3] > Sqrt(shpar[1] * shpar[1] + shpar[2] * shpar[2])) 
		Errmsg("Fatal", "Tarslabhole", "hole diam > block diag");
//
// Now populate lattice:
// Lowest layer will have JX=-NX
//
// Specify X0(3) = (x,y,z)/d in TF corresponding to IXYZ=(0,0,0)
// Top surface of slab (0.5d above topmost dipole layer) is assumed to have x=0
// Dipoles with JX=-1 are assumed to be at x=-dx/2 (top layer of target).
// Therefore X0(1)=0.5*DX(1)
// Set y=0 and z=0 to be at middle of slab
// JY runs from 1 to NY -> X0(2) = -(1+NY)/2
// JZ runs from 1 to NZ -> X0(3) = -(1+NZ)/2

	VectorX();

	real r2 = shpar[3] * shpar[0];
	r2 = r2 * r2;
	nat0 = 0;

	int curSize = nx * ny;
	ixyz.Dimension(curSize, 3);

	minJx = -nx;
	maxJx = -1;
	minJy = 1;
	maxJy = ny;
	minJz = 1;
	maxJz = nz;
	for(int jx=-1; jx>=-nx; --jx)
	{
		for(int jy=0; jy<ny; ++jy)
		{
			real rady2 = (jy + 1 + x0.data[1]) / dx.data[1];
			rady2 = rady2 * rady2;
			for(int jz=0; jz<nz; ++jz)
			{
				real rad2 = (jz + 1 + x0.data[2]) / dx.data[2];
				rad2 = rad2 * rad2;
				rad2 += rady2;
				if (rad2 > r2)
				{
					if (nat0 >= curSize)
					{
						ixyz.Extend(nz);
						curSize = ixyz.GetSize(0);
					}
					ixyz.Fill3(nat0, jx, jy+1, jz+1);
					int index = GetLinearAddress(nat0);
					Composer(index, 0);
					iocc[index] = true;
					++nat0;
				}
			}
		}
	}
}

void Tarslabhole::OutShpar(FILE *file)
{
	fprintf(file, "A/d,B/d/C/d,R/d = %8.4lf %8.4lf %8.4lf %8.4lf\n", shpar[0], shpar[1], shpar[2], shpar[3]);
}

const char *TargetVerboseDescriptor_Slabhole(int num)
{
	static const char *descr[4] = 
	{
		"a/d for block", "aspect ratio b/a", "aspect ratios c/a", "r/a for cylindrical hole"
	};
	if (num < 4)
		return descr[num];
	else
		return NULL;
}

REGISTER_TARGET(Slabhole,4,false,-1,0,"Rect. block with cylindrical hole")
void Target_Slabhole::SayHello(FILE *stream)
{
	fprintf(stream, "Slab dimensions = a  b  c ; cylindrical hole radius = r\n");
	fprintf(stream, "The target %s parameters are:\n", shortDescr.c_str());
	fprintf(stream, "a/d for block, s.b. > 0.5:\n");
	fprintf(stream, "%lf\n", shpar[0]);
	fprintf(stream, "aspect ratios b/a  c/a\n");
	fprintf(stream, "%lf %lf\n", shpar[1], shpar[2]);
	fprintf(stream, "r/a for cylindrical hole\n");
	fprintf(stream, "%lf\n", shpar[3]);
	if (manager->GetShpar(3) >= half_ * Sqrt(shpar[1] * shpar[1] + shpar[2] * shpar[2]))
	{
		fprintf(stream, "hole is larger than y-z dimensions of block\n");
	}
}

const char *TargetVerboseDescriptor_Slbholpbc(int num)
{

	return NULL;
}

REGISTER_TARGET(Slbholpbc,6,false,4,0,"block with cylindrical hole (one TUC)")
void Target_Slbholpbc::SayHello(FILE *stream)
{
	fprintf(stream, "TUC dimensions = a  b  c ; cylindrical hole radius = r\n");
	fprintf(stream, "The target %s parameters are:\n", shortDescr.c_str());
	fprintf(stream, "a/d for TUC, s.b. >= 0.5:\n");
	fprintf(stream, "%20.16lf\n", shpar[0]);
	fprintf(stream, "aspect ratios b/a  c/a, s.b. >= 0.5:\n");
	fprintf(stream, "%20.16lf %20.16lf\n", shpar[1], shpar[2]);
	fprintf(stream, "r/a for cylindrical hole\n");
	fprintf(stream, "%20.16lf\n", shpar[3]);
	if (manager->GetShpar(3) < half_ * Sqrt(shpar[1] * shpar[1] + shpar[2] * shpar[2]))
	{
		fprintf(stream, "next two parameters s.b. zeros, %20.16lf %20.16lf\n", shpar[4], shpar[5]);
	}
	else
	{
		fprintf(stream, "hole is larger than y-z dimensions of block\n");
	}
}
