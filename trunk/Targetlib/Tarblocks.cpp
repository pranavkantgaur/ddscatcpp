#include "StdAfx.h"

#include "Tarblocks.h"
#include "TargetManager.h"

/* **
Subroutine to construct a target out of "cubic" blocks.
Blocks are only cubic for case of a cubic lattice DX=(1.,1.,1.)
In case of a noncubic lattice, "blocks" are rectangular, with an equal
number (=BLKSIZ) of dipole spacings in each direction.

Given:
DX(1-3)     =(dx/d,dy/d,dz/d) where dx,dy,dz = x,y,z lattice spacing,
             and d = (dx*dy*dz)**(1/3) = effective lattice spacing
NBLOCKS     =number of cubic blocks in target
BLKSIZ      =(integer) ratio of block width/dipole spacing
             (there will be BLKSIZ**3 dipoles per block)
XYZB(1-3,J) =x,y,z coordinates of block J, in units of block width
IPRINAX     =0. to return A1=(1,0,0), A2=(0,1,0)
            =1. to return A1,A2=principal axes with largest, 2nd large
             moment of inertia
IOSHP       =device number for output file "target.out"
            =-1 to suppress generation of file "target.out"
MXNAT       =limit on largest allowed value of NAT=number of dipoles

Returns:
NAT         =number of dipoles in target
IXYZ(J,1-3) =x,y,z location of dipole J on lattice
ICOMP(J,1-3)=composition for dipole J; x,y,z directions
A1(1-3)     =principal axis A1 in target frame (normalized)
A2(1-3)     =principal axis A2 in target frame (normalized)
X0(1-3)     =location in TF of dipole with IXYZ=0 0 0
             TF origin is taken to be at centroid of target
CDESCR      =string describing target

where A1 = principal axis with largest moment of inertia
      A2 = principal axis with second-largest moment of inertia

B.T. Draine, Princeton Univ. Observatory, 95.12.11
Fortran history records removed.

Copyright (C) 1996,1997,1998,2003,2004,2007,2008, B.T. Draine and P.J. Flatau
Copyright (C) 2012, 2013, C++ version, Choliy V.
This code is covered by the GNU General Public License.
** */

Tarblocks::Tarblocks(TargetManager *man) : AbstractTarget(man) 
{
	shortDescr = string("Tarblocks");
	longDescr = string("Cubic building blocks");
	blksiz = nblocks = iprinax = 0;
}

void Tarblocks::Sizer(void)
{
	real miX, maX, miY, maY, miZ, maZ;
	Loader(miX, maX, miY, maY, miZ, maZ);
	dx = manager->CashedDx();
	AllocateArrays((maX - miX + 1) * blksiz, (maY - miY + 1) * blksiz, (maZ - miZ + 1) * blksiz);
	int joff = 1 - (blksiz + 1) / 2;
	minJx = miX * blksiz + joff;
	minJy = miY * blksiz + joff;
	minJz = miZ * blksiz + joff;
	maxJx = maX * blksiz - joff + 1;
	maxJy = maY * blksiz - joff + 1;
	maxJz = maZ * blksiz - joff + 1;
}

void Tarblocks::Descriptor(void)
{
	sprintf(freeDescr, "Target containing %d dipoles arranged in %d cubic blocks", nat0, nblocks);
}

void Tarblocks::Vector(void)
{
	if (iprinax != 0)
		eigval = Prinaxis();				// Call PRINAXIS to compute principal axes (A1, A2 changed)
	else
		VectorA();							// Do not compute principal axes: simply set, a1=(1,0,0) and a2=(0,1,0) in target frame.
}

void Tarblocks::VectorX(void)
{
	x0.Clear();
	for(int jd=0; jd<nat0; ++jd)
	{
		x0 += Vect3<real>(ixyz.Value(jd, 0), ixyz.Value(jd, 1), ixyz.Value(jd, 2));
	}
	x0 /= -(real)nat0;
}

void Tarblocks::Loader(real &miX, real &maX, real &miY, real &maY, real &miZ, real &maZ)
{
	FILE *idvshpFile = fopen("blocks.par", "r");
	if (!idvshpFile)
		return;

	char Buffer[256], Format[256];
	sprintf(Format, "%s%s%s", realFormat, realFormat, realFormat);
    fgets(freeDescr, 68, idvshpFile);
    fgets(Buffer, 255, idvshpFile);
	sscanf(Buffer, "%d", &iprinax);
	fgets(Buffer, 255, idvshpFile);
	sscanf(Buffer, "%d", &nblocks);
	fgets(Buffer, 255, idvshpFile);
	sscanf(Buffer, "%d", &blksiz);
	miX = maX = miY = maY = miZ = maZ = zero_;
	for(int i=0; i<nblocks; ++i)
	{
		real a, b, c;
		fgets(Buffer, 255, idvshpFile);
		sscanf(Buffer, Format, &a, &b, &c);
		xyzb.push_back(new Vect3<real>(a, b, c));
		miX = min_(miX, a);		maX = max_(maX, a);
		miY = min_(miY, b);		maY = max_(maY, b);
		miZ = min_(miZ, c);		maZ = max_(maZ, c);
	}
	fclose(idvshpFile);
}

void Tarblocks::Allocator(void)
{
//
// For the moment, let us assume that blocks do not overlap.
	int curSize = nx * ny;
	ixyz.Dimension(curSize, 3);

	int joff = 1 - (blksiz + 1) / 2;
	nat0 = 0;
	for(int jb=0; jb<nblocks; ++jb)
	{
		int joffx = joff + (int)(blksiz * xyzb.at(jb)->data[0]);
		int joffy = joff + (int)(blksiz * xyzb.at(jb)->data[1]);
		int joffz = joff + (int)(blksiz * xyzb.at(jb)->data[2]);
		for(int jx=0; jx<blksiz; ++jx)
		{
			int jjx = jx + joffx;
			for(int jy=0; jy<blksiz; ++jy)
			{
				int jjy = jy + joffy;
				for(int jz=0; jz<blksiz; ++jz)
				{
					int jjz = jz + joffz;
					if (nat0 >= curSize)
					{
						ixyz.Extend(nz);
						curSize = ixyz.GetSize(0);
					}
                    ixyz.Fill3(nat0, jjx, jjy, jjz);
					int index = GetLinearAddress(nat0);
					Composer(index, 0);									// Homogeneous target:
					iocc[index] = true;
					++nat0;
				}
			}
		}
	}
	ixyz.Close(nat0);
}

const char *TargetVerboseDescriptor_Dw1996tar(int num)
{

	return NULL;
}

REGISTER_TARGET(Dw1996tar,1,false,-1,0,"13 cube target used by Draine & Weingartner (1996)")
void Target_Dw1996tar::SayHello(FILE *stream)
{
	fprintf(stream, "The target %s parameters are:\n", shortDescr.c_str());
	fprintf(stream, "Enter block size (lattice units):\n");
	fprintf(stream, "%20.16lf\n", shpar[0]);
}

const char *TargetVerboseDescriptor_Mltblocks(int num)
{

	return NULL;
}

REGISTER_TARGET(Mltblocks,0,false,-1,0,"construct from cubic blocks input from file <blocks.par>")
void Target_Mltblocks::SayHello(FILE *stream)
{
	fprintf(stream, "MLTBLOCKS : target defined by file blocks.par\n");
	fprintf(stream, "The target %s has no parameters.\n", shortDescr.c_str());
}
