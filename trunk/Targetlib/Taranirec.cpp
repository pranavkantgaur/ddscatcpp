#include "StdAfx.h"

#include "Taranirec.h"
#include "TargetManager.h"

/* **
Subroutine TARANIREC
Purpose: to construct rectangular target with anisotropic dielectric tensor.
         Dielectric tensor is assumed to be diagonalized in target frame, with
         material 1 for x-direction
         material 2 for y-direction
         material 3 for z-direction
         Note: Must specify NCOMP=3 in ddscat.par, and provide
               file names for three materials
Input:
       XV=(x-length)/d    (d=lattice spacing)
       YV=(y-length)/d
       ZV=(z-length)/d
       DX(1-3)=(dx,dy,dz)/d where dx,dy,dz=lattice spacings in x,y,z
               directions, and d=(dx*dy*dz)**(1/3)=effective lattice spacing
       IOSHP=device number for "target.out" file
            =-1 to suppress printing of "target.out"
       MXNAT=dimensioning information (max number of atoms)
Output:
       A1(1-3)=unit vector (1,0,0) defining target axis 1 in Target Frame
       A2(1-3)=unit vector (0,1,0) defining target axis 2 in Target Frame
       NAT=number of atoms in target
       IXYZ(1-NAT,1-3)=x/d,y/d,z/d for atoms of target
       X0(1-3)=(location/d) in TF of dipole with IXYZ=(0,0,0).  This
               will henceforth be treated as the origin of physical
               coordinates in the TF.
               origin is assumed to be located at midpoint of upper surface
               normal to the A1 vector.
               Upper surface is presumed to be 0.5*d above uppermost
               layer of dipoles.

B.T.Draine, Princeton Univ. Obs.
Fortran history records removed.

Copyright (C) 1993,1995,1996,1997,1998,2000,2004,2007,2008,2010
              B.T. Draine and P.J. Flatau
Copyright (C) 2012, 2013, C++ version, Choliy V.
This code is covered by the GNU General Public License.
** */

Taranirec::Taranirec(TargetManager *man) : AbstractTarget(man) 
{ 
	shortDescr = string("Taranirec");
	longDescr = string("Anisotropic rectangular prism");
}

void Taranirec::Sizer(void)
{
	dx = manager->CashedDx();
	AllocateArrays((int)(shpar[0] / dx.data[0] + half_), (int)(shpar[1] / dx.data[1] + half_), (int)(shpar[2] / dx.data[2] + half_));
}

void Taranirec::Descriptor(void)
{
	sprintf(freeDescr, " %s; NX,NY,NZ=%4d%4d%4d", longDescr.c_str(), nx, ny, nz);
}

void Taranirec::Vector(void)
{
	VectorA();
}

void Taranirec::VectorX(void)
{
//
// Set TF origin of coordinates to be located at midpoint of upper surface normal to A1.  Upper surface is presumed to be located 0.5*d above uppermost dipole layer
	x0.Set(-nx-half_, -(real)(1 + ny) / 2, -(real)(1 + nz) / 2);
}

void Taranirec::Allocator(void)
{
//
// Now populate lattice:
	int curSize = nx * ny;
	ixyz.Dimension(curSize, 3);
//
	minJx = minJy = minJz = 1;
	maxJx = nx;
	maxJy = ny;
	maxJz = nz;
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
				ixyz.Fill3(nat0, jx+1, jy+1, jz+1);
				int index = GetLinearAddress(nat0);
				icomp.Fill3(index, (short)1, (short)2, (short)3);			// Homogeneous, anisotropic composition:
				iocc[index] = true;
				++nat0;
			}
		}
	}
	ixyz.Close(nat0);
}

void Taranirec::OutShpar(FILE *file)
{
	fprintf(file, " AX,AY,AZ=%8.4lf%8.4lf%8.4lf\n", shpar[0], shpar[1], shpar[2]); 
}

const char *TargetVerboseDescriptor_Anirctngl(int num)
{
	static const char *descr[3] = 
	{
		"x-length/d", "y-length/d", "z-length/d"
	};
	if (num < 3)
		return descr[num];
	else
		return NULL;
}

REGISTER_TARGET(Anirctngl,3,false,-1,0,"Anisotropic rectangular solid")
void Target_Anirctngl::SayHello(FILE *stream)
{
	fprintf(stream, "The target %s parameters are:\n", shortDescr.c_str());
	fprintf(stream, "length/d in x, y, z directions are:\n");
	fprintf(stream, "%20.16lf %20.16lf %20.16lf\n", shpar[0], shpar[1], shpar[2]);
}

void Target_Anirctngl::PrepareIaniso(void)
{
	ianiso = TargetIsAnisotropic;
}
