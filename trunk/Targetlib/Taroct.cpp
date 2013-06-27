#include "StdAfx.h"

#include "Taroct.h"
#include "TargetManager.h"

/* **
Routine to construct octagonal prism from "atoms"

Input:
       A = prism length/d
       B = (vertex-vertex diameter of hexagon face)/d
         = 2*(one hexagon side)/d
       DX(1-3)=(dx,dy,dz)/d where dx,dy,dz=lattice spacings in x,y,z
               directions, and d=(dx*dy*dz)**(1/3)=effective lattice
               spacing
       IOSHP=device number for "target.out" file
            =-1 to suppress printing of "target.out"
       MXNAT=dimensioning information (max number of atoms)
Output:
       A1(1-3)=unit vector (1,0,0) defining target axis 1 (prism axis)
       A2(1-3)=unit vector (0,1,0) defining target axis 2 (normal to
               one of the rectangular faces.
       NAT=number of atoms in target
       IXYZ(1-NAT,1-3)=x/d,y/d,z/d for atoms of target
       CDESCR=string describing target (up to 67 chars.)

Fortran history records removed.

Copyright (C) 1993,1995,1996,1997,1998,2000 B.T. Draine and P.J. Flatau
Copyright (c) 2012, C++ version, Choliy V.
This code is covered by the GNU General Public License.
** */

Taroct::Taroct(TargetManager *man) : AbstractTarget(man) 
{ 
	shortDescr = string("Taroct");
	longDescr = string("Octagonal prism");
}

void Taroct::Sizer(void)
{
	dx = manager->CashedDx();
//
// A = prism length
// B = prism "diameter"
//        for hex prism, volume=A*area=(3/8)*SQRT(3)*A*B**2
// Now determine limits for testing x,y,z values
// Along axis, run from 1 to NA=INT(A+0.5)
// Perp. to axis, run from 1 to NB=INT(B+0.5)
	int jx = (int)(shpar[0] + half_);
	xoff = (jx % 2) ? zero_ : half_;
	int jy = (int)(shpar[1] + half_);
	yoff = (jy % 2) ? zero_ : half_;
	zoff = yoff;

	minJx = -(int)(half_ * shpar[0] + half_);
	maxJx =  (int)(half_ * shpar[0] - (real)0.25);
	minJy = -(int)(half_ * shpar[1] + half_);
	maxJy =  (int)(half_ * shpar[1] - (real)0.25);
	minJz = minJy;
	maxJz = maxJy;

	AllocateArrays(maxJx - minJx + 1, maxJy - minJy + 1, maxJz - minJz + 1);
}

void Taroct::Descriptor(void)
{
	sprintf(freeDescr, " Octagonal prism of NAT=%7d dipoles", nat0);
}

void Taroct::Vector(void)
{
	VectorA();
}

void Taroct::VectorX(void)
{
//
// Set X0 = location in TF of IXYZ=0 0 0
// We assume that origin of TF is located at centroid of hex prism
// Set composition
	x0.Clear();
	for(int i=0; i<nat0; ++i)
	{
		x0 += Vect3<real>((real)ixyz.Value(i, 0), (real)ixyz.Value(i, 1), (real)ixyz.Value(i, 2));
	}
	x0 /= -(real)nat0;
}

void Taroct::Allocator(void)
{
//
// Current version of TARHEX is restricted to cubic lattices
	if ((dx.data[0] != onex_) || (dx.data[1] != onex_))
		Errmsg("Fatal", shortDescr.c_str(), " Taroct does not support noncubic lattice");
//
	const real limx = onex_ / Sqrt(twox_);
	int curSize = nx * ny;
	ixyz.Dimension(curSize, 3);

	nat0 = 0;
	for(int jx=minJx; jx<=maxJx; ++jx)
	{
		for(int jy=minJy; jy<=maxJy; ++jy)
		{
			real y = Fabs(jy + yoff);
			for(int jz=minJz; jz<=maxJz; ++jz)
			{
				real z = Fabs(jz + zoff);
				if ((y + z) / shpar[1] <= limx)
				{
					if (nat0 >= curSize)
					{
						ixyz.Extend(nz);
						curSize = ixyz.GetSize(0);
					}
					ixyz.Fill3(nat0, jx, jy, jz);
					int index = GetLinearAddress(nat0);
					Composer(index, 0);
					iocc[index] = true;
					++nat0;
				}
			}
		}
	}
	ixyz.Close(nat0);

/* **
//
// Now compute some geometric properties of pseudo-hexagonal prism
//   NLAY = number of layers = effective length of prism/d
//   NFAC = number of atoms in one hexagonal face
	int nlay = (int)shpar[0];
	int nfac = nat0 / nlay;
//
// BEFF = effective length of hexagonal size/d computed for hexagon of area NFAC*d**2 =(3/2)*SQRT(3)*BEFF**2
// Note: for hexagon, vertex-vertex diameter = 2*BEFF
	real beff = (real)0.6204032 * Sqrt((real)nfac);
//
// ASPR = effective aspect ratio of target = length/(2*BEFF)
	real aspr = half_ * (real)nlay / beff;
//
// Here print any additional target info which is desired by using subroutine WRIMSG
	sprintf(cmsgnm, " NAT=%7d NFAC=%4d NLAY=%4d aspect ratio=%7.4lf", nat0, nfac, nlay, aspr);
	Wrimsg(shortDescr, cmsgnm);
** */
}

void Taroct::OutShpar(FILE *file)
{
	fprintf(file, " A,B = %8.4lf%8.4lf\n", shpar[0], shpar[1]); 
}

const char *TargetVerboseDescriptor_OctPrism(int num)
{
	static const char *descr[2] = 
	{
		"length along symmetry axis/d", "distance between octagohal faces/d",
	};
	if (num < 2)
		return descr[num];
	else
		return NULL;
}

REGISTER_TARGET(OctPrism,2,false,-1,0,"Octagonal prism")
void Target_OctPrism::SayHello(FILE *stream)
{
	fprintf(stream, "The target %s parameters are:\n", shortDescr.c_str());
	fprintf(stream, "length (along symmetry axis):\n");
	fprintf(stream, "%20.16lf\n", manager->GetShpar(0));
	fprintf(stream, "distance between octagohal faces:\n");
	fprintf(stream, "%20.16lf\n", manager->GetShpar(1));
}
