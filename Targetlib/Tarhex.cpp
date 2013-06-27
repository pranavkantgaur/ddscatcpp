#include "StdAfx.h"

#include "Tarhex.h"
#include "TargetManager.h"

/* **
Routine to construct hexagonal prism from "atoms"

Input:
       A = prism length/d
       B = (vertex-vertex diameter of hexagon face)/d = 2*(one hexagon side)/d
       ORI = 1. for a1 in x direction, a2 in y direction
           = 2. for a1 in x direction, a2 in z direction
           = 3. for a1 in y direction, a2 in x direction
           = 4. for a1 in y direction, a2 in z direction
           = 5. for a1 in z direction, a2 in x direction
           = 6. for a1 in z direction, a2 in y direction (a1 is parallel to prism axis [normal to hexagonal face] a2 is normal to rectangular face
       DX(1-3)=(dx,dy,dz)/d where dx,dy,dz=lattice spacings in x,y,z directions, and d=(dx*dy*dz)**(1/3)=effective lattice spacing
       IOSHP=device number for "target.out" file
            =-1 to suppress printing of "target.out"
       MXNAT=dimensioning information (max number of atoms)

Output:
       A1(1-3)=unit vector (1,0,0) defining target axis 1 (prism axis)
       A2(1-3)=unit vector (0,1,0) defining target axis 2 (normal to one of the rectangular faces.
       X0(1-3)=location/d in TF of dipole with IXYZ=0 0 0
       NAT=number of atoms in target
       IXYZ(1-NAT,1-3)=x/d,y/d,z/d for atoms of target
       CDESCR=string describing target (up to 67 chars.)

Fortran history records removed.

Copyright (C) 1993,1995,1996,1997,1998,2000,2006,2007,2008 B.T. Draine and P.J. Flatau
Copyright (c) 2012, C++ version, Choliy V.
This code is covered by the GNU General Public License.
** */

Tarhex::Tarhex(TargetManager *man) : AbstractTarget(man) 
{ 
	shortDescr = string("Tarhex");
	longDescr = string("Hexagonal prism");
	iori = (int)shpar[2];
	acm = bcm = ccm = (real)0.;
}

void Tarhex::Sizer(void)
{
	a1.Clear();
	a2.Clear();
	dx = manager->CashedDx();
//
// A = prism length
// B = prism "diameter"
//       B/2 = length of one side of hexagon
//       B= distance between opposite vertices=max.diameter of hexagon
//       for hexagon, area=(3/8)*SQRT(3)*B**2
//       for hex prism, volume=A*area=(3/8)*SQRT(3)*A*B**2
// Now determine limits for testing x,y,z values
// Along axis, run from 1 to NA=INT(A+0.5)
// Perp. to axis, run from 1 to NB=INT(B+0.5)
	int na = (int)(shpar[0] + half_);
	int nb = (int)(shpar[1] + half_);
	int nc = (int)(half_ * sqrt(3.) * shpar[1] + half_);
//
// ACM="center" of figure in axial direction
// BCM="center" of figure in B direction
// CCM="center" of figure in C direction
	acm = (real)na / twox_ + half_;
	bcm = (real)nb / twox_ + half_;
	ccm = (real)nc / twox_ + half_;

	minJx = minJy = minJz = 1;
	switch(iori)
	{
	case 1:
		maxJx = na; maxJy = nc; maxJz = nb; a1.data[0] = a2.data[1] = onex_;
		break;

	case 2:
		maxJx = na; maxJz = nc; maxJy = nb;	a1.data[0] = a2.data[2] = onex_;
		break;

	case 3:
		maxJy = na; maxJx = nc; maxJz = nb;	a1.data[1] = a2.data[0] = onex_;
		break;
		
	case 4:
		maxJy = na; maxJz = nc; maxJx = nb;	a1.data[1] = a2.data[2] = onex_;
		break;
		
	case 5:
		maxJz = na; maxJx = nc; maxJy = nb;	a1.data[2] = a2.data[0] = onex_;
		break;
		
	case 6:
		maxJz = na; maxJy = nc; maxJx = nb;	a1.data[2] = a2.data[1] = onex_;
		break;

	default:
		Errmsg("Fatal", shortDescr.c_str(), " invalid value for SHPAR[2])");
		break;
	}
	AllocateArrays(maxJx - minJx + 1, maxJy - minJy + 1, maxJz - minJz + 1);
}

void Tarhex::Descriptor(void)
{
	sprintf(freeDescr, " Hexagonal prism of NAT=%7d dipoles", nat0);
}

void Tarhex::Vector(void)
{

}

void Tarhex::VectorX(void)
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

void Tarhex::Allocator(void)
{
//
// Current version of TARHEX is restricted to cubic lattices
	if ((dx.data[0] != onex_) || (dx.data[1] != onex_))
		Errmsg("Fatal", shortDescr.c_str(), " Tarhex does not support noncubic lattice");
//
	int curSize = nx * ny;
	ixyz.Dimension(curSize, 3);

	nat0 = 0;
	for(int jx=minJx; jx<=maxJx; ++jx)
	{
		for(int jy=minJy; jy<=maxJy; ++jy)
		{
			for(int jz=minJz; jz<=maxJz; ++jz)
			{
				bool in = false;
				switch(iori)
				{
				case 1:
					in = Testhex((real)jx, (real)jy, (real)jz);
					break;
               
				case 2:
					in = Testhex((real)jx, (real)jz, (real)jy);
					break;

				case 3:
					in = Testhex((real)jy, (real)jx, (real)jz);
					break;

				case 4:
					in = Testhex((real)jy, (real)jz, (real)jx);
					break;

				case 5:
					in = Testhex((real)jz, (real)jx, (real)jy);
					break;

				case 6:
					in = Testhex((real)jz, (real)jy, (real)jx);
					break;

				default:
					in = false;
					break;
				}
				if(in)
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

bool Tarhex::Testhex(real x, real y, real z)
{
// ACM,BCM,CCM=location of center in A,B,C directions
//    A direction is along axis
//    B direction is vertex to vertex of hexagon
//    C direction is face to face of hexagon
	real ax = x;
	real u = y - ccm;
	real v = z - bcm;
//
// Test along axis:
	if ((real)2. * Fabs(ax - acm) > shpar[0])
		return false;

	if (Fabs(u) > (real)0.4330127 * shpar[1]) 
		return false;
// 
// Now test whether closer to CM than line bounding edge 0.5773503=1/SQRT(3)
	real q = half_ * Fabs(u) + (real)0.8660254 * Fabs(v);
	if (q > (real)0.4330127 * shpar[1])
		return false;

	return true;
}

void Tarhex::OutShpar(FILE *file)
{
	fprintf(file, " A,B = %8.4lf%8.4lf\n", shpar[0], shpar[1]); 
}

const char *TargetVerboseDescriptor_Hexgonpbc(int num)
{
	return NULL;
}

REGISTER_TARGET(Hexgonpbc,5,false,3,0,"hexagonal slab (one TUC)")
void Target_Hexgonpbc::SayHello(FILE *stream)
{
	fprintf(stream, "The target %s parameters are:\n", shortDescr.c_str());
	fprintf(stream, "length (along symmetry axis):\n");
	fprintf(stream, "%20.16lf\n", manager->GetShpar(0));
	fprintf(stream, "2*length of 1 hexagon side/d = (max diameter)/d:\n");
	fprintf(stream, "%20.16lf\n", manager->GetShpar(1));
	fprintf(stream, "orientation of axes a1 (hex axis) and a2 is %d:\n", (int)shpar[2]);
	fprintf(stream, " 1 for a1=x,a2=y; 2 for a1=x,a2=z; 3 for a1=y,a2=x; 4 for a1=y,a2=z; 5 for a1=z,a2=x; 6 for a1=z,a2=y\n");	
	fprintf(stream, "next two parameters s.b. zeros, %20.16lf %20.16lf\n", shpar[3], shpar[4]);
}

const char *TargetVerboseDescriptor_HexPrism(int num)
{
	static const char *descr[3] = 
	{
		"length along symmetry axis", "length of 1 hehagon side/d", "orientation of axes a1 and a2 (1=xy, 2=xz, 3=yx, 4=yz, 5=zx, 6=zy)"
	};
	if (num < 3)
		return descr[num];
	else
		return NULL;
}

REGISTER_TARGET(HexPrism,3,false,-1,0,"Hexagonal prism")
void Target_HexPrism::SayHello(FILE *stream)
{
	fprintf(stream, "The target %s parameters are:\n", shortDescr.c_str());
	fprintf(stream, "length (along symmetry axis):\n");
	fprintf(stream, "%20.16lf\n", manager->GetShpar(0));
	fprintf(stream, "length of 1 hexagon side:\n");
	fprintf(stream, "%20.16lf\n", manager->GetShpar(1));
	fprintf(stream, "orientation of axes a1 (hex axis) and a2 is %d:\n", (int)shpar[2]);
	fprintf(stream, " 1 for a1=x,a2=y; 2 for a1=x,a2=z; 3 for a1=y,a2=x; 4 for a1=y,a2=z; 5 for a1=z,a2=x; 6 for a1=z,a2=y\n");
}
