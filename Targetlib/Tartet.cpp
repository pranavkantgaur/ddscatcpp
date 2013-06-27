#include "StdAfx.h"

#include "Tartet.h"
#include "TargetManager.h"

/* **
Routine to construct regular tetrahedral target array one of the faces is parallel to the y-z plane one of the edges of this face is parallel to x-y plane

Input:
       AX = length of one edge of tetrahedron
       DX(1-3)=(dx,dy,dz)/d where dx,dy,dz=lattice spacings in x,y,z directions, and d=(dx*dy*dz)**(1/3)=effective lattice spacing
       IOSHP = device number for "target.out" file
             = -1 to suppress printing of "target.out"
       MXNAT = dimensioning information

Returns:
       A1(1-3) = unit vector (1,0,0) (along one axis of tetrahedron)
       A2(1-3) = unit vector (0,1,0)
       CDESCR = string describing target (up to 67 chars.)
       NAT = number of atoms in target
       IXYZ(1-NAT,1-3)=x/d,y/d,z/d for atoms in target
       ICOMP(1-NAT,1-3)=1 (composition identifier)
       X0(1-3) = location/d in TF of site with IXYZ=0 0 0
Occupied array sites are those within tetrahedral surface
Size:
   S=length of each side of tetrahedron
   Volume = S**3/(6*sqrt(2))
Orientation:
   Center of mass is at origin.
   One face is parallel to yz plane.
   Projection onto yz plane has vertex on y axis.
   S=length of one side (in lattice units)
   Vertices are at
   A=( sqrt(3/8),         0,    0)*S
   B=(-sqrt(1/24), sqrt(1/3),    0)*S
   C=(-sqrt(1/24),-sqrt(1/12),-1/2)*S
   D=(-sqrt(1/24),-sqrt(1/12), 1/2)*S
Length in x direction = S*sqrt(2/3)
          y           = S*sqrt(3)/2
          z           = S
Angle(AOB)=arccos(-1/3)=109.4712 deg.
OA=OB=OC=OD=sqrt(3/8)*S
Occupied sites are assumed to be located at
(X,Y,Z)=(I+XOFF,J+YOFF,K+ZOFF)
where I,J,K are integers, and XOFF,YOFF,ZOFF are constants.
Program sets XOFF,YOFF,ZOFF depending on choice of parameter S.

Criterion for choosing XOFF:
   Base of tetrahedron is located at x = -sqrt(1/24)*S
   Let IMIN be value of I for this plane
   Choose XOFF so that IMIN+XOFF = -sqrt(1/24)*S + 0.5
   with -0.5 < XOFF < 0.5
Criterion for choosing YOFF:
   One edge of tetrahedron is located at y= -S/(2*sqrt(3))
   Let JMIN be value of J for this line
   Choose YOFF so that JMIN+YOFF = -S/sqrt(12) + 0.5
Criterion for choosing ZOFF:
   One edge of tetrahedron is parallel to z axis
   Choose ZOFF in order to have number of dipoles along
   this edge as close as possible to S
   e.g., if S=odd int, then take ZOFF=0 to place dipoles
         at integral locations
         if S=even int, then take ZOFF=1/2 to place dipoles
         at half-integral locations

Fortran history records removed.

Copyright (C) 1993,1995,1996,1997,1998,2000,2007,2008, B.T. Draine and P.J. Flatau
Copyright (c) 2012, 2013, C++ version, Choliy V.
This code is covered by the GNU General Public License.
** */

Tartet::Tartet(TargetManager *man) : AbstractTarget(man) 
{
	shortDescr = string("Tartet");
	longDescr = string("tetrahedral grain");
	xoff = yoff = zoff = (real)0.;
}

void Tartet::Descriptor(void)
{
	sprintf(freeDescr, " Tetrahedron of NAT=%7d dipoles", nat0);
}

void Tartet::Vector(void)
{
	VectorA();
}

void Tartet::VectorX(void)
{
//
// Initialize composition:
// Set X0 so that origin of TF is at centroid
	x0.Clear();
	for(int i=0; i<nat0; ++i)
	{
		x0 += Vect3<real>(ixyz.Value(i, 0), ixyz.Value(i, 1), ixyz.Value(i, 2));
	}
	x0 /= -(real)nat0;
}

void Tartet::Sizer(void)
{
	dx = manager->CashedDx();
//
// Current version of TARTET is restricted to cubic lattices
	real s = shpar[0];
//
// Set XOFF (and IMIN,IMAX):
	minJx = -(int)(s * Sqrt(onex_ / (real)24.));
	xoff = half_ - s * Sqrt(onex_ / (real)24.) - minJx;
	maxJx = minJx + (int)(s * Sqrt(twox_ / (real)3.) + half_) - 1;
//
// Set YOFF (and JMIN,JMAX):
	minJy = -(int)(s / Sqrt(12.));
	yoff = half_ - s / Sqrt(12.) - minJy;
	maxJy = minJy + (int)(s * Sqrt(0.75) + half_) - 1;
//
// Set ZOFF (and KMIN,KMAX): Determine whether S is closest to even or odd int. (Temporarily let KMIN be int which S is closest to)
	minJz = (int)(s + half_);
//
// If KMIN is even, then ZOFF=0.5
// If KMIN is odd, then ZOFF=0.
	zoff = zero_;
	if (minJz % 2 == 0) 
		zoff = half_;
	minJz = -(int)(half_ * s + zoff);
	maxJz = minJz + (int)(s + half_) - 1;
	AllocateArrays(maxJx - minJx + 1, maxJy - minJy + 1, maxJz - minJz + 1);
}

void Tartet::Allocator(void)
{
	if ((dx.data[0] != onex_) || (dx.data[1] != onex_))
		Errmsg("Fatal", "Tartet", " tartet does not support noncubic lattice");
//
// Determine list of occupied sites.
	int curSize = nx * ny;
	ixyz.Dimension(curSize, 3);

	nat0 = 0;
	for(int i=minJx; i<=maxJx; ++i)
	{
		real x = i + xoff;
// YMAX=largest value of Y which can occur for this X value
// YMIN=smallest value of Y which can occur for this X value
// ZMAX0=largest value of Z which can occur for this X value
		real ymax = (real)( shpar[0] * Sqrt(3./16.) - x / Sqrt(twox_));
		real ymin = (real)(-shpar[0] * Sqrt(3./64.) + x / Sqrt(8.));
		real zmax0 = (real)(3. * shpar[0] / 8. - x * Sqrt(3./8.));
		for(int j=minJy; j<=maxJy; ++j)
		{
			real y = j + yoff;
			if ((y >= ymin) && (y <= ymax))
			{
				real fy = (y - ymin) / (ymax - ymin);
				real zmax = (onex_ - fy) * zmax0;				// ! ZMAX=largest value of Z which can occur for this (X,Y)
				for(int k=minJz; k<=maxJz; ++k)
				{
					real z = k + zoff;
					if (Fabs(z) <= zmax)					// ! Site is occupied:
					{
						if (nat0 >= curSize)
						{
							ixyz.Extend(nz);
							curSize = ixyz.GetSize(0);
						}
						ixyz.Fill3(nat0, i, j, k);
						int index = GetLinearAddress(nat0);
						Composer(index, 0);
						iocc[index] = true;
						++nat0;
					}
				}
			}
		}
	}
	ixyz.Close(nat0);
}

void Tartet::OutShpar(FILE *file)
{
	fprintf(file, "S = %8.4lf\n", shpar[0]);
}

const char *TargetVerboseDescriptor_Tetrahdrn(int num)
{
	static const char *descr[1] = 
	{
		"side length"
	};
	if (num < 1)
		return descr[num];
	else
		return NULL;
}

REGISTER_TARGET(Tetrahdrn,1,false,-1,0,"Regular tetrahedron")
void Target_Tetrahdrn::SayHello(FILE *stream)
{
	fprintf(stream, "The target %s parameters are:\n", shortDescr.c_str());
	fprintf(stream, "tetrahedron side length, s.b. equal:\n");
	fprintf(stream, "%lf\n", shpar[0]);
}
