#include "StdAfx.h"

#include "Tarrctell.h"
#include "TargetManager.h"

/* **
Routine to construct target consisting of ellipsoid of dimensions
AX*d, AY*d, AZ*d
embedded in a concentric rectangular volume of dimensions
BX*d, BY*d, BZ*d

Input:
        AX=(x-length)/d of ellipsoid  (d=lattice spacing) of material 1
        AY=(y-length)/d "   "     "
        AZ=(z-length)/d "   "     "
        BX=(x-length)/d of larger rectangular volume of material 2
        BY=(y-length)/d "   "     "
        BZ=(z-length)/d "   "     "
        DX(1-3)=lattice spacing (dx/d,dy/d,dz/d), d=(dx*dy*dz)**(1/3)
        IOSHP=device number for "target.out" file =-1 to suppress printing of "target.out"
        MXNAT=dimensioning information (max number of atoms)

Output:
        A1(1-3)=(1,0,0)=unit vector defining target axis 1 in Target Frame
        A2(1-3)=(0,1,0)=unit vector defining target axis 2 in Target Frame
        NAT=number of atoms in target
        IXYZ(1-NAT,1-3)=x/d,y/d,z/d for atoms of target
        ICOMP(1-NAT,1-3)=1 for sites within inner ellipsoid,
                        =2 for sites between inner and outer ellipsoids
        CDESCR=description of target (up to 67 characters)
        X0(1-3)=(location/d) in Target Frame corresponding to dipole with IXYZ=(0,0,0).  
				This will be treated at the origin of physical coordinates in the TF.
                Here set to be centroid of the cube and enclosed ellipsoid 

Fortran history records removed.

Copyright (C) 2011 B.T. Draine and P.J. Flatau
Copyright (c) 2012, C++ version, Choliy V.
This code is covered by the GNU General Public License.
** */

Tarrctell::Tarrctell(TargetManager *man) : AbstractTarget(man) 
{ 
	shortDescr = string("Tarrctell");
	longDescr = string("in rect elipsoid");
	xoff = yoff = zoff = (real)0.;
}

void Tarrctell::Sizer(void)
{
	dx = manager->CashedDx();
//
// Check that input parameters have "inner" ellipsoid smaller than "outer" rectangular volume:
	if (shpar[0] > shpar[3]) 
		Errmsg("Fatal", shortDescr.c_str(), " ax > bx ");
	if (shpar[1] > shpar[4]) 
		Errmsg("Fatal", shortDescr.c_str(), " ay > by ");
	if (shpar[2] > shpar[5]) 
		Errmsg("Fatal", shortDescr.c_str(), " az > bz ");
//
	int jx = nint_(shpar[3]);
	xoff = jx % 2 ? zero_ : half_;
	int jy = nint_(shpar[4]);						// Similar criterion for YOFF:
	yoff = jy % 2 ? zero_ : half_;
	int jz = nint_(shpar[5]);						// Similar criterion for ZOFF:
	zoff = jz % 2 ? zero_ : half_;
//
// JX even: run from -JX/2 to JX/2-1      e.g. -2 to +1   if JX=4, XOFF=0.5
//    odd :          -(JX/2) to (JX/2)    e.g. -2 to +2   if JX=5
	minJx = -(jx/2);
	maxJx =  (jx/2) - nint_(2 * xoff);
	minJy = -(jy/2);
	maxJy =  (jy/2) - nint_(2 * yoff);
	minJz = -(jz/2);
	maxJz =  (jz/2) - nint_(2 * zoff);
	AllocateArrays(maxJx - minJx + 1, maxJy - minJy + 1, maxJz - minJz + 1);
}

void Tarrctell::Descriptor(void)			// Write target description into string CDESCR
{
	if((shpar[0] == shpar[1]) && (shpar[0] == shpar[2]))
		sprintf(freeDescr, " Spherical with%7d dipoles", nin);
	else
		sprintf(freeDescr, " Ellipsoid with%7d dipoles", nin);
}

void Tarrctell::Vector(void)
{
	VectorA();
}

void Tarrctell::VectorX(void)				// Locate centroid = origin of coordinates in TF
{
	x0.Clear();
	for(int jy=0; jy<nat0; ++jy)
	{
		x0 += Vect3<real>((real)ixyz.Value(jy, 0), (real)ixyz.Value(jy, 1), (real)ixyz.Value(jy, 2)); 
	}
	x0 /= (real)nat0;
}

void Tarrctell::Allocator(void)
{
//
// Current version of TARRCTELL is restricted to cubic lattices
	if((dx.data[0] != onex_) || (dx.data[1] != onex_))
		Errmsg("Fatal", "Tarrctell", " tarrctell does not support noncubic lattice");
//
// Inner ellipsoid is defined by (X/AX)**2+(Y/AY)**2+(Z/AZ)**2 = 0.25
// Outer rectangular volume has sides BX,BY,BZ
// Material within ellipsoid is of composition 1
// Material between ellipsoids and outer surface is of composition 2
//
// Dipoles are located at sites
// (x,y,z)=(I+XOFF,J+YOFF,Z+KOFF), I,J,K=integers XOFF,YOFF,ZOFF=constants
	nat0 = 0;
	nin = 0;
	int curSize = nx * ny;
	ixyz.Dimension(curSize, 3);
//
// Determine list of occupied sites.
	real ax2 = shpar[0] * shpar[0];
	real ay2 = shpar[1] * shpar[1];
	real az2 = shpar[2] * shpar[2];
	for(int jx=minJx; jx<=maxJx; ++jx)
	{
		real x = jx + xoff;
		real rx2 = x * x / ax2;
		for(int jy=minJy; jy<=maxJy; ++jy)
		{
			real y = jy + yoff;
			real rxy2 = rx2 + y * y / ay2;
			for(int jz=minJz; jz<=maxJz; ++jz)
			{
				real z = jz + zoff;
				real r = rxy2 + z * z / az2;
				if (nat0 >= curSize)
				{
					ixyz.Extend(nz);
					curSize = ixyz.GetSize(0);
				}
				ixyz.Fill3(nat0, jx, jy, jz);
				int index = GetLinearAddress(nat0);
				if (r < quat_)
				{
					Composer(index, 0);
					++nin;
				}
				else
					Composer(index, 1);
				iocc[index] = true;
				++nat0;
			}
		}
	}
	Descriptor();
}

const char *TargetVerboseDescriptor_Elinrct(int num)
{
	static const char *descr[6] = 
	{
		"Ellipsoid length/d in x", "Ellipsoid length/d in y", "Ellipsoid length/d in z", 
		"Length of larger rectangular volume in x", "Length of larger rectangular volume in y", "Length of larger rectangular volume in z"
	};
	if (num < 6)
		return descr[num];
	else
		return NULL;
}

REGISTER_TARGET(Elinrct,6,false,-1,0,"Ellipsiod embedded in rectangular volume")
void Target_Elinrct::SayHello(FILE *stream)
{
	fprintf(stream, "The target %s parameters are:\n", shortDescr.c_str());
	fprintf(stream, "ellipsoid diameters xv,yv,zv:\n");
	fprintf(stream, "%lf %lf %lf\n", shpar[0], shpar[1], shpar[2]);
	fprintf(stream, "rectangular block dimensions x,y,z:\n");
	fprintf(stream, "%lf %lf %lf\n", shpar[3], shpar[4], shpar[5]);
}
