#include "StdAfx.h"

#include "Tarcel.h"
#include "TargetManager.h"

/* **
Routine to construct target consisting of two materials, with outer
surface an ellipsoid of dimensions AX,AY,AZ,
and core/mantle interface a concentric ellipsoid of dimensions BX,BY,BZ
Input:
       AX=(x-length)/d of outer ellipsoid  (d=lattice spacing)
       AY=(y-length)/d "   "     "
       AZ=(z-length)/d "   "     "
       BX=(x-length)/d of inner ellipsoid
       BY=(y-length)/d "   "     "
       BZ=(z-length)/d "   "     "
       DX(1-3)=lattice spacing (dx/d,dy/d,dz/d), d=(dx*dy*dz)**(1/3)
       IOSHP=device number for "target.out" file
            =-1 to suppress printing of "target.out"
       MXNAT=dimensioning information (max number of atoms)
Output:
       A1(1-3)=(1,0,0)=unit vector defining target axis 1 in Target Frame
       A2(1-3)=(0,1,0)=unit vector defining target axis 2 in Target Frame
       NAT=number of atoms in target
       IXYZ(1-NAT,1-3)=x/d,y/d,z/d for atoms of target
       ICOMP(1-NAT,1-3)=1 for sites within inner ellipsoid,
                       =2 for sites between inner and outer ellipsoids
       CDESCR=description of target (up to 67 characters)
       X0(1-3)=(location/d) in Target Frame corresponding to dipole with
               IXYZ=(0,0,0).  This will be treated at the origin of physical
               coordinates in the TF.
               Here set to be centroid of the ellipsoids.

B.T.Draine, Princeton Univ. Obs.
Fortran history records removed.

Copyright (C) 1994,1996,1997,1998,2007,2008 B.T. Draine and P.J. Flatau
Copyright (C) 20120, 2013, C++ version, Choliy V.
This code is covered by the GNU General Public License.
** */

Tarcel::Tarcel(TargetManager *man) : AbstractTarget(man) 
{
	shortDescr = string("Tarcel");
	longDescr = string("Concentric ellipsoids");
	nin = 0;
}

void Tarcel::Sizer(void)
{
	dx = manager->CashedDx();
	AllocateArrays((int)(shpar[0] / dx.data[0] + half_), (int)(shpar[1] / dx.data[1] + half_), (int)(shpar[2] / dx.data[2] + half_));
}

void Tarcel::Descriptor(void)				// Write target description into string CDESCR
{
	if ((shpar[0] == shpar[1]) && (shpar[0] == shpar[2]))
		sprintf(freeDescr, " Spherical target containing %7d dipoles", nat0);
	else
		sprintf(freeDescr, " Concentric ellipsoids with %6d inner %7d total dipoles", nin, nat0);
}

void Tarcel::Vector(void)
{
	VectorA();
}

void Tarcel::VectorX(void)
{
//
// Locate centroid = origin of coordinates in TF
	x0.Clear();
	for(int jx=0; jx<nat0; ++jx)
	{
		x0 += Vect3<real>((real)ixyz.Value(jx, 0), (real)ixyz.Value(jx, 1), (real)ixyz.Value(jx, 2));
	}
	x0 /= (real)nat0;
}

void Tarcel::Allocator(void)
{
	const char *TarcelLabel = "Tarcel";
//
// Current version of TARCEL is restricted to cubic lattices
	if (dx.data[0] != onex_ || dx.data[1] != onex_)
		Errmsg("Fatal", TarcelLabel, " tarcel does not support noncubic lattice");
//
// Inner ellipsoid is defined by
//      (X/AX)**2+(Y/AY)**2+(Z/AZ)**2 = 0.25
// Outer ellipsoid is defined by
//      (X/BX)**2+(Y/BY)**2+(Z/BZ)**2 = 0.25
// Material within inner ellipsoid is of composition 1
// Material between inner and outer ellipsoids is of composition 2
//
// Ideal volume V=(pi/6)*BX*BY*BZ
//
// Dipoles are located at sites
// (x,y,z)=(I+XOFF,J+YOFF,Z+KOFF), I,J,K=integers XOFF,YOFF,ZOFF=constants
//
// For concentric spheres: call with AX=AY=AZ, BX=BY=BZ
// For spheroids: call with AX=AY (or AY=AZ or AX=AZ) and BX=BY (or BY=BZ or BX=BZ)
//
// Check that input parameters have "inner" ellipsoid smaller than "outer" ellipsoid:
	if (shpar[0] < shpar[3]) 
		Errmsg("Fatal", TarcelLabel, " ax < bx ");
	if (shpar[1] < shpar[4]) 
		Errmsg("Fatal", TarcelLabel, " ay < by ");
	if (shpar[2] < shpar[5]) 
		Errmsg("Fatal", TarcelLabel, " az < bz ");
//
// Criterion for choosing XOFF: try to optimize outer ellipsoid surface
// If AX is close to an even int, take XOFF=1/2 If AX is close to an odd int, take XOFF=0
	int jx = (int)(shpar[0] + half_);
	real xoff = (jx%2) ? zero_ : half_;
	int jy = (int)(shpar[1] + half_);
	real yoff = (jy%2) ? zero_ : half_;
	int jz = (int)(shpar[2] + half_);
	real zoff = (jz%2) ? zero_ : half_;
//
	minJx =-(int)(half_ * shpar[0] + half_);
	maxJx = (int)(half_ * shpar[0] - quat_);
	minJy =-(int)(half_ * shpar[1] + half_);
	maxJy = (int)(half_ * shpar[1] - quat_);
	minJz =-(int)(half_ * shpar[2] + half_);
	maxJz = (int)(half_ * shpar[2] - quat_);
	real ax2 = shpar[0] * shpar[0];
	real ay2 = shpar[1] * shpar[1];
	real az2 = shpar[2] * shpar[2];
	real bx2 = shpar[3] * shpar[3];
	real by2 = shpar[4] * shpar[4];
	real bz2 = shpar[5] * shpar[5];
	nat0 = 0;
	nin = 0;
//
// Determine list of occupied sites.
	int curSize  = nx * ny;
	ixyz.Dimension(curSize, 3);
//
	for(jx=minJx; jx<=maxJx; ++jx)
	{
		real x = jx + xoff;
		real rx2 = x*x / ax2;
		if (rx2 < quat_)
		{
			for(jy=minJy; jy<=maxJy; ++jy)
			{
				real y = jy + yoff;
				real ryz2 = rx2 + y*y / ay2;
				if (ryz2 < quat_)
				{
					for(jz=minJz; jz<=maxJz; ++jz)
					{
						real z = jz + zoff;
						real rz2 = ryz2 + z*z / az2;
						if (rz2 < quat_)							// Site is occupied:
						{
							real rr = z*z/bz2 + y*y/by2 + x*x/bx2;
							if (nat0 >= curSize)
							{
								ixyz.Extend(nz);
								curSize = ixyz.GetSize(0);
							}
							ixyz.Fill3(nat0, jx, jy, jz);
							int index = GetLinearAddress(nat0);
							Composer(index, (rr < quat_) ? 0 : 1);
							iocc[index] = true;
							++nat0;
						}
					}
				}
			}
		}
	}
	ixyz.Close(nat0);
}

void Tarcel::OutShpar(FILE *file)
{
	fprintf(file, " AX,AY,AZ=%8.4lf%8.4lf%8.4lf\n BX,BY,BZ=%8.4lf%8.4lf%8.4lf\n", shpar[0], shpar[1], shpar[2], shpar[3], shpar[4], shpar[5]);
}

const char *TargetVerboseDescriptor_Conellips(int num)
{
	static const char *descr[6] = 
	{
		"x-length/d of outer ellipsoid (icomp=2)", "y-length/d of outer ellipsoid (icomp=2)", "z-length/d of outer ellipsoid (icomp=2)",
		"x-length/d of inner ellipsoid (icomp=1)", "y-length/d of inner ellipsoid (icomp=1)", "z-length/d of inner ellipsoid (icomp=1)"
	};
	if (num < 6)
		return descr[num];
	else
		return NULL;
}

REGISTER_TARGET(Conellips,6,false,-1,2,"Two concentric ellipsoids")
void Target_Conellips::SayHello(FILE *stream)
{
	fprintf(stream, "The target %s parameters are:\n", shortDescr.c_str());
	fprintf(stream, "length of outer ellipsoid (icomp=2) in x,y,z directions (lattice units):\n");
	fprintf(stream, "%lf %lf %lf\n", shpar[0], shpar[1], shpar[2]);
	fprintf(stream, "length of inner ellipsoid (icomp=1) in x,y,z directions (lattice units):\n");
	fprintf(stream, "%lf %lf %lf\n", shpar[3], shpar[4], shpar[5]);
	if ((shpar[3] > shpar[0]) || (shpar[4] > shpar[1]) || (shpar[5] > shpar[2]))
	{
		fprintf(stream, "Input error: inner ellipsoid is not contained within outer ellipsoid.\n");
	}
}
