#include "StdAfx.h"

#include "Tarell.h"
#include "TargetManager.h"

/* **
Routine to construct ellipsoid from "atoms"
Input:
       AX=(x-length)/d    (d=lattice spacing)
       AY=(y-length)/d
       AZ=(z-length)/d
       DX(1-3)=(dx,dy,dz)/d where dx,dy,dz=lattice spacings in x,y,z
               directions, and d=(dx*dy*dz)**(1/3)=effective lattice spacing
       IOSHP=device number for "target.out" file
            =-1 to suppress printing of "target.out"
       MXNAT=dimensioning information (max number of atoms)
Output:
       A1(1-3)=(1,0,0)=unit vector defining target axis 1 in Target Fr
       A2(1-3)=(0,1,0)=unit vector defining target axis 2 in Target Fr
       NAT=number of atoms in target
       IXYZ(1-NAT,1-3)=(x-x0(1))/d,(y-x0(2))/d,(z-x0(3))/d for atoms of target
       CDESCR=description of target (up to 67 characters)
       ICOMP(1-NAT,1-3)=1 (composition identifier)
       X0(1-3)=(location/d) in Target Frame corresponding to dipole with
               IXYZ=(0,0,0).  This will be treated at the origin of physical
               coordinates in the TF.
               Here origin is set to be centroid of ellipsoid.

B.T.Draine, Princeton Univ. Obs.
Fortran history records removed.

Copyright (C) 1993,1996,1997,1998,2000,2007,2008 
              B.T. Draine and P.J. Flatau
Copyright (C) 2012, 2013, C++ version, V.Choliy.
This code is covered by the GNU General Public License.
** */

Tarell::Tarell(TargetManager *man) : AbstractTarget(man)
{ 
	shortDescr = string("Tarell");
	longDescr = string("Ellipsoidal grain");
}

void Tarell::Sizer(void)
{
	dx = manager->CashedDx();
	AllocateArrays((int)(shpar[0] / dx.data[0] + half_), (int)(shpar[1] / dx.data[1] + half_), (int)(shpar[2] / dx.data[2] + half_));
}

void Tarell::Descriptor(void)
{
//
// Write target description into string CDESCR
	const char *sphere = "Sphere";
	const char *ellips = "Ellipsoid";
	const char *FormatX = " %s,%7d dipoles, %6.3lf%6.3lf%6.3lf=x,y,z lattice spacing";
	const char *what = ellips;
	if ((shpar[0] == shpar[1]) && (shpar[0] == shpar[2]))
		what = sphere;
	sprintf(freeDescr, FormatX, what, nat0, dx.data[0], dx.data[1], dx.data[2]);
}

void Tarell::Vector(void)
{
	VectorA();
}

/* **
Routine to construct pseudo-ellipsoidal target aray.
With occupied array sites contained within ellipsoidal surface
defined by (X/AX*d)**2+(Y/AY*d)**2+(Z/AZ*d)**2=0.25
Ideal volume V=(pi/6)*AX*AY*AZ*d**3
where d = effective lattice spacing

Dipoles are located on lattice at sites
(x,y,z)=(I+XOFF,J+YOFF,Z+KOFF), I,J,K=integers XOFF,YOFF,ZOFF=constants

For sphere: call with AX=AY=AZ
For spheroid: call with AX=AY (or AY=AZ or AX=AZ)
B.T.Draine, Princeton Univ. Obs., 88.08.12

Criterion for choosing XOFF:
If AX is close to an even integer, take XOFF=1/2
If AX is close to an odd integer, take XOFF=0
** */
void Tarell::Allocator(void)
{
	int jx = (int)(shpar[0] / dx.data[0] + half_);
	real xoff = (jx % 2) ? zero_ : half_;
	int jy = (int)(shpar[1] / dx.data[1] + half_);
	real yoff = (jy % 2) ? zero_ : half_;
	int jz = (int)(shpar[2] / dx.data[2] + half_);
	real zoff = (jz % 2) ? zero_ : half_;
//
	minJx =-(int)(half_ * shpar[0] / dx.data[0] + half_);
	maxJx = (int)(half_ * shpar[0] / dx.data[0] - quat_);
	minJy =-(int)(half_ * shpar[1] / dx.data[1] + half_);
	maxJy = (int)(half_ * shpar[1] / dx.data[1] - quat_);
	minJz =-(int)(half_ * shpar[2] / dx.data[2] + half_);
	maxJz = (int)(half_ * shpar[2] / dx.data[2] - quat_);
	real ax2 = shpar[0] * shpar[0];
	real ay2 = shpar[1] * shpar[1];
	real az2 = shpar[2] * shpar[2];
	nat0 = 0;
//
// Determine list of occupied sites.
	int curSize  = nx * ny;
	ixyz.Dimension(curSize, 3);

	for(jx=minJx; jx<=maxJx; ++jx)
	{
		real x = (jx + xoff) * dx.data[0];
		real rx2 = x*x / ax2;
		if (rx2 < quat_)
		{
			for(jy=minJy; jy<=maxJy; ++jy)
			{
				real y = (jy + yoff) * dx.data[1];
				real rxy2 = rx2 + y*y / ay2;
				if (rxy2 < quat_)
				{
					for(jz=minJz; jz<=maxJz; ++jz)
					{
						real z = (jz + zoff) * dx.data[2];
						real r = rxy2 + z*z / az2;
						if (r < quat_)
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
		}
	}
	ixyz.Close(nat0);
	Descriptor();
}

void Tarell::VectorX(void)
{
	x0.Clear();
	for(int jx=0; jx<nat0; ++jx)
	{
		x0 += Vect3<real>((real)ixyz.Value(jx, 0), (real)ixyz.Value(jx, 1), (real)ixyz.Value(jx, 2));
	}
	x0 = -x0 / (real)nat0;
}

void Tarell::OutShpar(FILE *file)
{
	fprintf(file, " AX,AY,AZ=%8.4lf%8.4lf%8.4lf\n", shpar[0], shpar[1], shpar[2]); 
}

const char *TargetVerboseDescriptor_Aniellips(int num)
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

void Target_Aniellips::Composer(int index, int item)
{
	AnisotropicComposer(index, item);
}

REGISTER_TARGET(Aniellips,3,false,-1,3,"Anisotropic ellipsoid")
void Target_Aniellips::SayHello(FILE *stream)
{
	fprintf(stream, "The target = %s parameters are:\n", shortDescr.c_str());
	fprintf(stream, "x-length/d, y-length/d, z-length/d:\n");
	fprintf(stream, "%lf %lf %lf\n", shpar[0], shpar[1], shpar[2]);
}

void Target_Aniellips::PrepareIaniso(void)
{
	ianiso = TargetIsAnisotropic;
}

const char *TargetVerboseDescriptor_Ellipsoid(int num)
{
	static const char *descr[3] = 
	{
		"Ellipsoid diameters xv", "Ellipsoid diameters yv", "Ellipsoid diameters zv"
	};
	if (num < 3)
		return descr[num];
	else
		return NULL;
}

REGISTER_TARGET(Ellipsoid,3,false,-1,0,"Homogenous, isotropic ellipsoid")
void Target_Ellipsoid::SayHello()
{
	fprintf(stdout, "The target %s parameters are:\n", shortDescr.c_str());
	fprintf(stdout, "ellipsoid diameters xv,yv,zv:\n");
	fprintf(stdout, "%lf %lf %lf\n", shpar[0], shpar[1], shpar[2]);
}
