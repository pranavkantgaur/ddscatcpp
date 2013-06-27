#include "StdAfx.h"

#include "TarNel.h"
#include "TargetManager.h"

/* **
Routine to construct N (un)touching ellipsoids from "atoms"
Input:
       AX=(x-length of one ellipsoid)/d    (d=lattice spacing)
       AY=(y-length of one ellipsoid)/d
       AZ=(z-length of one ellipsoid)/d
       DX(1-3)=x,y,z lattice spacing/d
       IOSHP=device number for "target.out" file
            =-1 to suppress printing of "target.out"
       MXNAT=dimensioning information (max number of atoms)
Output:
       A1(1-3)=(1,0,0)=unit vector defining target axis 1 in Target Fr
       A2(1-3)=(0,1,0)=unit vector defining target axis 2 in Target Fr
       NAT=number of dipoles in target
       IXYZ(1-NAT,1-3)=x/d,y/d,z/d for dipole locations
       ICOMP(1-NAT,1-3)=x,y,z composition at dipole locations
                       = 1 for locations in first ellipsoid
                         2 for locations in second ellipsoid
       CDESCR=description of target (up to 67 characters)
       X0(1-3)=(location/d) in TF of dipole with IXYZ=(0,0,0)
               This determines physical origin of coordinates of TF.
               We set origin to be at point symmetrically between the
               two ellipses.

Copyright (C) 2013, C++ version, Choliy V.
This code is covered by the GNU General Public License.
** */

TarNel::TarNel()
{
	shortDescr.clear();
	longDescr.clear();
	mySize = 0;
	delta = 0;
}

TarNel::~TarNel()
{

}

TarNel::TarNel(TargetManager *man) : AbstractTarget(man) 
{
	shortDescr = string("TarNel");
	longDescr = string("Set of ellipsoidal grain");
	mySize = (int)man->GetShpar(3);
	delta = (int)man->GetShpar(4);
}

void TarNel::Sizer(void)
{
	dx = manager->CashedDx();
	AllocateArrays(mySize * (int)(shpar[0] / dx.data[0] + half_) + (mySize - 1) * delta, (int)(shpar[1] / dx.data[1] + half_), (int)(shpar[2] / dx.data[2] + half_));
}

void TarNel::Vector(void)
{
	VectorA();
}

void TarNel::VectorX(void)				// Now find X0, assuming TF origin to be at centroid
{
	x0.Clear();
	for(int jx=0; jx<nat0; ++jx)
	{
		x0 += Vect3<real>((real)ixyz.Value(jx, 0), (real)ixyz.Value(jx, 1), (real)ixyz.Value(jx, 2));
	}
	x0 = -x0 / (real)nat0;
}

void TarNel::OutShpar(FILE *file)
{
	fprintf(file, " AX,AY,AZ=%8.4lf%8.4lf%8.4lf mySize=%d delta=%d\n", shpar[0], shpar[1], shpar[2], (int)shpar[3], (int)shpar[4]); 
}

void TarNel::Descriptor(void)						// Write target description into string CDESCR
{
	const char *sphere = "Sphere";
	const char *ellips = "Ellipsoid";
	const char *what = ellips;
	const char *FormatX = " %s %ss, each containing %7d dipoles";
	if ((shpar[0] == shpar[1]) && (shpar[0] == shpar[2]))
		what = sphere;
	sprintf(freeDescr, FormatX, NumberInWords(mySize), what, nat0/mySize);
}

void TarNel::Allocator(void)
{
	if ((dx.data[0] != onex_) || (dx.data[1] != onex_))
		Errmsg("Fatal", shortDescr.c_str(), " the code does not support noncubic lattice");
//
// Routine to construct pseudo-ellipsoidal target aray.
// With occupied array sites contained within ellipsoidal surface defined by (X/AX)**2+(Y/AY)**2+(Z/AZ)**2=0.25
// Ideal volume V=(pi/6)*AX*AY*AZ
//
// Dipoles are located at sites
// (x,y,z)=(I+XOFF,J+YOFF,Z+KOFF), I,J,K=integers XOFF,YOFF,ZOFF=constants
//
// For sphere: call with AX=AY=AZ
// For spheroid: call with AX=AY (or AY=AZ or AX=AZ)
//
// Criterion for choosing XOFF:
// If AX is close to an even int, take XOFF=1/2
// If AX is close to an odd int, take XOFF=0
//
// ChB: mySize is the amount of ellipsoids in line
	int jx = (int)(shpar[0] + half_);
	real xoff = (jx % 2) ? zero_ : half_;
	int jy = (int)(shpar[1] + half_);
	real yoff = (jy % 2) ? zero_ : half_;
	int jz = (int)(shpar[2] + half_);
	real zoff = (jz % 2) ? zero_ : half_;
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
	nat0 = 0;
//
// Determine list of occupied sites in first ellipsoid
	int curSize  = nx * ny;
	ixyz.Dimension(curSize, 3);
	int jxmax = minJx;
	int jxmin = maxJx;
	for(jx=minJx; jx<=maxJx; ++jx)
	{
		real x = jx + xoff;
		real rx2 = x*x / ax2;
		if (rx2 < quat_)
		{
			for(jy=minJy; jy<=maxJy; ++jy)
			{
				real y = jy + yoff;
				real rxy2 = rx2 + y*y / ay2;
				if (rxy2 < quat_)
				{
					for(jz=minJz; jz<=maxJz; ++jz)
					{
						real z = jz + zoff;
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
		if (jx > jxmax) 
			jxmax = jx;
		if (jx < jxmin)
			jxmin = jx;
	}
	ixyz.Close(mySize * nat0);
//
// Now create duplicate second ...  ellipsoids
	int ii, ja;
	for(ii=1; ii<mySize; ++ii)
	{
		ixyz.CopyItself(ii * nat0, 0, nat0);
	}
	jxmax = jxmax - jxmin + 1;
	for(ii=1; ii<mySize; ++ii)
	{
		int shift = ii * nat0;
		int deltt = ii * (jxmax + delta);
		for(ja=0; ja<nat0; ++ja)
		{
			ixyz.Value(ja + shift, 0) += deltt;
			int index = GetLinearAddress(ja + shift);
			Composer(index, ii);
			iocc[index] = true;
		}
	}
	nat0 *= mySize;
	Descriptor();
//
	VectorX();
}

const char *TargetVerboseDescriptor_AniEllN(int num)
{
	static const char *descr[5] = 
	{
		"Length of ellipsoid in x direction", "Length of ellipsoid in y direction", "Length of ellipsoid in z direction", "Number of ellipsoids", "Distance between surfaces along x"
	};
	if (num < 5)
		return descr[num];
	else
		return NULL;
}

REGISTER_TARGET(AniEllN,5,false,-1,0,"Set of (un)touching identical anisotropic ellipsoids")
void Target_AniEllN::SayHello(FILE *stream)
{
	fprintf(stream, "The target %s parameters are:\n", shortDescr.c_str());
	fprintf(stream, "Length of one ellipsoid in x, y, z directions are:\n");
	fprintf(stream, "%lf %lf %lf\n", shpar[0], shpar[1], shpar[2]);
	fprintf(stream, "Number of ellipsoids and distance between their surfacex\n");
	fprintf(stream, "%lf %lf\n", shpar[3], shpar[4]);
}

void Target_AniEllN::Composer(int index, int item)
{
	AnisotropicComposer(index, item);
}

void Target_AniEllN::PrepareIaniso(void)
{
	ianiso = TargetIsAnisotropic;
}

const char *TargetVerboseDescriptor_EllipsoN(int num)
{
	static const char *descr[5] = 
	{
		"Length of ellipsoid in x direction", "Length of ellipsoid in y direction", "Length of ellipsoid in z direction", "Number of ellipsoids", "Distance between surfaces along x"
	};
	if (num < 5)
		return descr[num];
	else
		return NULL;
}

REGISTER_TARGET(EllipsoN,5,false,-1,0,"Set of (un)touching ellipsoids aligned along x-direction")
void Target_EllipsoN::SayHello(FILE *stream)
{
	fprintf(stream, "The target %s parameters are:\n", shortDescr.c_str());
	fprintf(stream, "x,y,z length/d for one ellipsoid:\n");
	fprintf(stream, "%lf %lf %lf\n", shpar[0], shpar[1], shpar[2]);
	fprintf(stream, "Number of ellipsoids and distance between their surfacex\n");
	fprintf(stream, "%lf %lf\n", shpar[3], shpar[4]);
}

void Target_EllipsoN::Composer(int index, int item)
{
	IsotropicComposer(index, item);
}
