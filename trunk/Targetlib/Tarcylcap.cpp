#include "StdAfx.h"

#include "Tarcylcap.h"
#include "TargetManager.h"

/* **
Purpose: 
to construct cylindrical target with hemispherical caps from "atoms", with cylinder axis along x axis.
cylinder diameter is b, length of cylinder itself is a length of cylinder with caps is (a+b)

Input:
       A = cylinder length (in units of lattice spacing d) *NOT* including caps
       B = cylinder diameter (in units of lattice spacing d)
       DX(1-3)=(dx,dy,dz)/d where dx,dy,dz=lattice spacings in x,y,z directions, and d=(dx*dy*dz)**(1/3)=effective lattice spacing
       IOSHP=device number for "target.out" file =-1 to suppress printing of "target.out"
       MXNAT = dimensioning information

Returns:
       A1(1-3) = unit vector along cylinder axis
       A2(1-3) = unit vector perpendicular to cylinder axis
       CDESCR = string describing target (up to 67 chars.)
       NAT = number of atoms in target
       IXYZ(1-NAT,1-3)=x/d,y/d,z/d for atoms in target
       ICOMP(1-NAT,1-3)=1 (composition identifier)
       X0(1-3)=location/d) in Target Frame corresponding to dipole with IXYZ=(0,0,0).  
	   This will be treated at the origin of physical coordinates in the TF. Here set to be centroid of capped cylinder

Fortran history records removed.

Copyright (C) 2006,2007,2008 B.T. Draine and P.J. Flatau
Copyright (c) 2012, C++ version, V.Choliy
This code is covered by the GNU General Public License.
** */

Tarcylcap::Tarcylcap(TargetManager *man) : AbstractTarget(man) 
{
	shortDescr = string("Tarcylcap");
	longDescr = string("Cylinder with end-cups");
	xcm = ycm = zcm = (real)0.;
}

void Tarcylcap::Sizer(void)
{
	dx = manager->CashedDx();
	jxu = (int)(shpar[0] / dx.data[0] + half_);
	maxJx = (int)(half_ * (jxu + 1 + shpar[0] + shpar[1]));
	minJx = jxu + 1 - maxJx;
	jxl = 1;
	xcm = half_ * (real)(jxu + 1);

	int nb = (int)(shpar[1] / dx.data[1] + half_);
	if (nb % 2)
	{
		ycm = zero_;
		minJy = -nb/2;
		maxJy =  nb/2;
	}
	else
	{
		ycm = half_;
		minJy = -nb/2 + 1;
		maxJy =  nb/2;
	}

	nb = (int)(shpar[1] / dx.data[2] + half_);
	if (nb % 2)
	{
		zcm = zero_;
		minJz = -nb/2;
		maxJz =  nb/2;
	}
	else
	{
		zcm = half_;
		minJz = -nb/2 + 1;
		maxJz =  nb/2;
	}
	AllocateArrays(maxJx - minJx + 1, maxJy - minJy + 1, maxJz - minJz + 1);
}

void Tarcylcap::Descriptor(void)
{

}

void Tarcylcap::Vector(void)
{
	VectorA();
}

void Tarcylcap::VectorX(void)
{
//
// Now redetermine location of centroid
	Vect3<real> xyzcm;
	for(int ja=0; ja<nat0; ++ja)
	{
		xyzcm += Vect3<real>((real)ixyz.Value(ja, 0), (real)ixyz.Value(ja, 1), (real)ixyz.Value(ja, 2));
	}
	x0 = -xyzcm / (real)nat0;
}

void Tarcylcap::Allocator(void)
{
//
// A = cylinder length/d, B = cylinder diameter/d
//
// Determine centroid XCM,YCM,ZCM
//
// If A/d_x is near an odd int, XCM=0
// If A/d_x is near an even int, XCM=0.5
//
// If B/d_y is near an odd or even int, YCM=0 or 0.5
// If B/d_z is near an odd or even int, ZCM=0 or 0.5
//
// Determine limits for testing x,y,z values
// In case of disk axis, run from I=1 to I=INT(A+0.5)
// In radial directions, place atoms as follows:
//     If B is close to even number, y = j + 0.5, with j running from -int((b+.5)/2) to int((b+.5)/2)-1
//     If B is close to odd number, x = j, with j running from -int((b+.5)/2) to int((b+.5)/2)
//
	int jx, jy, jz;
	nat0 = 0;
	int curSize  = nx * ny;
	ixyz.Dimension(curSize, 3);
//
	for(jx=jxl; jx<=jxu; ++jx)
	{
		for(jy=minJy; jy<=maxJy; ++jy)
		{
			for(jz=minJz; jz<=maxJz; ++jz)
			{
				real z2 = (jz - zcm) * dx.data[2];
				z2 = z2 * z2;
				real y2m = (quat_ * shpar[1] * shpar[1] - z2) / (dx.data[1] * dx.data[1]);
				if((jy - ycm) * (jy - ycm) < y2m)
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
//
// NLAY = number of layers in cylinder
// NFAC = number of atoms in slice
	int nlay = maxJx;
	int nfac = nat0 / nlay;
//
// Now add hemispherical caps
// XCM = 0.5*REAL(NX2+1)
// XU  = XCM+0.5*(A+B)
// JXU = INT(XU)
// JXL = NX2+1-JXU
	real r2 = quat_ * shpar[1] * shpar[1];
//	int jxu = (int)(half_ * (maxJx + 1 + shpar[0] + shpar[1]));
//	int jxl = maxJx + 1 - jxu;
//
	bool bBegin = true;
	for(jx=minJx; jx<=maxJx; ++jx)
	{
		if ((jx >= jxl) && (jx <= jxu))
		{
			bBegin = false;
			continue;
		}
		real x2 = (bBegin == true) ? (jx - xcm + half_ * shpar[0]) * dx.data[0] : (jx - xcm - half_ * shpar[0]) * dx.data[0];
		x2 = x2 * x2;
		if (x2 <= r2)
		{
			for(jy=minJy; jy<=maxJy; ++jy)
			{
				real y2 = (jy - ycm) * dx.data[1];
				y2 = y2 * y2;
				if (x2 + y2 < r2)
				{
					for(jz=minJz; jz<=maxJz; ++jz)
					{
						real z2 = (jz - zcm) * dx.data[2];
						z2 = z2 * z2;
						if (x2 + y2 + z2 < r2)
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
//
// REFF2=effective radius**2 of disk
	real reff2 = (real)nfac / Pi;
//
// ASPR=aspect ratio (length/diameter)
	real aspr = half_ * (real)nlay / Sqrt(reff2);
	sprintf(freeDescr, " Cyl.prism, NAT=%7d NFAC=%4d NLAY=%4d asp.ratio=%7.4lf", nat, nfac, nlay, aspr);
}

void Tarcylcap::OutShpar(FILE *file)
{
	fprintf(file, " Length=%8.4lf Diameter=%8.4lf\n", shpar[0], shpar[1]);
}

const char *TargetVerboseDescriptor_Cylndrcap(int num)
{
	static const char *descr[2] = 
	{
		"cylinder length", "cylinder diameter",
	};
	if (num < 2)
		return descr[num];
	else
		return NULL;
}

REGISTER_TARGET(Cylndrcap,2,false,-1,0,"homogenous, isotropic cylinder with hemispherical endcaps")
void Target_Cylndrcap::SayHello(FILE *stream)
{
	fprintf(stream, " CYLNDRCAP: cylinder with hemispherical endcaps\n");
	fprintf(stream, "The target %s parameters are:\n", shortDescr.c_str());
	fprintf(stream, "cylinder length (not including end-caps) and diameter:\n");
	fprintf(stream, "%lf %lf\n", shpar[0], shpar[1]);
}
