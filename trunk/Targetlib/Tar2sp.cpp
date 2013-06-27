#include "StdAfx.h"

#include "Tar2sp.h"
#include "TargetManager.h"

/* **
Routine to construct two touching spheroids from "atoms"
   First spheroid has symmetry axis in y direction.
   Second spheroid is displaced in x direction, with symmetry
   axis in yz plane, in direction ey*cos(phi)+ez*sin(phi).
   Separation between centroids=(B_1+B_2)/2

Input:
       A_1=length/d of 1st spheroid along symm.axis  (d=lattice spacing)
       A_2=length/d of 2nd spheroid along symm.axis
       B_1=diameter/d of 1st spheroid
       B_2=diameter/d of 2nd spheroid
       PHI=angle (deg) specifying relative orientation of spheroids
       PRINAX=0. to set A1=(1,0,0),A2=(0,1,0) in Target Framee
             =1. to set A1,A2=principal axes of largest and second
                              largest moment of inertia
       DX(1-3)=(dx/d,dy/d,dz/d) where dx,dy,dz=lattice spacing in x,y,z
           directions, and d=(dx*dy*dz)**(1/3)=effective lattice spacing
       IOSHP=device number for "target.out" file
            =-1 to suppress printing of "target.out"
       MXNAT=dimensioning information (max number of atoms)
Output:
       A1(1-3)=(1,0,0)=unit vector defining target axis 1 in Target Frame
       A2(1-3)=(0,1,0)=unit vector defining target axis 2 in Target Frame
       X0(1-3)=location in TF of dipole with IXYZ=0 0 0
       NAT=number of atoms in target
       IXYZ(1-NAT,1-3)=x/d,y/d,z/d for atoms of target
       ICOMP(1-NAT,1-3)=1 for sites within first spheroid
                       =2 for sites within second spheroid
       CDESCR=description of target (up to 67 characters)

B.T.Draine, Princeton Univ. Obs.
Fortran history records removed.

Copyright (C) 1996,1997,1998,1999,2000,2003,2004,2007,2008
               B.T. Draine and P.J. Flatau
Copyright (C) 2012, 2013, C++ version, Choliy V.
This code is covered by the GNU General Public License.
** */

Tar2sp::Tar2sp(TargetManager *man) : AbstractTarget(man) 
{
	shortDescr = string("Tar2sp");
	longDescr = string("Touching spheroids");
	xc1 = yc1 = zc1 = (real)0.;
}

void Tar2sp::Sizer(void)
{
	dx = manager->CashedDx();
// Note that other criteria are obviously possible. For example, could have chosen YC1 and ZC1 to maximize number of lattice sites falling within spheroidal surfaces.
	xc1 = half_ * (onex_ - shpar[1]);
	if ((xc1 - (int)(xc1) < quat_) || (xc1 - (int)(xc1) > (real)0.75))
         yc1 = zc1 = zero_;
	else
         yc1 = zc1 = half_;
//
// Determine list of occupied sites in first spheroid
	int lmx1 = (int)(xc1 - half_ * shpar[1]);
	int lmx2 = (int)(xc1 + half_ * shpar[1]);
	int lmy1 = (int)(yc1 - half_ * shpar[0]);
	int lmy2 = (int)(yc1 + half_ * shpar[0]);
	int lmz1 = (int)(zc1 - half_ * shpar[1]);
	int lmz2 = (int)(zc1 + half_ * shpar[1]);
//
	real xc2 = xc1 + half_ * (shpar[1] + shpar[3]);
	minJx = min_(lmx1, (int)(xc2 - half_ * shpar[3]));
	maxJx = max_(lmx2, (int)(xc2 + half_ * shpar[3]));
	minJy = min_(lmy1, (int)(yc1 - half_ * max(shpar[2], shpar[3])));
	maxJy = max_(lmy2, (int)(yc1 + half_ * max(shpar[2], shpar[3])));
	minJz = min_(lmz1, (int)(zc1 - half_ * max(shpar[2], shpar[3])));
	maxJz = max_(lmz2, (int)(zc1 + half_ * max(shpar[2], shpar[3])));
	AllocateArrays(maxJx - minJx + 1, maxJy - minJy + 1, maxJz - minJz + 1);
}

void Tar2sp::Descriptor(void)						// Write target description into string CDESCR
{
	sprintf(freeDescr, " Spheroids with%7d +%7d=%7d dipoles", nat1, nat2, nat0);
}

void Tar2sp::Vector(void)
{
//
// Specify target axes A1 and A2
// If PRINAX=0. then
//    A1=(1,0,0) in target frame
//    A2=(0,1,0) in target frame
// If PRINAX=1. then A1,A2 are principal axes of largest,second largest moment of inertia
	if (shpar[5] <= (real)0.)
		VectorA();
	else
		Prinaxis();
}

void Tar2sp::VectorX(void)					// Determine X0(1-3)=TF coordinates/d of dipole with IXYZ=0 0 0
{
	x0.Set(-half_, -yc1, -zc1);
}

void Tar2sp::Allocator(void)
{
//
// Current version of TAR2SP is restricted to cubic lattices
	if ((dx.data[0] != onex_) || (dx.data[1] != onex_))
		Errmsg("Fatal", "tar2sp", " tar2sp does not yet support noncubic lattice");
//
// Routine to construct target array representing two spheroids in contact.
// line connecting spheroid centers is parallel to x
// symmetry axis of first spheroid is parallel to y
// symmetry axis of second spheroid is in y-z plane, at angle PHI to y
//
// First spheroid:
// Array sites contained within spheroidal surface defined by ((X-XC1)/B_1)**2+((Y-YC1)/A_1)**2+((Z-ZC1)/B_1)**2=0.25
// Ideal volume V=(pi/6)*A_1*B_1*B_1
//
// Second spheroid:
//    ((X-XC2)/B_2)**2+((U-UC1)/A_1)**2+((V-VC1)/B_1)**2=0.25
// where
//    U  =  Y*cos(PHI)  +  Z*sin(PHI)
//    V  = -Y*sin(PHI)  +  Z*cos(PHI)
//   UC1 = YC1*cos(PHI) + ZC1*sin(PHI)
//   VC1 =-YC1*sin(PHI) + ZC1*cos(PHI)
//
// Dipoles are located at sites (x,y,z)=(I,J,K)*d, I,J,K=integers
//
// Criterion for choosing XC1:
// Point where spheroids contact each other should be midway between lattice points.  
// Thus take XC1+0.5*B_1 = 0.5 or XC1 = 0.5(1.-B_1)
//
// Criterion for choosing YC1 and ZC1:
// If XC1 is close to an int, take YC1=ZC1=0
// If XC1 is close to a half-int, take YC1=ZC1=0.5
//
// Note that other criteria are obviously possible. For example, could have chosen YC1 and ZC1 to maximize number of lattice sites falling within spheroidal surfaces.
	nat0 = 0;
//
// Determine list of occupied sites in first spheroid
	int lmx1 = (int)(xc1 - half_ * shpar[1]);
	int lmx2 = (int)(xc1 + half_ * shpar[1]);
	int lmy1 = (int)(yc1 - half_ * shpar[0]);
	int lmy2 = (int)(yc1 + half_ * shpar[0]);
	int lmz1 = (int)(zc1 - half_ * shpar[1]);
	int lmz2 = (int)(zc1 + half_ * shpar[1]);
	real ax2 = shpar[1] * shpar[1];
	real ay2 = shpar[0] * shpar[0];
	real az2 = shpar[1] * shpar[1];

	int curSize = nx * ny;
	ixyz.Dimension(curSize, 3);

	real x, y, z, r, rx2, rxy2;
	int jx, jy, jz;
	for(jx=lmx1; jx<=lmx2; ++jx)
	{
		x = jx - xc1;
		rx2 = x * x / ax2;
		if (rx2 < quat_)
		{
			for(jy=lmy1; jy<=lmy2; ++jy)
			{
				y = jy - yc1;
				rxy2 = rx2 + y * y / ay2;
				if (rxy2 < quat_)
				{
					for(jz=lmz1; jz<=lmz2; ++jz)
					{
						z = jz - zc1;
						r = rxy2 + z * z / az2;
						if (r < quat_)
						{
							if (nat0 >= curSize)
							{
								ixyz.Extend(nz);
								curSize = ixyz.GetSize(0);
							}
							ixyz.Fill3(nat0, jx, jy, jz);							// Positions are expressed starting from one, not zero 
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
	nat1 = nat0;
//
// Determine occupied sites in second spheroid
    real xc2 = xc1 + half_ * (shpar[1] + shpar[3]);
    real yc2 = yc1;
    real zc2 = zc1;
	lmx1 = (int)(xc2 - half_ * shpar[3]);
	lmx2 = (int)(xc2 + half_ * shpar[3]);
//
// This is not time-consuming, so no need to optimize choices of LMY1,LMY2,LMZ1,LMZ2
	lmy1 = (int)(yc1 - half_ * max(shpar[2], shpar[3]));
	lmy2 = (int)(yc1 + half_ * max(shpar[2], shpar[3]));
	lmz1 = (int)(zc1 - half_ * max(shpar[2], shpar[3]));
	lmz2 = (int)(zc1 + half_ * max(shpar[2], shpar[3]));
	real sinphi = (real)Sin(Pi * shpar[4] / 180.);
	real cosphi = (real)Cos(Pi * shpar[4] / 180.);
//
// Transform (y,z) -> (u,v) with u=y*cos(phi)-z*sin(phi) v=y*sin(phi)+z*cos(phi)
	for(jx=lmx1; jx<=lmx2; ++jx)
	{
		x = jx - xc2;
		real rx2 = x / shpar[3];
		rx2 = rx2 * rx2;
		if (rx2 < quat_)
		{
			for(jy=lmy1; jy<=lmy2; ++jy)
			{
				y = jy - yc2;
				for(jz=lmz1; jz<=lmz2; ++jz)
				{
					z = jz - zc2;
					real u = (cosphi*y - sinphi*z) / shpar[2];
					real v = (sinphi*y + cosphi*z) / shpar[3];
					r = rx2 + u*u + v*v;
					if (r < quat_)										// ! Site is occupied:
					{
						if (nat0 >= curSize)
						{
							ixyz.Extend(nz);
							curSize = ixyz.GetSize(0);
						}
						ixyz.Fill3(nat0, jx, jy, jz);
						int index = GetLinearAddress(nat0);
						Composer(index, 1);
						iocc[index] = true;
						++nat0;
					}
				}
			}
		}
	}
	ixyz.Close(nat0);
	nat2 = nat0 - nat1;
}

void Tar2sp::OutShpar(FILE *file)
{
	fprintf(file, "A_1,B_2,A_2,B_2,PHI, PRINAXIS= %8.flf%8.flf%8.flf%8.flf%8.flf%d\n", shpar[0], shpar[1], shpar[2], shpar[3], shpar[4], (int)shpar[5]);
}

const char *TargetVerboseDescriptor_Sphroid2(int num)
{
	static const char *descr[6] = 
	{
		"length of first spheroid", "diameter of first spheroid", "length of second spheroid", "diameter of second spheroid", 
		"angle phi (deg) between axes a_1 and a_2", "PRINAX=0 for a_1=(1,0,0),a_2=(0,1,0) in TF, =1 for a_1,a_2 = principal axes"
	};
	if (num < 6)
		return descr[num];
	else
		return NULL;
}

REGISTER_TARGET(Sphroid2,6,false,-1,2,"Two touching spheroids (with angle phi between symm. axes)")
void Target_Sphroid2::SayHello(FILE *stream)
{
	fprintf(stream, "The target %s parameters are:\n", shortDescr.c_str());
	fprintf(stream, "length a_1 and diameter b_1 of 1st spheroid:\n");
	fprintf(stream, "%lf %lf\n", shpar[0], shpar[1]);
	fprintf(stream, "length a_2 and diameter b_2 of 2nd spheroid:\n");
	fprintf(stream, "%lf %lf\n", shpar[2], shpar[3]);
	fprintf(stream, "angle phi (deg) between axes a_1 and a_2:\n");
	fprintf(stream, "%lf\n", shpar[4]);
	fprintf(stream, "PRINAX=0 for a_1=(1,0,0),a_2=(0,1,0) in TF, =1 for a_1,a_2 = principal axes\n");
	fprintf(stream, "%lf\n", shpar[5]);
}
