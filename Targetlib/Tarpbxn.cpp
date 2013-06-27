#include "StdAfx.h"

#include "Tarpbxn.h"
#include "TargetManager.h"

/* **
Purpose: to construct "pillbox" target: disk on top of a bilayer slab disk has composition 1
	slab has compositions 2 and 3: 2 above, 3 below disk-slab interface is in the y-z plane (x=0) (in TF)
	disk symmetry axis is in x-direction (in TF)

Input:
        XD=(disk thickness)/d
        YD=(disk diameter)/d
        XS2=(x thickness of composition 2 of slab)/d   (d=lattice spacing)
        XS3=(x thickness of composition 3 of slab)/d   
        YS=(y length of slab)/d
        ZS=(z length of slab)/d
        DX(1-3)=(dx,dy,dz)/d where dx,dy,dz=lattice spacings in x,y,z directions, and d=(dx*dy*dz)**(1/3)=effective lattice spacing
        IOSHP=device number for "target.out" file =-1 to suppress printing of "target.out"
        MXNAT=dimensioning information (max number of atoms)

Output:
       A1(1-3)=unit vector (1,0,0) defining target axis 1 in Target Frame
       A2(1-3)=unit vector (0,1,0) defining target axis 2 in Target Frame
       NAT=number of atoms in target
       IXYZ(1-NAT,1-3)=x/d,y/d,z/d for atoms of target
       X0(1-3)=offset vector, such that x/d=IXYZ(J,1)+X0(1) y/d=IXYZ(J,2)+X0(2) z/d=IXYZ(J,3)+X0(3)
		where (x,y,z)=location of dipole in Target Frame with (0,0,0)=point where axis of disk
		intersects upper surface of slab upper surface of slab is assumed to located 0.5*d above dipoles defining slab

Target consists of a rectangular slab of composition 1
       extent XS x YS x ZS
       center of upper surface of slab is at (0,0,0)
       upper surface of slab is in the x=0 plane
       lower surface of slab is in the x=-XS plane
       slab runs from -XS   to 0
                      -YS/2 to +YS/2
                      -ZS/2 to +ZS/2
       plus a disk of composition 2, diameter YD, thickness XD
       with symmetry axis in x direction
       centroid of disk located at (XD/2,0,0)

       XS2,XS3,YS,ZS should ideally be integers
       if YS is even, lattice will be located at y = +/-0.5, +/-1.5, . odd y = 0, +/-1, +/-2, ..
       if ZS is even, lattice will be located at z = +/-0.5, +/-1.5, . odd z = 0, +/-1, +/-2, ..

Fortran history records removed.

Copyright (C) 2006,2007,2008 B.T. Draine and P.J. Flatau
Copyright (C) 2012, 2013, C++ version, Choliy V.
This code is covered by the GNU General Public License.
** */

Tarpbxn::Tarpbxn(TargetManager *man) : AbstractTarget(man) 
{
	shortDescr = string("Tarpbxn");
}

void Tarpbxn::Sizer(void)
{
	dx = manager->CashedDx();
	VectorX();
//
	const real half = (real)0.5;
	real r2 = quat_ * shpar[1] * shpar[1];
	minJx =  ihuge_;
	maxJx = -ihuge_;
	minJy =  ihuge_;
	maxJy = -ihuge_;
	minJz =  ihuge_;
	maxJz = -ihuge_;
//
    int jx, jy, jz;
	int jxmin = 0;
	int jxmax = nint_(shpar[0]) - 1;
	int jymin = nint_(-half * shpar[1] - x0.data[1]);
	int jymax = nint_( half * shpar[1] - x0.data[1]);
	int jzmin = nint_(-half * shpar[1] - x0.data[2]);
	int jzmax = nint_( half * shpar[1] - x0.data[2]);

	for(jx=jxmin; jx<=jxmax; ++jx)
	{
		for(jy=jymin; jy<=jymax; ++jy)
		{
			real y = jy + x0.data[1];
			y = y * y;
			for(jz=jzmin; jz<=jzmax; ++jz)
			{
				real z = jz + x0.data[2];
				z = z * z;
				if (y + z < r2)
				{
					InternalMinMax(jx, jy, jz);
				}
			}
		}
	}
//
// Completed definition of disk
// begin definition of slab
// layer 1: extends from     -xs1   to  0         [composition 2] dipoles from    JXMIN1  to -1
// layer 2: extends from -(xs1+xs2) to -xs1       [composition 3] dipoles from   JXMIN2   to (JXMIN1-1)
	jymin = nint_(-half_ * shpar[4] - x0.data[1]);
	jymax = nint_( half_ * shpar[4] - x0.data[1]);
	jzmin = nint_(-half_ * shpar[5] - x0.data[2]);
	jzmax = nint_( half_ * shpar[5] - x0.data[2]);

	int jxmin2 = -nint_(shpar[2] + shpar[3]);
	int jxmax2 = -nint_(shpar[2]) - 1;
	int jxmin1 = -nint_(shpar[2]);
	int jxmax1 = -1;
//
// layer 1 = composition 2:
	if (jxmin1 <= jxmax1)
	{
		for(jx=jxmin1; jx<=jxmax1; ++jx)
		{
			for(jy=jymin; jy<=jymax; ++jy)
			{
				if (Fabs(jy + x0.data[1]) < half_ * shpar[4])
				{
					for(jz=jzmin; jz<=jzmax; ++jz)
					{
						if (Fabs(jz + x0.data[2]) < half_ * shpar[5])
						{
							InternalMinMax(jx, jy, jz);
						}
					}
				}
			}
		}
	}
//
	if (jxmin2 <= jxmax2)
	{
		for(jx=jxmin2; jx<=jxmax2; ++jx)
		{
			for(jy=jymin; jy<=jymax; ++jy)
			{
				if (Fabs(jy + x0.data[1]) < half_ * shpar[4])
				{
					for(jz=jzmin; jz<=jzmax; ++jz)
					{
						if ((Fabs(jz + x0.data[2])) < half_ * shpar[5])
						{
							InternalMinMax(jx, jy, jz);
						}
					}
				}
			}
		}
	}
// 
	AllocateArrays(maxJx - minJx + 1, maxJy - minJy + 1, maxJz - minJz + 1);
}

void Tarpbxn::Descriptor(void)
{
	if ((shpar[2] > zero_) && (shpar[3] > zero_))
	{
		longDescr = string(" Two-layer slab+disk");
		sprintf(freeDescr, " Bilayer slab+disk with NAT = %7d dipoles", nat0);
	}
	else
	{
		longDescr = string(" Slab+disk");
		sprintf(freeDescr, " Slab+disk with NAT = %7d dipoles", nat0);
	}
}

void Tarpbxn::Vector(void)
{
	VectorA();
}

void Tarpbxn::VectorX(void)
{
//
// Target reference point = intersection of disk axis with upper surface of slab: (x,y,z)=(0,0,0) in Target Frame
// Set lattice offset in x direction: Dipole with IX=0 is located at x/d=0.5
	x0.data[0] = half_;
//
// Determine lattice offset in y and z directions
//   if NY is even, dipole with IY=0 is located at y/d=0.5 odd 0
//   if NZ is even, dipole with IZ=0 is located at z/d=0.5 odd 0
	x0.data[1] = nint_(shpar[4]) % 2 ? zero_ : half_;
	x0.data[2] = nint_(shpar[5]) % 2 ? zero_ : half_;
}

//
// ChB first init slabs and them disk
void Tarpbxn::Allocator(void)
{
// layer 1: extends from     -xs1   to  0         [composition 2] dipoles from    JXMIN1  to -1
// layer 2: extends from -(xs1+xs2) to -xs1       [composition 3] dipoles from   JXMIN2   to (JXMIN1-1)

	nat0 = 0;
	int curSize  = nx * ny;
	ixyz.Dimension(curSize, 3);

	int jxmin = 0;
	int jxmax = nint_(shpar[0]) - 1;
	int jxmin2 = -nint_(shpar[2] + shpar[3]);
	int jxmax2 = -nint_(shpar[2]) - 1;
	int jxmin1 = -nint_(shpar[2]);
	int jxmax1 = -1;

	if (jxmin2 <= jxmax2)				// layer 2 = composition 3:
		InitSlab(jxmin2, jxmax2, 2, curSize);
	if (jxmin1 <= jxmax1)				// layer 1 = composition 2:
		InitSlab(jxmin1, jxmax1, 1, curSize);
	InitDisk(jxmin, jxmax, 0, curSize);
	ixyz.Close(nat0);
}

void Tarpbxn::InitSlab(int jxmin, int jxmax, int element, int &curSize)
{
	for(int jx=jxmin; jx<=jxmax; ++jx)
	{
		for(int jy=minJy; jy<=maxJy; ++jy)
		{
			if (Fabs(jy + x0.data[1]) < half_ * shpar[4])
			{
				for(int jz=minJz; jz<=maxJz; ++jz)
				{
					if (Fabs(jz + x0.data[2]) < half_ * shpar[5])
					{
						if (nat0 >= curSize)
						{
							ixyz.Extend(nz);
							curSize = ixyz.GetSize(0);
						}
						ixyz.Fill3(nat0, jx, jy, jz);
						int index = GetLinearAddress(nat0);
						Composer(index, element);
						iocc[index] = true;
						++nat0;
					}
				}
			}
		}
	}
}

void Tarpbxn::InitDisk(int jxmin, int jxmax, int element, int &curSize)
{
	const real r2 = quat_ * shpar[1] * shpar[1];
	for(int jx=jxmin; jx<=jxmax; ++jx)
	{
		for(int jy=minJy; jy<=maxJy; ++jy)
		{
			real y = jy + x0.data[1];
			y = y * y;
			for(int jz=minJz; jz<=maxJz; ++jz)
			{
				real z = jz + x0.data[2];
				z = z * z;
				if (y + z < r2)
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

void Tarpbxn::ShiftDipolesAndX(void)
{
	int dx = minJx - 1;
	int dy = minJy - 1;
	int dz = minJz - 1;
	if (dx || dy || dz)
	{
		for(int j=0; j<nat0; ++j)
		{
			ixyz.Value(j, 0) -= dx;
			ixyz.Value(j, 1) -= dy;
			ixyz.Value(j, 2) -= dz;
		}
	}
	x0 += Vect3<real>(dx, dy, dz);
	minJx -= dx;
	minJy -= dy;
	minJz -= dz;
	maxJx -= dx;
	maxJy -= dy;
	maxJz -= dz;
}

void Tarpbxn::OutShpar(FILE *file)
{
	fprintf(file, "slab+disk; xd,yd,xs2,xs3,ys,zs = %8.4lf %8.4lf %8.4lf %8.4lf %8.4lf %8.4lf\n", shpar[0], shpar[1], shpar[2], shpar[3], shpar[4], shpar[5]);
}

const char *TargetVerboseDescriptor_Dskblypbc(int num)
{

	return NULL;
}

REGISTER_TARGET(Dskblypbc,8,false,6,0,"disk on bilayer slab (one TUC)")
void Target_Dskblypbc::SayHello(FILE *stream)
{
	fprintf(stream, "The target %s parameters are:\n", shortDescr.c_str());
	fprintf(stream, "TUC for array of (disk on top of bilayer slab)\n");
	fprintf(stream, "With disk (comp. 1) located on x=0 surface\n");
	fprintf(stream, "disk thickness/d and diameter/d (0 0 to stop)\n");
	fprintf(stream, "%lf %lf\n", shpar[0], shpar[1]);
	fprintf(stream, "thickness/d of upper and lower layers of slab (comp. 2)\n");
	fprintf(stream, "%lf %lf\n", shpar[2], shpar[3]);
	fprintf(stream, "slab extent Ly/d  Lz/d\n");
	fprintf(stream, "%lf %lf\n", shpar[4], shpar[5]);
	fprintf(stream, "next two parameters s.b. zeros, %lf %lf\n", shpar[6], shpar[7]);
}

const char *TargetVerboseDescriptor_Dskrctngl(int num)
{
	static const char *descr[6] = 
	{
		"disk thickness/d", "disk diameter/d", "should be zero, please enter zero", "Lx/d", "Ly/d", "Lz/d"
	};
	if (num < 6)
		return descr[num];
	else
		return NULL;
}

REGISTER_TARGET(Dskrctngl,5,false,-1,0,"Disk on homogeneous rectangular block")
void Target_Dskrctngl::SayHello(FILE *stream)
{
	fprintf(stream, "The target %s parameters are:\n", shortDescr.c_str());
	fprintf(stream, "Rectangular slab with dimension Lx x Ly x Lz\n");
	fprintf(stream, "With disk located on x=0 surface\n");
	fprintf(stream, "disk thickness/d and diameter/d:\n");
	fprintf(stream, "%lf %lf\n", shpar[0], shpar[1]);
	fprintf(stream, "Lx/d  Ly/d  Lz/d:\n");
	fprintf(stream, "%lf %lf %lf\n", shpar[2], shpar[4], shpar[5]);
}

const char *TargetVerboseDescriptor_Dskrctpbc(int num)
{
	static const char *descr[7] = 
	{
		"(disk thickness)/d", "(disk diameter)/d", "(x thickness)/d", "(y length of slab)/d", "(z length of slab)/d", "(period in y direction)/d", "(period in z direction)/d"
	};
	if (num < 7)
		return descr[num];
	else
		return NULL;
}

REGISTER_TARGET(Dskrctpbc,7,false,5,0,"pillbox on slab with periodic boundary conditions")
void Target_Dskrctpbc::SayHello(FILE *stream)
{
	fprintf(stream, "The target %s parameters are:\n", shortDescr.c_str());
	fprintf(stream, "TUC for array of (disk on slab)\n");
	fprintf(stream, "With disk located on x=0 surface\n");
	fprintf(stream, "disk thickness/d and diameter/d:\n");
	fprintf(stream, "%lf %lf\n", shpar[0], shpar[1]);
	fprintf(stream, "slab dimensions Lx/d  Ly/d  Lz/d:\n");
	fprintf(stream, "%lf %lf %lf\n", shpar[2], shpar[3], shpar[4]);
	fprintf(stream, "next two parameters s.b. zeros, %lf %lf\n", shpar[5], shpar[6]);
}

void Target_Dskrctpbc::PreparePyzd()
{
	VectorA();
	pyd = shpar[6];
	pzd = shpar[7];
}