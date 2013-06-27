#include "StdAfx.h"

#include "Tarrctblk3.h"
#include "TargetManager.h"

/* **
Purpose: to construct target consisting of 3 stacked rectangular blocks: boundary between blocks 1 and 2 lies in x=0 plane
	block 1, of composition 1, extends from x/d= 0 to XS1 y/d=-YS1/2 to +YS1/2 z/d=-ZS1/2 to +ZS1/2
	block 2, of composition 2, extends from x/d=-XS2 to 0 y/d=-YS2/2 to +YS2/2 z/d=-ZS2/2 to +ZS2/2
	block 3, of composition 3, extends from x/d=-(XS2+XS3) to -XS2 y/d=-YS3/2 to +YS3/2 z/d=-ZS3/2 to +ZS3/2

Input:
       XS1=(x thickness of composition 1 of slab)/d   (d=lattice spacing)
       YS1=(y extent of composition 1 of slab)/d   (d=lattice spacing)
       ZS1=(z extent of composition 1 of slab)/d   (d=lattice spacing)
       XS2=(x thickness of composition 2 of slab)/d
       YS2=(y extent of composition 2 of slab)/d
       ZS2=(z extent of composition 2 of slab)/d
       XS3=(x thickness of composition 3 of slab)/d
       YS3=(y extent of composition 3 of slab)/d
       ZS3=(z extent of composition 3 of slab)/d
       DX(1-3)=(dx,dy,dz)/d where dx,dy,dz=lattice spacings in x,y,z directions, and d=(dx*dy*dz)**(1/3)=effective lattice spacing
       IOSHP=device number for "target.out" file =-1 to suppress printing of "target.out"
       MXNAT=dimensioning information (max number of atoms)

Output:
       A1(1-3)=unit vector (1,0,0) defining target axis 1 in Target Frame
       A2(1-3)=unit vector (0,1,0) defining target axis 2 in Target Frame
       NAT=number of atoms in target
       IXYZ(1-NAT,1-3)=x/d,y/d,z/d for atoms of target
       ICOMP(1-NAT,1-3)=composition number for locations 1-NAT and directions 1-3
       X0(1-3)=offset vector, such that
               x/d=(IXYZ(J,1)+X0(1))*DX(1)
               y/d=(IXYZ(J,2)+X0(2))*DX(2)
               z/d=(IXYZ(J,3)+X0(3))*DX(3)
               where (x,y,z)=location of dipole in Target Frame
               with (0,0,0)=point where axis of disk intersects upper surface of slab

Target consists of three stacked rectangular slabs with x-axis running
       through center of each slab
       slab 1 is of extent XS1 x YS1 x ZS1 (in units of d)
       slab 2 is of extent XS2 x YS2 x ZS2 (in units of d)
       slab 3 is of extent XS3 x YS3 x ZS3 (in units of d)

       center of upper surface of slab 1 is at (XS1*d*DX(1), 0 , 0)
       center of upper surface of slab 2 is at (0,0,0)
       center of lower surface of slab 2 = upper surface of slab 3 is at (-XS2*d*DX(2) , 0 , 0 )
       boundary between compositions 1 and 2 is in the x=0 plane
       boundary between compositions 2 and 3 is in the x=-XS2*d*DX(1) plane
       lower surface of slab 3 is in the x=-(XS2+XS3)*d*DX(1) plane

       if max(YS1,YS2,YS3)/DX(2) is even: lattice will be located at y/(d*DX(2)) = +/-0.5, +/-1.5, . odd: y/(d*DX(2)) = 0, +/-1, +/-2, ..
       if max(ZS1,ZS2,ZS3)/DX(3) is even: lattice will be located at z/(d*DX(3)) = +/-0.5, +/-1.5, . odd: z/(d*DX(3)) = 0, +/-1, +/-2, ..

Fortran history records removed.

Copyright (C) 2008,2010 B.T. Draine and P.J. Flatau
Copyright (c) 2012, C++ version, Choliy V.
This code is covered by the GNU General Public License.
** */

Tarrctblk3::Tarrctblk3(TargetManager *man) : AbstractTarget(man) 
{ 
	shortDescr = string("Tarrctblk3");
	longDescr = string("Layered block structure");
}

void Tarrctblk3::Sizer(void)
{
	dx = manager->CashedDx();
// 
// Target reference point =intersection of x_TF axis with plane separating slabs 1 and 2
// Set lattice offset in x direction: Dipole with IX=0 is located at x/(d*DX(1))=0.5
	minJx = -nint_((shpar[3] + shpar[6]) / dx.data[0]);
	maxJx =  nint_(shpar[0] / dx.data[0]) - 1;
	x0.data[0] = half_;
//
// Determine lattice offset in y and z directions
//   if NY is even, dipole with IY=0 is located at y/d=0.5 odd  0
//   if NZ is even, dipole with IZ=0 is located at z/d=0.5 odd  0
	real term = max_(shpar[1], shpar[4]);
	term = max_(term, shpar[7]);
	int ny1 = nint_(term / dx.data[1]);
	x0.data[1] = ny1 % 2 ? zero_ : half_;
	minJy = -ny1*half_ - x0.data[1] + half_;
	maxJy =  ny1*half_ - x0.data[1] - half_;
//
	term = max_(shpar[2], shpar[5]);
	term = max_(term, shpar[8]);
	int nz1 = nint_(term / dx.data[2]);
	x0.data[2] = nz1 % 2 ? zero_ : half_;
	minJz = -nz1*half_ - x0.data[2] + half_;
	maxJz =  nz1*half_ - x0.data[2] - half_;
	AllocateArrays(maxJx - minJx + 1, maxJy - minJy + 1, maxJz - minJz + 1);
}

void Tarrctblk3::Descriptor(void)
{
	if (shpar[0] * shpar[3] * shpar[6] > zero_)
		strcpy(freeDescr, " Three stacked rectangular slabs");
	else
		strcpy(freeDescr, " Two stacked rectangular slabs");
}

void Tarrctblk3::Vector(void)
{
	VectorA();
}

//
// Creates elementary block
int Tarrctblk3::ElementaryBlock(int jxmin, int jxmax, int jymin, int jymax, int jzmin, int jzmax, int comp, int &curSize)
{
	for(int jx=jxmin; jx<=jxmax; ++jx)
	{
		for(int jy=jymin; jy<=jymax; ++jy)
		{
			for(int jz=jzmin; jz<=jzmax; ++jz)
			{
				if (nat0 >= curSize)
				{
					ixyz.Extend(nz);
					curSize = ixyz.GetSize(0);
				}
				ixyz.Fill3(nat0, jx, jy, jz);
				int index = GetLinearAddress(nat0);
				Composer(index, comp);
				iocc[index] = true;
				++nat0;
			}
		}
	}
	return nat0;
}

//
// ChB: from negative to positive x
void Tarrctblk3::Allocator(void)
{
//	int jymin, jymax, jzmin, jzmax;
	int nats1 = 0, nats2 = 0, nats3 = 0;
	int curSize  = nx * ny;
	ixyz.Dimension(curSize, 3);
	nat0 = 0;
//
// Lower slab layer = composition 3
	int jxmin1 = 0;
	int jxmax1 = nint_(shpar[0] / dx.data[0]) - 1;
	int jxmin2 = -nint_(shpar[3] / dx.data[0]);
	int jxmax2 = -1;
	int jxmin3 = -nint_((shpar[3] + shpar[6]) / dx.data[0]);
	int jxmax3 = -nint_(shpar[3] / dx.data[0]) - 1;
//
//	jymin = nint_(-half_ * shpar[7] / dx.data[1] - x0.data[1] + half_);
//	jymax = nint_( half_ * shpar[7] / dx.data[1] - x0.data[1] - half_);
//	jzmin = nint_(-half_ * shpar[8] / dx.data[2] - x0.data[2] + half_);
//	jzmax = nint_( half_ * shpar[8] / dx.data[2] - x0.data[2] - half_);
	if (jxmin3 < jxmin2)
		nats3 = ElementaryBlock(jxmin3, jxmax3, 
			nint_(-half_ * shpar[7] / dx.data[1] - x0.data[1] + half_), nint_( half_ * shpar[7] / dx.data[1] - x0.data[1] - half_),
			nint_(-half_ * shpar[8] / dx.data[2] - x0.data[2] + half_), nint_( half_ * shpar[8] / dx.data[2] - x0.data[2] - half_), 2, curSize);
//
// Now do block 2 (material 2)
// layer 1: extends from 0 to  xs1 [composition 2] dipoles from 1 to NAT1
// layer 2: extends from -(xs1+xs2) to -xs1 [composition 3] dipoles from   JXMIN2   to (JXMIN1-1)
//	jymin = nint_(-half_ * shpar[4] / dx.data[1] - x0.data[1] + half_);
//	jymax = nint_( half_ * shpar[4] / dx.data[1] - x0.data[1] - half_);
//	jzmin = nint_(-half_ * shpar[5] / dx.data[2] - x0.data[2] + half_);
//	jzmax = nint_( half_ * shpar[5] / dx.data[2] - x0.data[2] - half_);
//
// Middle layer = composition 2:
	if (jxmin2 < jxmin1)
		nats2 = ElementaryBlock(jxmin2, jxmax2, 
			nint_(-half_ * shpar[4] / dx.data[1] - x0.data[1] + half_), nint_( half_ * shpar[4] / dx.data[1] - x0.data[1] - half_),
			nint_(-half_ * shpar[5] / dx.data[2] - x0.data[2] + half_), nint_( half_ * shpar[5] / dx.data[2] - x0.data[2] - half_), 1, curSize);
//
// First do the top block (material 1)
//	jymin = nint_(-half_ * shpar[1] / dx.data[1] - x0.data[1] + half_);
//	jymax = nint_( half_ * shpar[1] / dx.data[1] - x0.data[1] - half_);
//	jzmin = nint_(-half_ * shpar[2] / dx.data[2] - x0.data[2] + half_);
//	jzmax = nint_( half_ * shpar[2] / dx.data[2] - x0.data[2] - half_);
//
// number of dipoles in y direction 
//     NY1 = jymax1-jymin1+1 
//         = nint(0.5*ys1/dx(2)-x0(2)+delta)-nint(-0.5*ys1/dx(2)-x0(2)-delta)+1
//         = nint(0.5*ys1/dx(2)-x0(2)+delta)+nint(0.5*ys1/dx(2)+x0(2)+delta)+1
// number of dipoles in Z direction 
//     NZ1 = jzmax1-jzmin1+1 
//         = nint(0.5*zs1/dx(3)-x0(3)+delta)-nint(-0.5*zs1/dx(3)-x0(3)-delta)+1
//         = nint(0.5*zs1/dx(3)-x0(3)+delta)+nint(0.5*zs1/dx(3)+x0(3)+delta)+1
	if (jxmin1 <= jxmax1)
		nats1 = ElementaryBlock(jxmin1, jxmax1, 
			nint_(-half_ * shpar[1] / dx.data[1] - x0.data[1] + half_), nint_( half_ * shpar[1] / dx.data[1] - x0.data[1] - half_), 
			nint_(-half_ * shpar[2] / dx.data[2] - x0.data[2] + half_), nint_( half_ * shpar[2] / dx.data[2] - x0.data[2] - half_), 0, curSize);
//
//	nats2 -= nats1;
//	nats3 -= (nats1 + nats2);
//
	if (nats1 * nats2 * nats3 > 0)
		sprintf(freeDescr, " Trilayer block structure with %7d %7d %7d = %7d dipoles", nats1, nats2, nats3, nat0);
	else
		sprintf(freeDescr, " Bilayer block structure with %7d %7d %7d = %7d dipoles", nats1, nats2, nats3, nat0);
}

void Tarrctblk3::OutShpar(FILE *file)
{
	fprintf(file, "XS1,YS1,ZS1,XS2,YS2,ZS2,XS3,YS3,ZS3=%8.4lf%8.4lf%8.4lf%8.4lf%8.4lf%8.4lf%8.4lf%8.4lf%8.4lf\n",
		shpar[0], shpar[1], shpar[2], shpar[3], shpar[4], shpar[5], shpar[6], shpar[7], shpar[8]);
}

const char *TargetVerboseDescriptor_Bislinpbc(int num)
{
	static const char *descr[6] = 
	{
		"x-thickness/d of line", "y-width/d of line", "x-thickness/d of upper layer", "x-thickness/d of lower layer", "L_y/d dimension of slab", NULL
	};
	if (num < 6)
		return descr[num];
	else
		return NULL;
}

REGISTER_TARGET(Bislinpbc,6,false,5,0,"bilayer slab with parallel lines (one TUC)")
void Target_Bislinpbc::SayHello(FILE *stream)
{
	fprintf(stream, "TUC for bilayer slab with parallel line on top\n");
	fprintf(stream, "The target %s parameters are:\n", shortDescr.c_str());
	fprintf(stream, "x-thickness/d of line %20.16lf\n", shpar[0]);
	fprintf(stream, "y-width/d of line %20.16lf\n", shpar[1]);
	fprintf(stream, "x-thickness/d of upper layer %20.16lf\n", shpar[2]);
	fprintf(stream, "x-thickness/d of lower layer %20.16lf\n", shpar[3]);
	fprintf(stream, "L_y/d dimension of slab %20.16lf\n", shpar[4]);
	fprintf(stream, "Paramater 5 should be zero, %20.16lf\n", shpar[5]);
}

const char *TargetVerboseDescriptor_Rctglblk3(int num)
{
	static const char *descr[9] = 
	{
		"x length/d for first block", "y length/d for first block", "z length/d for first block", 
		"x length/d for second block", "y length/d for second block", "z length/d for second block", 
		"x length/d for third block", "y length/d for third block", "z length/d for third block"
	};
	if (num < 9)
		return descr[num];
	else
		return NULL;
}

REGISTER_TARGET(Rctglblk3,9,false,-1,0,"Stack of 3 rectangular blocks")
void Target_Rctglblk3::SayHello(FILE *stream)
{
	fprintf(stream, "The target %s parameters are:\n", shortDescr.c_str());
	fprintf(stream, "x,y,z length/d for first block:\n");
	fprintf(stream, "%lf %lf %lf\n", shpar[0], shpar[1], shpar[2]);
	fprintf(stream, "x,y,z length/d for second block:\n");
	fprintf(stream, "%lf %lf %lf\n", shpar[3], shpar[4], shpar[5]);
	fprintf(stream, "x,y,z length/d for third block:\n");
	fprintf(stream, "%lf %lf %lf\n", shpar[6], shpar[7], shpar[8]);
}

const char *TargetVerboseDescriptor_Trilyrpbc(int num)
{

	return NULL;
}

REGISTER_TARGET(Trilyrpbc,11,false,9,0,"3 stacked rect. blocks (one TUC)")
void Target_Trilyrpbc::SayHello(FILE *stream)
{
	fprintf(stream, "TRILYRPBC TUC=three stacked rectangular blocks\n");
	fprintf(stream, "The target %s parameters are:\n", shortDescr.c_str());
	fprintf(stream, "size x/d, y/d, z/d for upper block:\n");
	fprintf(stream, "%20.16lf %20.16lf %20.16lf\n", shpar[0], shpar[1], shpar[2]);
	fprintf(stream, "size x/d, y/d, z/d for middle block:\n");
	fprintf(stream, "%20.16lf %20.16lf %20.16lf\n", shpar[3], shpar[4], shpar[5]);
	fprintf(stream, "size x/d, y/d, z/d for bottom block:\n");
	fprintf(stream, "%20.16lf %20.16lf %20.16lf\n", shpar[6], shpar[7], shpar[8]);
	fprintf(stream, "next two parameters s.b. zeros, %20.16lf %20.16lf\n", shpar[9], shpar[10]);
}
