#include "StdAfx.h"

#include "Tarprsm.h"
#include "TargetManager.h"

/* **
Routine to construct triangular prism from "atoms".
triangle side lengths = a,b,c ; prism length = L
prism axis is assumed to be in x direction
normal to prism face of width a : (0,1,0)
normal to prism face of width b : (0,-cos(gamma),sin(gamma))
normal to prism face of width c : (0,-cos(beta),-sin(beta))
angles alpha,beta,gamma are opposite sides a,b,c

first layer of target array is at x=0 (target boundary at x=-0.5)
layer along side a is at int(b*sin(gamma))
     boundary is at ymax=int(b*sin(gamma))+0.5
Input:
       A      = a/d
       BA     = b/a
       CA     = c/a
       LA     = L/a
       DX(1-3)=(dx,dy,dz)/d where dx,dy,dz=lattice spacings in x,y,z
               directions, and d=(dx*dy*dz)**(1/3)=effective lattice spacing
       IOSHP  =device number for "target.out" file
              =-1 to suppress printing of "target.out"
       MXNAT  =dimensioning information (max number of atoms)
Output:
       A1(1-3)=unit vector (1,0,0) defining target axis 1 (prism axis)
       A2(1-3)=unit vector (0,1,0) defining target axis 2 (normal to
               rectangular faces with sides a,L
       NAT=number of atoms in target
       IXYZ(1-NAT,1-3)=x/d,y/d,z/d for atoms of target
       ICOMP(1-NAT,1-3)=dielectric function identifier for dipole
                        locations and 3 directions in space
       CDESCR=string describing target (up to 67 chars.)
       X0(3)=location/d in TF of dipole with IXYZ=0 0 0
             we set the TF origin to be the centroid of the prism

B.T.Draine, Princeton Univ. Obs., 2002.02.12
Fortran history records removed.

Copyright (C) 2002,2004,2006,2007,2008 B.T. Draine and P.J. Flatau
Copyright (C) 2012,2012 C++ version, Choliy V.
This code is covered by the GNU General Public License.
** */
Tarprsm::Tarprsm(TargetManager *man) : AbstractTarget(man) 
{ 
	shortDescr = string("Tarprsm");
	longDescr = string("triangular prism");
}

void Tarprsm::Sizer(void)
{
	dx = manager->CashedDx();
//
// Current version of TARPRSM is restricted to cubic lattices
	if ((dx.data[0] != onex_) || (dx.data[1] != onex_))
		Errmsg("Fatal", shortDescr.c_str(), " tarprsm does not support noncubic lattice");
//
// A = a/d
// B = b/d
// C = c/d
// L = L/d
	real a = shpar[0];
	real b = shpar[1] * shpar[0];
	real c = shpar[2] * shpar[0];
	real l = shpar[3] * shpar[0];
//
	if ((a >= b + c) || (b >= a + c) || (c >= a + b))
		Errmsg("Fatal", shortDescr.c_str(), "prism sides a,b,c do not satisfy triangle inequality:check input");
//
	real beta  = Acos((a*a + c*c - b*b) / (twox_ * a * c));
	real gamma = Acos((a*a + b*b - c*c) / (twox_ * a * b));
	cotbeta = Cos(beta) / Sin(beta);
	cotgamma = Cos(gamma) / Sin(gamma);
//
// Ideal prism extent:
//    x = xmin to x=xmax
//    triangle vertices at (0,ymax,zmin),(0,ymax,zmax),(0,ymin,za)
//    where xmin=-0.5
//          xmax=xmin+L
//          ymax=int[b*sin(gamma)]+0.5
//          ymin=ymax-b*sin(gamma)
//          zmax=a/2
//          zmin=-a/2
//          za=-a/2+c*cos(beta)
	minJx = 0; //-half_;
	maxJx = minJx + l;
	minJy = 0;
	maxJy = (int)(b * Sin(gamma)) + half_;
	maxJz = half_ * shpar[0];
	minJz = -maxJz;
//
// Now determine limits for testing x,y,z values
// Along axis (x), run from 0 to NX=INT(XMAX)
// Along y direction, run from 0 to NY=INT(YMAX)
// Along z direction, run from -NZ to NZ=INT(ZMAX)
	AllocateArrays(maxJx - minJx + 1, maxJy - minJy + 1, maxJz - minJz + 1);
}

void Tarprsm::Descriptor(void)
{
	sprintf(freeDescr, " Hexagonal prism of NAT%7d dipoles", nat0);
}

void Tarprsm::Vector(void)
{
	VectorA();
}

void Tarprsm::VectorX(void)
{
	x0.Clear();
	for(int ja=0; ja<nat0; ++ja)
	{
		x0 += Vect3<real>((real)ixyz.Value(ja, 0), (real)ixyz.Value(ja, 1), (real)ixyz.Value(ja, 2));
	}
	x0 /= -(real)nat0;
}

void Tarprsm::Allocator(void)
{
	int curSize = nx * ny;
	ixyz.Dimension(curSize, 3);

	nat0 = 0;
	for(int jx=minJx; jx<=maxJx; ++jx)
	{
		for(int jy=minJy; jy<=maxJy; ++jy)
		{
			for(int jz=minJz; jz<=maxJz; ++jz)
			{
				if ((jz <= maxJz + (jy - maxJy) * cotgamma) && (jz >= minJz - (jy - maxJy) * cotbeta))
				{
					if (nat0 >= curSize)
					{
						ixyz.Extend(2*nz + 1);
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
}

void Tarprsm::OutShpar(FILE *file)
{
	fprintf(file, "a/d = %7.4lf  b/a = %7.4lf  c/a = %7.4lf  L/a = %7.4lf\n", shpar[0], shpar[1], shpar[2], shpar[3]);
}

const char *TargetVerboseDescriptor_Trnglprsm(int num)
{
	static const char *descr[4] = 
	{
		"a/d = first triangle side/d", "b/a", "c/a", "L/a"
	};
	if (num < 4)
		return descr[num];
	else
		return NULL;
}

REGISTER_TARGET(Trnglprsm,4,false,-1,0,"Triangular prism")
void Target_Trnglprsm::SayHello(FILE *stream)
{
	fprintf(stream, "The target %s parameters are:\n", shortDescr.c_str());
	fprintf(stream, "a/d = first triangle side/d:\n");
	fprintf(stream, "%lf\n", shpar[0]);
	fprintf(stream, "b/a , c/a, L/a\n:\n");
	fprintf(stream, "%lf %lf %lf\n", shpar[1], shpar[2], shpar[3]);
}
