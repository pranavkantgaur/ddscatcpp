#include "StdAfx.h"

#include "AbstractTarget.h"
#include "TargetManager.h"

const real quat_ = (real)0.25;
const real half_ = (real)0.5;
const real zero_ = (real)0.;
const real onex_ = (real)1.;
const real twox_ = (real)2.;
const real four_ = (real)4.;

AbstractTarget::AbstractTarget(void)
{
	manager = NULL;
	Init();
}

AbstractTarget::~AbstractTarget(void)
{
	if (iocc)
		free(iocc);
	iocc = NULL;
}

AbstractTarget::AbstractTarget(TargetManager *man)
{
	manager = man;
	Init();
}

void AbstractTarget::Init(void)
{
	a1.Clear();
	a2.Clear();
	dx.Clear();
	x0.Clear();
	iocc = NULL;
	ianiso = TargetIsIsotropic;
	nat = nat0 = 0;
	shpar = manager ? manager->GetShpar() : NULL;
	minJx = minJy = minJz =  (real)ihuge_;
	maxJx = maxJy = maxJz = -(real)ihuge_;
	freeDescr[0] = '\0';
	pyd = pzd = (real)0.;
	nx = ny = nz = 0;
}

void AbstractTarget::Build(void)
{
	Sizer();
	Descriptor();
	Allocator();
	Vector();
	VectorX();
	Printer();
	ShiftDipolesAndX();
	PrepareIaniso();
	PreparePyzd();
}

void AbstractTarget::ShiftDipolesAndX(void)
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
}

void AbstractTarget::PreparePyzd()
{
	pyd = pzd = (real)0.;
}

int AbstractTarget::GetLinearAddress(int index)
{
	return ((ixyz.Value(index, 0) - minJx) * ny + (ixyz.Value(index, 1) - minJy)) * nz + (ixyz.Value(index, 2) - minJz);
}

int AbstractTarget::GetRelativeLinearAddress(int ax, int ay, int az, int mX, int mY, int mZ)
{
	return ((ax - mX) * ny + (ay - mY)) * nz + (az - mZ);
}

void AbstractTarget::Vector(void)
{
	VectorA();
}

void AbstractTarget::VectorA(void)					// Specify target axes A1 and A2
{
	a1.Set(onex_, zero_, zero_);			// A1 will be (1,0,0) in target frame
	a2.Set(zero_, onex_, zero_);			// A2 will be (0,1,0) in target frame
}

void AbstractTarget::OutVectorsAndArrays(FILE *file)
{
	a1.Fprintf(file, "%10.6lf", NULL, " = A_1 vector\n");
	a2.Fprintf(file, "%10.6lf", NULL, " = A_2 vector\n");
	dx.Fprintf(file, "%10.6lf", NULL, " = lattice spacings (d_x,d_y,d_z)/d\n");		
	x0.Fprintf(file, "%10.6lf", NULL, " = lattice offset x0(1-3) = (x_TF,y_TF,z_TF)/d ");		 
	fprintf(file, "for dipole 0 0 0\n     JA  IX  IY  IZ ICOMP(x,y,z)\n");
	for(int ja=0; ja<nat0; ++ja)
	{
		fprintf(file, "%7d", ja);
		ixyz.PrintLine(file, "%4d", ja);
		int index = GetLinearAddress(ja);
		icomp.PrintLine(file, "%2d", index);
		fprintf(file, " %d\n", index);
	}
}

void AbstractTarget::OutShpar(FILE *file)
{
	fprintf(file, "OutShpar for the target %s is not implemented yet.\n", shortDescr.c_str());
}

real AbstractTarget::DipoleScalar(const Vect3<real> &op, int absIndex)
{
	int *xyz = ixyz.Row(absIndex);
	return ((x0 + Vect3<real>((real)xyz[0], (real)xyz[1], (real)xyz[2])).Magnify(dx).Scalar(op));
}

void AbstractTarget::DipoleVector(int absIndex, real *res)
{
	int *xyz = ixyz.Row(absIndex);
	memcpy(res, (x0 + Vect3<real>((real)xyz[0], (real)xyz[1], (real)xyz[2])).Magnify(dx).data, 3*sizeof(real));
}

void AbstractTarget::Printer(void)					// Technical output
{
	if (manager->Ioshp())
	{
		FILE *ioshpFile = fopen("target.out", "w");
		if (ioshpFile)
		{
			fprintf(ioshpFile, " >%s %s;", shortDescr.c_str(), longDescr.c_str());
			OutShpar(ioshpFile);
			fprintf(ioshpFile, "%10d = NAT\n", nat0);
			OutVectorsAndArrays(ioshpFile);
			fclose(ioshpFile);
		}
		else
			fprintf(stderr, "Cannot open ioshpFile in AbstractFile::Printer for %s.", shortDescr.c_str());
	}
}

void AbstractTarget::Composer(int index, int item)
{
	IsotropicComposer(index, item);
}

void AbstractTarget::IsotropicComposer(int index, int item)
{
	icomp.Fill3(index, (short)(item + 1));
}

void AbstractTarget::AnisotropicComposer(int index, int item)
{
	int dx = 3*item;
	icomp.Fill3(index, (short)(dx + 1), (short)(dx + 2), (short)(dx + 3));
}

void AbstractTarget::Debug(const char *fileName)
{
	FILE *fx = fopen(fileName, "w");
	if (!fx)
	{
		Wrimsg("AbstractTarget::Debug", "Cannot open Debug.txt file fx\n");
		return;
	}
	fprintf(fx, "Nat0 = %d, nat = %d, nx,ny,nz = %d %d %d\n", nat0, nat, nx, ny, nz);
	fprintf(fx, "Minmax = %d %d %d %d %d %d\n", minJx, maxJx, minJy, maxJy, minJz, maxJz);
	x0.Fprintf(fx, "%lf");
	fprintf(fx, "\n");
    for(int ii=0; ii<nat0; ++ii)
	{
		fprintf(fx, "%5d", ii);
		ixyz.PrintLine(fx, "%5d", ii);
		int jj = GetLinearAddress(ii);
		fprintf(fx, "%5d", iocc[jj]);
		icomp.PrintLine(fx, "%5d", jj);
		fprintf(fx, " %d\n", jj);
	}
	fclose(fx);
}

void AbstractTarget::AllocateArrays(int sizeX, int sizeY, int sizeZ, bool bWithIxyz)
{
	if ((sizeX != nx) || (sizeY != ny) || (sizeZ != nz))
	{
		int del = Analyze(sizeX);
		if (del > 0) maxJx += del;
		nx = sizeX;
		del = Analyze(sizeY);
		if (del > 0) maxJy += del;
		ny = sizeY;
		del = Analyze(sizeZ);
		if (del > 0) maxJz += del;
		nz = sizeZ;
		nat = nx * ny * nz;
		iocc = (bool *)realloc(iocc, nat * sizeof(bool));
		memset(iocc, 0, nat * sizeof(bool));
		icomp.Deallocate();
		icomp.Dimension(nat, 3);
		icomp.Clear();
		if (bWithIxyz)
		{
			ixyz.Deallocate();
			ixyz.Dimension(nat, 3);
			ixyz.Clear();
		}
	}
}

Vect6<int> AbstractTarget::IxyzMinmax(void)
{
	Vect6<int>iMinmax;
	iMinmax.Set(ihuge_, -ihuge_, ihuge_, -ihuge_, ihuge_, -ihuge_);
	for(int i=0; i<nat0; ++i)
	{
		if (ixyz.Value(i, 0) < iMinmax.data[0]) iMinmax.data[0] = ixyz.Value(i, 0);
		if (ixyz.Value(i, 0) > iMinmax.data[1]) iMinmax.data[1] = ixyz.Value(i, 0);
		if (ixyz.Value(i, 1) < iMinmax.data[2]) iMinmax.data[2] = ixyz.Value(i, 1);
		if (ixyz.Value(i, 1) > iMinmax.data[3]) iMinmax.data[3] = ixyz.Value(i, 1);
		if (ixyz.Value(i, 2) < iMinmax.data[4]) iMinmax.data[4] = ixyz.Value(i, 2);
		if (ixyz.Value(i, 2) > iMinmax.data[5]) iMinmax.data[5] = ixyz.Value(i, 2);
	}
	return iMinmax;
}

void AbstractTarget::PrepareIaniso(void)
{
	ianiso = TargetIsIsotropic;
}

void AbstractTarget::PrepareIanisoSpecial(void)
{
	ianiso = TargetIsAnisotropic;
}

void AbstractTarget::SayHello(void)
{
	fprintf(stdout, "The target %s has no SayHello yet.\n", shortDescr.c_str());
}

int AbstractTarget::Analyze(int &number)
{
	int num = GetNearestFactor(number);
	if (num != number)
	{
		char Buffer[256];
		sprintf(Buffer, "Dimension %d changed to %d", number, num);
		Wrimsg("AbstractTarget", Buffer);
		int delta = num - number;
		number = num;
		return delta;
	}
	return 0;
}

Vect3<real> AbstractTarget::Prinaxis(void)
{
/* **
Subroutine to determine principal axes of moment of inertia tensor of target.
Target is assumed to be homogeneous.

Given:
MXNAT       =limit on largest allowed value of NAT=number of dipoles
NAT         =number of dipoles in target
IXYZ(J,1-3) =x,y,z location of dipole Jon lattice
DX(1-3)     =(dx/d,dy/d,dz/d) where (dx,dy,dz)=lattice spacing in x,y,z directions, and d = (dx*dy*dz)**(1/3) = effective lattice spacing
ICOMP(J,1-3)=composition for dipole J; x,y,z directions

Returns:
A1(1-3)     =principal axis A1 in target frame (normalized)
A2(1-3)     =principal axis A2 in target frame (normalized)
EIGV(1-3)   =(eigenvalues of moment of inertia tensor)/(0.4*M*aeff^2)
             definition of aeff = (3*V/4*pi)^{1/3} , V = solid volume
             EIGV are in decreasing order:
             EIGV(1) is largest, EIGV(3) is smallest

             A1 = principal axis with largest moment of inertia
             A2 = principal axis with second-largest moment of inertia

B.T. Draine, Princeton Univ. Observatory, 96.01.26
History:
Fortran versions history removed.
Copyright (C) 1993,1994,1995,1996,1997,1998,2000,2004,2007
              B.T. Draine and P.J. Flatau
Copyright (C) 2012,2013, C++ versions, Choliy V.
This code is covered by the GNU General Public License.
** */
//
// Compute moment of inertia tensor.
	int ja, jd, j, k;
	Vect3<real> cm;
	for(ja=0; ja<nat0; ++ja)
		for(jd=0; jd<3; ++jd)
			cm.data[jd] += (ixyz.Value(ja, jd) * dx.data[jd]);
	cm /= (real)nat0;
//
// Compute moment of inertia tensor
	Matc3<real> a, q;
	real dsum = zero_;
	for(ja=0; ja<nat0; ++ja)
	{
		for(k=0; k<3; ++k)
		{
			real temp = (ixyz.Value(ja, k) - cm.data[k]) * dx.data[k];
			dsum += (temp * temp);
		}
	}

	for(j=0; j<3; ++j)
		a.Data(j, j) = (real)(dsum + (nat0/(real)18.) * dx.ModSquared());
	for(j=0; j<3; ++j)
	{
		for(k=0; k<3; ++k)
		{
			real temp = dx.data[j] * dx.data[k];
			dsum = zero_;
			for(ja=0; ja<nat0; ++ja)
				dsum += (ixyz.Value(ja, j) - cm.data[j]) * (ixyz.Value(ja, k) - cm.data[k]) * temp;
			a.Data(j, k) -= dsum;
		}
	}
//
// Normalize moment of inertia tensor to units of 0.4*mass*a_eff**2
	real denom = (real)0.4 * nat0 * Pow((real)(0.75 * nat0 / Pi), (real)(2./3.));
	a /= denom;
//
// Compute eigenvalues and eigenvectors of A(J,K)
	Vect3<real> w;
	Dsyevj3(a, q, w);
//
// On return from DSYEVJ3:
//     W(i) = eigenvalues i=1,2,3
//     Q(i,j=1,3) = eigenvector corresponding to eigenvalue W(i) 
//
// Normalize the eigenvectors:
	for(j=0; j<3; ++j)
	{
		real znorm2 = zero_;
		for(k=0; k<3; ++k)
			znorm2 += q.Data(k,j) * q.Data(k,j);
//
// Make first component of each eigenvector positive
		real znorm = Sqrt(znorm2);
		if (q.Data(0, j) < zero_)
			znorm = -znorm;
		else
		{
			if (q.Data(0, j) == zero_)
				if(q.Data(1, j) < zero_)
					znorm = -znorm;
		}
		for(k=0; k<3; ++k)
			q.Data(k, j) /= znorm;
	}
//
// Order the eigenvectors by eigenvalue, largest first:
	int j1 = 0, j2 = 0;
	if (w.data[1] > w.data[0])
		j1 = 1;
	if (w.data[2] > w.data[j1])
		j1 = 2;
	switch(j1)
	{
	case 0:
		j2 = 1;
		if (w.data[2] > w.data[1]) 
			j2 = 2;
		break;

	case 1:
		j2 = 0;
		if (w.data[2] > w.data[0])
			j2 = 2;
		break;

	case 2:
		j2 = 0;
		if (w.data[1] > w.data[0])
			j2 = 1;
		break;
	}

	for(k=0; k<3; ++k)
	{
		a1.data[k] = q.Data(k, j1);
		a2.data[k] = q.Data(k, j2);
	}
//
// Set J3=3 if J1+J2=1+2=3
//        2         =1+3=4
//        1         =2+3=5
	int j3 = 3 - j1 - j2;
//
// Having fixed eigenvectors a_1 and a_2, compute a_3 = a_1 x a_2
	Vect3<real> a3 = a1 * a2;
//
	char Buffer[256];
	sprintf(Buffer, "(I_1/.4Ma_eff^2 =%8.5lf)", w.data[j1]);
	Wrimsg("Prinaxis", Buffer);
	sprintf(Buffer, "(I_2/.4Ma_eff^2 =%8.5lf)", w.data[j2]);
	Wrimsg("Prinaxis", Buffer);
	sprintf(Buffer, "(I_3/.4Ma_eff^2 =%8.5lf)", w.data[j3]);
	Wrimsg("Prinaxis", Buffer);
	a1.Sprintf(Buffer, "%9.5lf", "(     axis a_1 = (", ")");
	Wrimsg("Prinaxis", Buffer);
	a2.Sprintf(Buffer, "%9.5lf", "(     axis a_2 = (", ")");
	Wrimsg("Prinaxis", Buffer);
	a3.Sprintf(Buffer, "%9.5lf", "(     axis a_3 = (", ")");
	Wrimsg("Prinaxis", Buffer);

	return w;
}

void AbstractTarget::Dsyevj3(Matc3<real> &a, Matc3<real> &q, Vect3<real> &w)
{
/* ** 
----------------------------------------------------------------------
Numerical diagonalization of 3x3 matrcies
Copyright (C) 2006  Joachim Kopp
----------------------------------------------------------------------
This library is free software; you can redistribute it and/or
modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation; either
version 2.1 of the License, or (at your option) any later version.

This library is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public
License along with this library; if not, write to the Free Software
Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-13
----------------------------------------------------------------------

----------------------------------------------------------------------
Calculates the eigenvalues and normalized eigenvectors of a symmetric
matrix A using the Jacobi algorithm.
The upper triangular part of A is destroyed during the calculation,
the diagonal elements are read but not destroyed, and the lower
triangular elements are not referenced at all.
----------------------------------------------------------------------
Parameters:
  A: The symmetric input matrix
  Q: Storage buffer for eigenvectors
  W: Storage buffer for eigenvalues
----------------------------------------------------------------------
C++ version was created by Choliy V., 2012.
** */
	const int n = 3;
//
// Initialize Q to the identitity matrix 
// --- This loop can be omitted if only the eigenvalues are desired -
	q.Unity();
//
// Initialize W to diag(A)
	int x, y;
	for(x=0; x<n; ++x)
		w.data[x] = a.Data(x, x);
//
// Calculate SQR(tr(A))
	real sd = zero_;
	for(x=0; x<n; ++x)
		sd += Fabs(w.data[x]);
	sd = sd*sd;
//
// Main iteration loop
	int i;
	for(i=0; i<50; ++i)
	{
//
// Test for convergence
		real so = zero_;
		for(x=0; x<n; ++x)
			for(y=x+1; y<n; ++y)
				so += Fabs(a.Data(x, y));
		if (so == zero_)
		{
			fprintf(stderr, "Dsyevj3 converged.\n");
			return;
		}

		real thresh = (i < 3) ? (real)0.2*so/(n*n) : zero_;
//        
// Do sweep
		real t, theta;
		for(x=0; x<n; ++x)
		{
			for(y=x+1; y<n; ++y)
			{
				real g = (real)100.0*(Fabs(a.Data(x,y)));
				if ((i > 3) && (Fabs(w.data[x])+g == Fabs(w.data[x])) && (Fabs(w.data[y])+g == Fabs(w.data[y])))
					a.Data(x,y) = zero_;
				else
				{
					if (Fabs(a.Data(x,y)) > thresh)							// Calculate Jacobi transformation
					{
						real h = w.data[y] - w.data[x];
						if (Fabs(h)+g == Fabs(h))
							t = a.Data(x,y) / h;
						else
						{
							theta = half_*h/a.Data(x,y);
							t = (theta < zero_) ? -(real)1./(Sqrt((real)1. + theta*theta) - theta) : (real)1./(sqrt((real)1. + theta*theta) + theta);
						}
						real c = (real)1./Sqrt((real)1. + t*t);
						real s = t*c;
						real z = t*a.Data(x,y);
// Apply Jacobi transformation
						a.Data(x,y) = zero_;
						w.data[x] -= z;
						w.data[y] += z;
						int r;
						for(r=0; r<x; ++r)
						{
							t = a.Data(r, x);
							a.Data(r, x) = c*t - s*a.Data(r, y);
							a.Data(r, y) = s*t + c*a.Data(r, y);
						}
						for(r=x+1; r<y; ++r)
						{
							t = a.Data(x, r);
							a.Data(x, r) = c*t - s*a.Data(r, y);
							a.Data(r, y) = s*t + c*a.Data(r, y);
						}
						for(r=y+1; r<n; ++r)
						{
							t = a.Data(x, r);
							a.Data(x, r) = c*t - s*a.Data(y, r);
							a.Data(y, r) = s*t + c*a.Data(y, r);
						}
// Update eigenvectors
// --- This loop can be omitted if only the eigenvalues are desired
						for(r=0; r<n; ++r)
						{
							t = q.Data(r,x);
							q.Data(r,x) = c*t - s*q.Data(r,y);
							q.Data(r,y) = s*t + c*q.Data(r,y);
						}
					}
				}
			}
		}
	}
	fprintf(stderr, ">DSYEVJ3: No convergence.\n");
}

void AbstractTarget::Evale(Vect3<Complex> &cxe00, Vect3<real> &akd, Complex *cxe)
{
/* **
Given:   CXE00(1-3)=Incident E field at origin (Complex) at t=0
         AKD(1-3)=(kx,ky,kz)*d for incident wave (d=effective lattice spacing)
         DX(1-3)=(dx/d,dy/d,dz/d) for lattice (dx,dy,dz=lattice spacings in x,y,z directions, d=(dx*dy*dz)**(1/3)
         X0(1-3)=(x,y,z)location/(d*DX(1-3)) in TF of lattice site with IX=0,IY=0,IZ=0
         IXYZ0(1-NAT0,3)=[x-x0(1)]/dx, [y-x0(2)]/dy, [z-x0(3)]/dz for each of NAT0 physical dipoles
         MXNAT,MXN3=dimensioning information
         NAT0=number of dipoles in physical target
         NAT=number of locations at which to calculate CXE

Returns: CXE(1-NAT,3)=incident E field at NAT locations at t=0

B.T.Draine, Princeton Univ. Obs., 88.05.09
History:
Fortran history records removed.
Copyright (C) 1993,1997,2007 B.T. Draine and P.J. Flatau
Copyright (C) 2012, C++ version, V.Choliy
This code is covered by the GNU General Public License.
** */
//
// Evaluate electric field vector at each dipole location.
// If NAT=NAT0, then evaluate E only at occupied sites.
// If NAT>NAT0, then evaluate E at all sites.

	real *akdData = akd.data;
	real *dxData = dx.data;
	real *x0Data = x0.data;

	int ia, m, index;
	if (nat == nat0)
	{
		int *ixyzData = ixyz.GetData();
		for(ia=0; ia<nat0; ++ia)
		{
			index = 3*ia;
			real x = zero_;
			for(m=0; m<3; ++m)
				x += (akdData[m] * dxData[m] * ((real)ixyzData[index+m] + x0Data[m]));
			Complex cxfac = Complex(Cos(x), Sin(x));
			for(m=0; m<3; ++m)
				cxe[index + m] = cxe00.data[m] * cxfac;
		}
	}
	else
	{
		ia = 0;
		for(int ix=1; ix<=nx; ++ix)
		{
			real x1 = akdData[0] * dxData[0] * (x0Data[0] + ix);
			for(int iy=1; iy<=ny; ++iy)
			{
				real x2 = x1 + akdData[1] * dxData[1] * (x0Data[1] + iy);
				for(int iz=1; iz<=nz; ++iz)
				{
					real x = x2 + akdData[2] * dxData[2] * (x0Data[2] + iz);
					Complex cxfac = Complex(Cos(x), Sin(x));
					for(m=0; m<3; ++m)
						cxe[ia + m] = cxe00.data[m] * cxfac;
					ia += 3;
				}
			}
		}
	}
}

void AbstractTarget::Reduce(Complex *cxv)
{
/* **
Given:
     CXV(1-3*NAT) defined for NAT lattice sites with ordering
         (v_x1,...,v_xj,...,v_xNAT,v_y1,...,v_yj,...,v_yNAT,
          v_z1,...,v_zj,...,v_zNAT)
         where index j runs over occupied and unoccupied lattice sites
     IOCC(1-NAT) = 0 or 1 depending on whether lattice site is vacant or occupied
     MXN3,MXNAT = dimensioning information
     NAT = number of lattice sites
     NAT0 = number of occupied lattice sites

Returns
     CXV(1-3*NAT0) defined for NAT0 occupied lattice sites with ordering
        (v_x1,...,v_xj,...,v_xNAT0,v_y1,...,v_yj,...,v_yNAT0, v_z1,...,v_zj,...,v_zNAT0)
        where index j now runs over only occupied lattice sites

B.T.Draine, Princeton Univ. Observatory, 90.11.29
Fortran history records removed.

Copyright (C) 1993, 2006 B.T. Draine and P.J. Flatau
Copyright (C) 2012, C++ versions, Choliy V.
This code is covered by the GNU General Public License.
** */
	if (nat == nat0)
		return;

	int joc = 0;
	for(int j=0; j<nat; ++j)
	{
		if (iocc[j] == true)
		{
			memcpy(cxv+joc, cxv+3*j, 3*sizeof(Complex));
			joc += 3;
		}
	}
}

void AbstractTarget::Unreduce(Complex *cxv)
{
/* **
Given:
   CXV(1:MXN3)   = complex vector in "reduced" form
   IOCC(1:MXNAT) = 0/1 if site is vacant/occupied
   MXN3          = dimensioning information
   MXNAT         = dimensioning information
   NAT           = number of lattice sites
   NAT0          = number of occupied lattice sites

Returns:
   CXV(1:MXN3)   = complex vector in "natural" form

History records of Fortran versions removed.
Copyright (C) 2011, B.T. Draine and P.J. Flatau
Copyright (C) 2012,2013 C++ versions, Choliy V.
This code is covered by the GNU General Public License.
** */

	if (nat <= nat0)
		return;

	memset(cxv + 3*nat0, 0, 3*(nat - nat0)*sizeof(Complex));

	int jocc = 3*(nat0 - 1);
	for(int j=nat-1; j>=0; j--)
	{
		if(iocc[j] == true)
		{
			if (3*j != jocc)
			{
				memcpy(cxv+3*j, cxv+jocc, 3*sizeof(Complex));
				memset(cxv+jocc, 0, 3*sizeof(Complex));
			}
			jocc -= 3;
		}
	}
}

void AbstractTarget::InternalMinMax(int x, int y, int z)
{
	if (minJx > x) minJx = x;
	if (maxJx < x) maxJx = x;
	if (minJy > y) minJy = y;
	if (maxJy < y) maxJy = y;
	if (minJz > z) minJz = z;
	if (maxJz < z) maxJz = z;
}

void AbstractTarget::SetMin(const Vect6<int> &op)
{
	minJx = op.data[0];
	minJy = op.data[2];
	minJz = op.data[4];
}

bool AbstractTarget::IsDipolePhysicallyAnisotropic(unsigned int index)
{
	if ((icomp.Value(index, 0) != icomp.Value(index, 1)) || (icomp.Value(index, 0) != icomp.Value(index, 2)))
		return true;
	else
		return false;
}