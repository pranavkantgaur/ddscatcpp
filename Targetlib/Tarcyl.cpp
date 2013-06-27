#include "StdAfx.h"

#include "Tarcyl.h"
#include "TargetManager.h"

/* **
Purpose: to construct cylindrical target by populating sites on a rectangular lattice.

Input:
       A = cylinder length (in units of lattice spacing d)
       B = cylinder diameter (in units of lattice spacing d)
       ORI = 1. for cylinder axis in x_TF direction a1=(1,0,0),a2=(0,1,0)
             2. for cylinder axis in y_TF direction a1=(0,1,0),a1=(0,0,1)
             3. for cylinder axis in z_TF direction a1=(0,0,1),a2=(1,0,0)
       DX(1-3)=(dx,dy,dz)/d where dx,dy,dz=lattice spacings in x,y,z directions, and d=(dx*dy*dz)**(1/3)=effective lattice spacing
       IOSHP=device number for "target.out" file
            =-1 to suppress printing of "target.out"
       MXNAT = dimensioning information

Returns:
       A1(1-3) = unit vector along cylinder axis
       A2(1-3) = unit vector perpendicular to cylinder axis
       CDESCR = string describing target (up to 67 chars.)
       NAT = number of atoms in target
       IXYZ(1-NAT,1-3)=[x-X0(1)]/d,[y-X0(2)]/d,[z-X0(3)]/d for dipoles in target where x=y=z=0 at target centroid
       X0(1-3) = offset vector defined by above condition
       ICOMP(1-NAT,1-3)=1 (composition identifier)

Fortran history records removed.

Copyright (C) 1993,1995,1996,1997,1998,2000,2006,2007,2008, B.T. Draine and P.J. Flatau
Copyright (c) 2012, 2013, C++ version, V.Choliy
This code is covered by the GNU General Public License.
** */

Tarcyl::Tarcyl(TargetManager *man) : AbstractTarget(man) 
{
	shortDescr = string("Tarcyl");
	longDescr = string("Cylindrical target");
	xcm = ycm = zcm = (real)0.;
}

// determine centroid XCM,YCM,ZCM
// If A/d_x is near an odd int, XCM=0, even, 0.5
// If B/d_y is near an odd int, YCM=0, even, 0.5
// If B/d_z is near an odd int, ZCM=0, even, 0.5
void Tarcyl::PreSizer(int jx, int jy, int jz)
{
	xcm = (jx%2) ? zero_ : half_;
	ycm = (jy%2) ? zero_ : half_;
	zcm = (jz%2) ? zero_ : half_;

	real r2 = (real)0.25 * shpar[1] * shpar[1];
	minJx =  jx + 1;
	maxJx = -jx;
	minJy =  jy + 1;
	maxJy = -jy;
	minJz =  jz + 1;
	maxJz = -jz;

	switch((int)shpar[2])
	{
		case 1:
			for(int x=-jx; x<=jx+1; ++x)
			{
				real rx2 = twox_ * Fabs(x * dx.data[0] - xcm);
				if (rx2 <= shpar[0])
				{
					for(int y=-jy; y<=jy+1; ++y)
					{
						real ry2 = y * dx.data[1] - ycm;
						ry2 = ry2 * ry2;
						for(int z=-jz; z<=jz+1; ++z)
						{
							real rz2 = z * dx.data[2] - zcm;
							rz2 = rz2 * rz2;
							if (ry2 + rz2 <= r2)
							{
								InternalMinMax(x, y, z);
							}
						}
					}
				}
			}
			break;

		case 2:
			for(int x=-jx; x<=jx+1; ++x)
			{
				real rx2 = x * dx.data[0] - xcm;
				rx2 = rx2 * rx2;
				for(int y=-jy; y<=jy+1; ++y)
				{
					real ry2 = twox_ * Fabs(y * dx.data[1] - ycm);
					if (ry2 <= shpar[0])
					{
						for(int z=-jz; z<=jz+1; ++z)
						{
							real rz2 = z * dx.data[2] - zcm;
							rz2 = rz2 * rz2;
							if (rx2 + rz2 <= r2)
							{
								InternalMinMax(x, y, z);
							}
						}
					}
				}
			}
			break;

		case 3:
			for(int x=-jx; x<=jx+1; ++x)
			{
				real rx2 = x * dx.data[0] - xcm;
				rx2 = rx2 * rx2;
				for(int y=-jy; y<=jy+1; ++y)
				{
					real ry2 = y * dx.data[1] - ycm;
					ry2 = ry2 * ry2;
					if (rx2 + ry2 <= r2)
					{
						for(int z=-jz; z<=jz+1; ++z)
						{
							real rz2 = twox_ * Fabs(z * dx.data[2] - zcm);
							if (rz2 <= shpar[0])
							{
								InternalMinMax(x, y, z);
							}
						}
					}
				}
			}
			break;

		default:
			break;
	}
}

void Tarcyl::Sizer(void)
{
	dx = manager->CashedDx();
	switch((int)shpar[2])
	{
	case 1:						// Cylinder axis in x direction
		PreSizer(nint_(shpar[0] / dx.data[0]), nint_(shpar[1] / dx.data[1]), nint_(shpar[1] / dx.data[2]));
		break;

	case 2:						// Cylinder axis in y direction
		PreSizer(nint_(shpar[1] / dx.data[0]), nint_(shpar[0] / dx.data[1]), nint_(shpar[1] / dx.data[2]));
		break;

	case 3:						// Cylinder axis in z direction
		PreSizer(nint_(shpar[1] / dx.data[0]), nint_(shpar[1] / dx.data[1]), nint_(shpar[0] / dx.data[2]));
		break;

	default:
		break;
	}
	AllocateArrays(maxJx - minJx + 1, maxJy - minJy + 1, maxJz - minJz + 1);
}

void Tarcyl::Descriptor(void)
{
// really empty
}

void Tarcyl::Vector(void)
{
	a1.Clear();
	a2.Clear();
	a1.data[(int)shpar[2] - 1] = onex_;
	a2.data[(int)shpar[2] % 3] = onex_;
}

void Tarcyl::VectorX(void)						// Now redetermine location of centroid
{
	Vect3<real> xyzcm;
	for(int ja=0; ja<nat0; ++ja)
	{
		xyzcm += Vect3<real>((real)ixyz.Value(ja, 0), (real)ixyz.Value(ja, 1), (real)ixyz.Value(ja, 2));
	}
	x0 = -xyzcm / (real)nat0;
}

void Tarcyl::Allocator(void) // A=cylinder length/d, B=cylinder diameter/d
{
	int jlo, jhi;
	int curSize  = nx * ny;
	ixyz.Dimension(curSize, 3);
	real r2 = (real)0.25 * shpar[1] * shpar[1];
//
	switch((int)shpar[2])
	{
		case 1:
			jlo = maxJx;
			jhi = minJx;
			nat0 = 0;
			for(int jx=minJx; jx<=maxJx; ++jx)
			{
				real rx2 = twox_ * Fabs(jx * dx.data[0] - xcm);
				if (rx2 <= shpar[0])
				{
					for(int jy=minJy; jy<=maxJy; ++jy)
					{
						real ry2 = jy * dx.data[1] - ycm;
						ry2 = ry2 * ry2;
						for(int jz=minJz; jz<=maxJz; ++jz)
						{
							real rz2 = jz * dx.data[2] - zcm;
							rz2 = rz2 * rz2;
							if (ry2 + rz2 <= r2)
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
								if (jx < jlo) 
									jlo = jx;
								if (jx > jhi) 
									jhi = jx;
								++nat0;
							}
						}
					}
				}
			}
			break;

		case 2:
			jlo = maxJy;
			jhi = minJy;
			nat0 = 0;
			for(int jx=minJx; jx<=maxJx; ++jx)
			{
				real rx2 = jx * dx.data[0] - xcm;
				rx2 = rx2 * rx2;
				for(int jy=minJy; jy<=maxJy; ++jy)
				{
					real ry2 = twox_ * Fabs(jy * dx.data[1] - ycm);
					if (ry2 <= shpar[0])
					{
						for(int jz=minJz; jz<=maxJz; ++jz)
						{
							real rz2 = jz * dx.data[2] - zcm;
							rz2 = rz2 * rz2;
							if (rx2 + rz2 <= r2)
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
								if (jy < jlo)
									jlo = jy;
								if (jy > jhi) 
									jhi = jy;
								++nat0;
							}
						}
					}
				}
			}
			break;

		case 3:
			jlo = maxJz;
			jhi = minJz;
			nat0 = 0;
			for(int jx=minJx; jx<=maxJx; ++jx)
			{
				real rx2 = jx * dx.data[0] - xcm;
				rx2 = rx2 * rx2;
				for(int jy=minJy; jy<=maxJy; ++jy)
				{
					real ry2 = jy * dx.data[1] - ycm;
					ry2 = ry2 * ry2;
					if (rx2 + ry2 <= r2)
					{
						for(int jz=minJz; jz<=maxJz; ++jz)
						{
							real rz2 = twox_ * Fabs(jz * dx.data[2] - zcm);
							if (rz2 <= shpar[0])
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
								if (jz < jlo) 
									jlo = jz;
								if (jz > jhi) 
									jhi = jz;
			                    ++nat0;
							}
						}
					}
				}
			}
			break;

		default:
			break;
	}
	int nlay = jhi - jlo + 1;
	ixyz.Close(nat0);
//
// NLAY = number of layers in cylinder
// NFAC = number of atoms in slice
	int nfac = nat0 / nlay;
//
// REFF2=effective radius**2 of disk
	real reff2 = (real)nfac / Pi;
//
// ASPR=aspect ratio (length/diameter)
	real aspr = half_ * (real)nlay / Sqrt(reff2);
//
// Description
	sprintf(freeDescr, " Cyl.prism, NAT=%7d NFAC=%4d NLAY=%4d asp.ratio=%lf", nat0, nfac, nlay, aspr);
}

void Tarcyl::OutShpar(FILE *file)
{
	fprintf(file, " Length=%8.4lf Diameter=%8.4lf\n", shpar[0], shpar[1]);
}

const char *TargetVerboseDescriptor_Cylinder1(int num)
{
	static const char *descr[3] = 
	{
		"cylinder length", "cylinder diameter", "cylinder orientation, (1,2,3)<-(x,y,z)"
	};
	if (num < 3)
		return descr[num];
	else
		return NULL;
}

REGISTER_TARGET(Cylinder1,3,false,-1,0,"Homogenous isotropic cylinder")
void Target_Cylinder1::SayHello(FILE *stream)
{
	fprintf(stream, "The target %s parameters are:\n", shortDescr.c_str());
	fprintf(stream, "cylinder length and diameter:\n");
	fprintf(stream, "%20.16lf %20.16lf\n", shpar[0], shpar[1]);
    fprintf(stream, "cylinder orientation, %d\n", (int)shpar[2]);
}

const char *TargetVerboseDescriptor_Cylndrpbc(int num)
{
	static const char *descr[5] = 
	{
		"cylinder length (normally 1)", "cylinder diameter", "cylinder orientation, (1,2,3)<-(x,y,z)", "periodicity in y", "periodicity in z"
	};
	if (num < 5)
		return descr[num];
	else
		return NULL;
	return NULL;
}

REGISTER_TARGET(Cylndrpbc,5,false,3,0,"Circular slice (TUC for infinite cylinder)")
void Target_Cylndrpbc::SayHello(FILE *stream)
{
	fprintf(stream, " CYLNDRPBC: TUC for infinite cylinder\n");
	fprintf(stream, "The target %s parameters are:\n", shortDescr.c_str());
	fprintf(stream, "cylinder length/d (normally 1) and diameter:\n");
    fprintf(stream, "%20.16lf %20.16lf\n", shpar[0], shpar[1]);
	fprintf(stream, "next parameter s.b. (1,2,3) for axis (x,y,z): %d\n", (int)(shpar[2]));
	fprintf(stream, "next two parameters are periodicities in y and z directions, %20.16lf %20.16lf\n", shpar[3], shpar[4]);
}

void Target_Cylndrpbc::PreparePyzd()
{
	VectorA();
	pyd = shpar[3];
	pzd = shpar[4];
}

const char *TargetVerboseDescriptor_Uniaxicyl(int num)
{
	static const char *descr[3] = 
	{
		"cylinder length", "cylinder diameter", "cylinder orientation, (1,2,3)<-(x,y,z)"
	};
	if (num < 3)
		return descr[num];
	else
		return NULL;
}

REGISTER_TARGET(Uniaxicyl,3,false,-1,2,"Homogenous finite cylinder with uniaxial anisotropic dielectric tensor")
void Target_Uniaxicyl::SayHello(FILE *stream)
{
	fprintf(stream, "The target %s parameters are:\n", shortDescr.c_str());
	fprintf(stream, "cylinder length/d (normally 1) and diameter:\n");
	fprintf(stream, "%20.16lf %20.16lf\n", shpar[0], shpar[1]);
}

void Target_Uniaxicyl::PrepareIaniso()
{
	ianiso = TargetIsAnisotropic;
}

void Target_Uniaxicyl::Composer(int index, int item)
{
	icomp.Fill3(index, (short)(item + 1), (short)(item + 2), (short)(item + 2));
}