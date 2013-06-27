#include "StdAfx.h"

#include "Tarnsp.h"
#include "TargetManager.h"

/* **
Routine to construct multisphere target

Input:
       DIAMX =max extent in X direction/d
       PRINAX=0 to set A1=(1,0,0), A2=(0,1,0)
             =1 to use principal axes for vectors A1 and A2
       DX(1-3)=(dx,dy,dz)/d where dx,dy,dz=lattice spacings in x,y,z
               directions, and d=(dx*dy*dz)**(1/3)=effective lattice spacing
       CFLSHP= name of file containing locations and radii of spheres
       IOSHP =device number for "target.out" file
             =-1 to suppress printing of "target.out"
       MXNAT =dimensioning information (max number of atoms)

and, from input file CFLSHP:

       NSPH              = number of spheres 
	   *four* descriptive lines which will be ignored (may be blank)
       XS(1) YS(1) ZS(1) AS(1)
       XS(2) YS(2) ZS(2) AS(2)
       ...
       XS(NSPH) YS(NSPH) ZS(NSPH) AS(NSPH)

where XS(J),YS(J),ZS(J) = x,y,z coordinates in target frame, in arbitrary units, of center of sphere J
       AS(J)            = radius of sphere J (same arbitrary units)

units used for XS,YS,ZS, and AS are really arbitrary: actual size of
overall target, in units of lattice spacing d, is controlled by
parameter SHPAR(1) in ddscat.par

Output:
       A1(1-3)=(1,0,0)=unit vector defining target axis 1 in Target Frame
       A2(1-3)=(0,1,0)=unit vector defining target axis 2 in Target Frame
       X0(1-3)=location/d in TF of lattice site IXYZ=0 0 0
               TF origin is taken to be volume-weighted centroid of N spheres
       NAT=number of atoms in target
       IXYZ(1-NAT,1-3)=x/d,y/d,z/d for atoms of target
       ICOMP(1-NAT,1-3)=composition identifier (currently 1 1 1)
       CDESCR=description of target (up to 67 characters)

Fortran history records removed.

Copyright (C) 2000,2001,2003,2004,2007,2008,2009,2010,2012, B.T. Draine and P.J. Flatau
Copyright (c) 2012, C++ version, V.Choliy
This code is covered by the GNU General Public License.
** */

Tarnsp::Tarnsp(TargetManager *man) : LoadableTarget(man) 
{
	shortDescr = string("Tarnsp");
	longDescr = string("Cluster of spheres");
}

void Tarnsp::Sizer(void)
{
	Reader();
//
// Determine max extent in X direction:
	int ja;
	Item0 *work = elements->at(0);
	real xmin = work->Xs() - work->As();
	real xmax = work->Xs() + work->As();
	real ymin = work->Ys() - work->As();
	real ymax = work->Ys() + work->As();
	real zmin = work->Zs() - work->As();
	real zmax = work->Zs() + work->As();
	if (nsph > 1)
	{
		for(ja=1; ja<nsph; ++ja)
		{
			work = elements->at(ja);
			if (work->Xs() - work->As() < xmin) xmin = work->Xs() - work->As();
			if (work->Xs() + work->As() > xmax) xmax = work->Xs() + work->As();
			if (work->Ys() - work->As() < ymin) ymin = work->Ys() - work->As();
			if (work->Ys() + work->As() > ymax) ymax = work->Ys() + work->As();
			if (work->Zs() - work->As() < zmin) zmin = work->Zs() - work->As();
			if (work->Zs() + work->As() > zmax) zmax = work->Zs() + work->As();
		}
	}
	real scale = shpar[0] / (xmax - xmin);
//
// Now determine min,max values of I,J,K:
	dx = manager->CashedDx();
	const real delta = (real)0.0001;
	int lmx1 = nint_(scale * xmin / dx.data[0] - delta);
	int lmx2 = nint_(scale * xmax / dx.data[0] + delta);
	int lmy1 = nint_(scale * ymin / dx.data[1] - delta);
	int lmy2 = nint_(scale * ymax / dx.data[1] + delta);
	int lmz1 = nint_(scale * zmin / dx.data[2] - delta);
	int lmz2 = nint_(scale * zmax / dx.data[2] + delta);

	for(ja=0; ja<nsph; ++ja)
	{
		work = elements->at(ja);
		work->As2() = (scale * work->As()) * (scale * work->As());
		work->Xs() *= scale;
		work->Ys() *= scale;
		work->Zs() *= scale;
	}
// 
// Do min/max analysis
	minJx = lmx2;
	maxJx = lmx1;
	minJy = lmy2;
	maxJy = lmy1;
	minJz = lmz2;
	maxJz = lmz1;
	for(int jx=lmx1; jx<=lmx2; ++jx)
	{
		real x = (real)jx * dx.data[0];
		for(int jy=lmy1; jy<=lmy2; ++jy)
		{
			real y = (real)jy * dx.data[1];
			for(int jz=lmz1; jz<=lmz2; ++jz)
			{
				real z = (real)jz * dx.data[2];
				for(int ja=0; ja<nsph; ++ja)
				{
					real r2 = elements->at(ja)->DistSquared(x, y, z);
					if (r2 < elements->at(ja)->As2())
					{
						InternalMinMax(jx, jy, jz);
					}
				}
			}
		}
	}
	AllocateArrays(maxJx - minJx + 1, maxJy - minJy + 1, maxJz - minJz + 1);
}

void Tarnsp::Descriptor(void)			// Write target description into string CDESCR
{
	sprintf(freeDescr, " Multisphere cluster containing%7d dipoles", nat0);
}

void Tarnsp::Vector(void)
{
//
// Specify target axes A1 and A2
// If PRINAX=0, then
//     A1=(1,0,0) in target frame
//     A2=(0,1,0) in target frame
//     set ALPHA(1)=ALPHA(2)=ALPHA(3)=0
// If PRINAX=1., then
//     A1,A2 are principal axes of largest, second largest moment of inertia
	if (shpar[1] <= zero_)
	{
		VectorA();
		alpha.Clear();
	}
	else
		alpha = Prinaxis();
}

void Tarnsp::VectorX(void)
{
// Find volume-weighted centroid of the N spheres in IJK space
// Take this to be TF origin. Set X0 = -(IJK centroid)
	int ja;
	real z = zero_;
	for(ja=0; ja<nsph; ++ja)
	{
		z += elements->at(ja)->GetAsCube();
	}
	x0.Clear();
	for(ja=0; ja<nsph; ++ja)
	{
		Item0 *item = elements->at(ja);
		x0 -= (Vect3<real>(item->Xs(), item->Ys(), item->Zs()) * item->GetAsCube() / z);
	}
}

void Tarnsp::Reader()
{
	fprintf(stderr, ">TARNSP open file =%s\n", manager->GetCashedFileName());

	FILE *idvshpFile = fopen(manager->GetCashedFileName(), "r");
	if (idvshpFile)
	{
		char Buffer[256], Format[256];
		sprintf(Format, "%s%s%s%s", realFormat, realFormat, realFormat, realFormat); 
		fgets(Buffer, 255, idvshpFile);
		sscanf(Buffer, "%d", &nsph);
		int ja;
		for(ja=0; ja<4; ++ja)
		{
			fgets(Buffer, 255, idvshpFile);
		}
		if (!elements && nsph)
			elements = new vector<Item0 *>;
		for(ja=0; ja<nsph; ++ja)
		{
			fgets(Buffer, 255, idvshpFile);
			Item0 *item = new Item0;
			item->Sscanf(Buffer, Format);
			elements->push_back(item);
		}
		fclose(idvshpFile);
		fprintf(stderr, ">TARNSP close file=%s\n", manager->GetCashedFileName());
	}
	else
	{
		Wrimsg("Tarnas", "Error: Cannot open shape file");
	}
}

void Tarnsp::Allocator(void)				// Determine list of occupied sites
{
	int curSize = nx * ny;
	ixyz.Dimension(curSize, 3);

	int iwarn = 0;
	nat0 = 0;
	for(int jx=minJx; jx<=maxJx; ++jx)
	{
		real x = (real)jx * dx.data[0];
		for(int jy=minJy; jy<=maxJy; ++jy)
		{
			real y = (real)jy * dx.data[1];
			for(int jz=minJz; jz<=maxJz; ++jz)
			{
				real z = (real)jz * dx.data[2];
				bool occ = false;
				for(int ja=0; ja<nsph; ++ja)
				{
					real r2 = elements->at(ja)->DistSquared(x, y, z);
					if(r2 < elements->at(ja)->As2())
					{
						if(occ)
						{
							char cmsgnm[64];
							iwarn++;
							if(iwarn <= 10)
							{
								sprintf(cmsgnm, "overlap for sphere%6d", ja);
								Wrimsg("Tarnsp", cmsgnm);
							}
							if(iwarn == 10)
							{
								sprintf(cmsgnm, "iwarn=10: further overlap warnings suppressed...");
								Wrimsg("Tarnsp", cmsgnm);
							}
						}
						occ = true;
					}
				}
				if (occ)			// Site is occupied:
				{
					if (nat0 >= curSize)
					{
						ixyz.Extend(nz);
						curSize = ixyz.GetSize(0);
					}
					ixyz.Fill3(nat0, jx, jy, jz);
					int index = GetLinearAddress(nat0);
					Composer(index, 0);					// Homogeneous target
					iocc[index] = true;
					++nat0;
				}
			}
		}
	}
	ixyz.Close(nat0);
	Descriptor();
}

const char *TargetVerboseDescriptor_SpheresN(int num)
{
	static const char *descr[3] = 
	{
		"target extent in x-direction", "0 for (1,0,0)&(0,1,0); 1 for principal axes", "sphere parameter file name"
	};
	if (num < 3)
		return descr[num];
	else
		return NULL;
}

REGISTER_TARGET(SpheresN,2,true,-1,1,"union of N spheres")
void Target_SpheresN::SayHello(FILE *stream)
{
	fprintf(stream, "The target %s parameters are:\n", shortDescr.c_str());
	fprintf(stream, "DIAMX = target extent in x-direction:\n");
	fprintf(stream, "%20.16lf\n", manager->GetShpar(0));
	fprintf(stream, "PRINAX=0 for a_1=(1,0,0),a_2=(0,1,0) in TF, =1 for a_1,a_2 = principal axes\n");
	fprintf(stream, "%20.16lf\n", shpar[1]);
	fprintf(stream, "name of sphere parameter file (e.g., tarnsp.par [in quotes])\n");
	fprintf(stream, "%s\n", cflshp);
}
