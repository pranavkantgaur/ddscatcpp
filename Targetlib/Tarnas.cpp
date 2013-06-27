#include "StdAfx.h"

#include "Tarnas.h"
#include "TargetManager.h"

/* **
Routine to construct multisphere target of general anisotropic materials (target option 'NANSPH')

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

       N              = number of spheres
       descriptive line which will be ignored (may be blank)
       descriptive line which will be ignored (may be blank)
       descriptive line which will be ignored (may be blank)
       descriptive line which will be ignored (may be blank)
       XS(1) YS(1) ZS(1) AS(1) IC1(1) IC2(1) IC3(1) TH(1) PH(1) BE(1)
       XS(2) YS(2) ZS(2) AS(2) IC1(2) IC2(2) IC3(2) TH(2) PH(2) BE(2)
       ...
       XS(N) YS(N) ZS(N) AS(N) IC1(N) IC2(N) IC3(N) TH(N) PH(N) BE(N)

where XS(J),YS(J),ZS(J) = x,y,z coordinates in target frame, in arbitrary units, of center of sphere J
       AS(J)            = radius of sphere J (same arbitrary units)
       IC1(J),IC2(J),IC3(J) = dielectric function identifier for directions 1,2,3 in DF
       TH(J),PH(J),BE(J) = angles theta_epsilon,phi_epsilon,beta_epsil specifying orientation of DF in TF

units used for XS,YS,ZS, and AS are really arbitrary: actual size of
overall target, in units of lattice spacing d, is controlled by parameter SHPAR(1) in ddscat.par

Output:
       A1(1-3)=(1,0,0)=unit vector defining target axis 1 in Target Frame
       A2(1-3)=(0,1,0)=unit vector defining target axis 2 in Target Frame
       X0(1-3)=location/d in TF of lattice site IXYZ=0 0 0 TF origin is taken to be located at centroid of N spheres.
       NAT=number of atoms in target
       IXYZ(1-NAT,1-3)=x/d,y/d,z/d for atoms of target
       ICOMP(1-NAT,1-3)=composition identifier
       BETADF(1-NAT)=angle beta_epsilon (radians) specifying orientation of Dielectric Frame (DF) rel to Target Frame (TF)
			See UserGuide for definition of angles theta_epsilon phi_epsilon, beta_epsilon and discussion.
       PHIDF(1-NAT)=angle phi_epsilon (radians) specifying orientation DF relative to TF
       THETADF(1-NAT)=angle theta_epsilon (radians) specifying orientation of DF relative to TF
       CDESCR=description of target (up to 67 characters)

Fortran history records removed.

Copyright (C) 2004,2007,2008,2009,2010 B.T. Draine and P.J. Flatau
Copyright (C) 2012, C++ version, Choliy V.

This code is covered by the GNU General Public License.
** */

Tarnas::Tarnas(TargetManager *man) : LoadableTarget(man)			// Write target description into string CDESCR
{
	shortDescr = string("Tarnas");
	longDescr = string(" Multisphere cluster");
	nsph = 0;
}

void Tarnas::Sizer(void)
{
	Reader();
//
// Check for overlap and, if so, possible conflict in composition or (if anisotropic) orientation
	int ja, jy;
	if (nsph >= 2)
	{
		for(ja=1; ja<nsph; ++ja)
		{
			Item *one = (Item *)(elements->at(ja));
			for(jy=0; jy<ja; ++jy)
			{
				Item *two = (Item *)(elements->at(jy));
				real r2 = one->Dist2(two);
//
// Exact R2 should be nonnegative, but allow for roundoff error
				if (r2 < -0.0001 * (one->As() + two->As()))
				{
					if (one->Overlap(two))
					{
						Errmsg("Fatal", "Tarnas", " Sphere overlap but differing composition");
					}
					if (two->IsAniso())
					{
						if (one->IsDisoriented(two))
							Errmsg("Fatal", "Tarnas", " Sphere overlap but have differing orientation");
					}
				}
			}
		}
	}
//
// Determine max extent in X direction:
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
	const real delta = (real)0.01;
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
		real x = jx * dx.data[0];
		for(int jy=lmy1; jy<=lmy2; ++jy)
		{
			real y = jy * dx.data[1];
			for(int jz=lmz1; jz<=lmz2; ++jz)
			{
				real z = jz * dx.data[2];
				for(int ja=0; ja<nsph; ++ja)
				{
					real r2 = elements->at(ja)->DistSquared(x, y, z);
					if (r2 < elements->at(ja)->As2())
					{
						InternalMinMax(jx, jy, jz);
						break;
					}
				}
			}
		}
	}
//
	AllocateArrays(maxJx - minJx + 1, maxJy - minJy + 1, maxJz - minJz + 1);
}

void Tarnas::Descriptor(void)
{
	sprintf(freeDescr, "%s containing %d dipoles", longDescr.c_str(), nat0);
}

void Tarnas::Vector(void)
{
//
// Specify target axes A1 and A2
// If PRINAX=0, then
//     A1=(1,0,0) in target frame
//     A2=(0,1,0) in target frame
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

void Tarnas::VectorX(void)
{
//
// Find volume-weighted centroid of the N spheres in IJK space
// Take this to be TF origin. Set X0  = -(IJK centroid)
	int j;
	x0.Clear();
	real z = zero_;
	for(j=0; j<nsph; ++j)
	{
		z += elements->at(j)->GetAsCube();
	}
	for(j=0; j<nsph; ++j)
	{
		Item0 *item = elements->at(j);
		x0 -= Vect3<real>(item->Xs(), item->Ys(), item->Zs()) * item->GetAsCube() / z;
	}
}

void Tarnas::Reader(void)
{
	fprintf(stderr, " >TARNAS open file=%s\n", manager->GetCashedFileName());

	FILE *idvshpFile = fopen(manager->GetCashedFileName(), "r");
	if (idvshpFile)
	{
		char Buffer[256], Format[256];
		sprintf(Format, "%s%s%s%s%%d%%d%%d%s%s%s", realFormat, realFormat, realFormat, realFormat, realFormat, realFormat, realFormat);
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
			Item *item = new Item;
			item->Sscanf(Buffer, Format);
			elements->push_back(item);
		}
		fclose(idvshpFile);
		fprintf(stderr, " >TARNAS close file=%s\n", manager->GetCashedFileName());
	}
	else
	{
		Wrimsg("Tarnas", "Error: Cannot open shape file");
	}
}

void Tarnas::Allocator(void)					// Determine list of occupied sites
{
	int curSize = nx * ny;
	ixyz.Dimension(curSize, 3);

	nat0 = 0;
	for(int jx=minJx; jx<=maxJx; ++jx)
	{
		real x = jx * dx.data[0];
		for(int jy=minJy; jy<=maxJy; ++jy)
		{
			real y = jy * dx.data[1];
			for(int jz=minJz; jz<=maxJz; ++jz)
			{
				real z = jz * dx.data[2];
				bool occ = false;
				int js = -1;
				for(int ja=0; ja<nsph; ++ja)
				{
					real r2 = elements->at(ja)->DistSquared(x, y, z);
					if (r2 < elements->at(ja)->As2())
					{
						occ = true;
						js = ja;
					}
				}
				if (occ)									// Site is occupied:
				{
					if (nat0 >= curSize)
					{
						ixyz.Extend(nz);
						curSize = ixyz.GetSize(0);
					}
					ixyz.Fill3(nat0, jx, jy, jz);
					int index = GetLinearAddress(nat0);
					Item *item = (Item *)(elements->at(js));
					icomp.Fill3(index, (short)item->Ic1(), (short)item->Ic2(), (short)item->Ic3());
					iocc[index] = true;
					if (dfdata == NULL)
					{
						dfdata = new AbstractDFData(elements);
						dfdata->Allocate(nx*ny*nz);
					}
					dfdata->SetAngles(index, js);
					++nat0;
				}
			}
		}
	}
	ixyz.Close(nat0);
	Descriptor();
}

const char *TargetVerboseDescriptor_Sphrn_pbc(int num)
{
	static const char *descr[5] = 
	{
		"maximum target extent in x-direction", "0 for (1,0,0)&(0,1,0); 1 for principal axes", "y/d period", "z/d period", "sphere parameter file name"
	};
	if (num < 5)
		return descr[num];
	else
		return NULL;
}

REGISTER_TARGET(Sphrn_pbc,3,true,1,0,"N possibly anisotropic spheres (one TUC)")
void Target_Sphrn_pbc::SayHello(FILE *stream)
{
	fprintf(stream, "The target %s parameters are:\n", shortDescr.c_str());
	fprintf(stream, "Enter diam_x/d = target extent in x-direction:\n");
	fprintf(stream, "%20.16lf\n", shpar[0]);
	fprintf(stream, "Prinaxis = 0\n");
	fprintf(stream, "filename for sphere locations (in quotes)\n");
	fprintf(stream, "%s\n", cflshp);
	fprintf(stream, "next two parameters are y and z periods, %20.16lf %20.16lf\n", shpar[2], shpar[3]);
}

void Target_Sphrn_pbc::PreparePyzd()
{
	pyd = shpar[2] * Pow((real)(0.75 * nat0 / Pi), (real)(1. / 3.)) / dx.data[1];
	pzd = shpar[3] * Pow((real)(0.75 * nat0 / Pi), (real)(1. / 3.)) / dx.data[2];
}

const char *TargetVerboseDescriptor_SphAniN(int num)
{
	static const char *descr[3] = 
	{
		"maximum target extent in x-direction", "0 for (1,0,0)&(0,1,0); 1 for principal axes", "sphere parameter file name"
	};
	if (num < 3)
		return descr[num];
	else
		return NULL;
}

REGISTER_TARGET(SphAniN,2,true,-1,0,"N-spheres of anisotropic materials")
void Target_SphAniN::SayHello(FILE *stream)
{
	fprintf(stream, "The target %s parameters are:\n", shortDescr.c_str());
	fprintf(stream, "PRINAX=0 for a_1=(1,0,0),a_2=(0,1,0) in TF, =1 for a_1,a_2 = principal axes\n");
	fprintf(stream, "%20.16lf\n", shpar[1]);
	fprintf(stream, "name of sphere parameter file (e.g., tarnsp.par [in quotes])\n");
	fprintf(stream, "%s\n", cflshp);
	fprintf(stream, "DIAMX = target extent in x-direction (0 to stop)\n");
	fprintf(stream, "%20.16lf\n", shpar[0]);
}

void Target_SphAniN::PrepareIaniso()
{
	PrepareIanisoSpecial();
}
