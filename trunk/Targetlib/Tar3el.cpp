#include "StdAfx.h"

#include "Tar3el.h"
#include "TargetManager.h"

/* **
Routine to construct three collinear touching ellipsoids from "atoms"
Input:
       AX=(x-length of one ellipsoid)/d    (d=lattice spacing)
       AY=(y-length of one ellipsoid)/d
       AZ=(z-length of one ellipsoid)/d
       DX(1-3)=lattice spacing(dx,dy,dz)/d [d=(dx*dy*dz)**(1/3)]
       IOSHP=device number for "target.out" file
            =-1 to suppress printing of "target.out"
       MXNAT=dimensioning information (max number of atoms)
Output:
       A1(1-3)=(1,0,0)=unit vector defining target axis 1 in Target Fr
       A2(1-3)=(0,1,0)=unit vector defining target axis 2 in Target Fr
       NAT=number of dipoles in target
       IXYZ(1-NAT,1-3)=x/d,y/d,z/d for target dipole locations
       ICOMP(1-NAT,1-3)=x,y,z composition for dipole locations
                       = 1 for locations in first ellipsoid
                         2                  second
                         3                  third
       CDESCR=description of target (up to 67 characters)
Note: atoms 1       - NAT/3   are in first ellipsoid
            NAT/3+1 - 2*NAT/3 are in second ellipsoid
            2*NAT/3+1 - NAT   are in third ellipsoid
Ellipsoids are displaced from one another in x-direction
First ellipsoid is at smallest x values, third at largest

B.T.Draine, Princeton Univ. Obs.
Fortran history records removed.

Copyright (C) 1993,1996,1997,1998,2000,2007,2008 B.T. Draine and P.J. Flatau
Copyright (C) 2012, 2013, C++ version, Choliy V.
This code is covered by the GNU General Public License.
** */

Tar3el::Tar3el(TargetManager *man) : TarNel(man) 
{
	shortDescr = string("Tar3el");
	longDescr = string("Ellipsoidal grain");
	mySize = 3;
	delta = 0;
}

const char *TargetVerboseDescriptor_AniEll3(int num)
{
	static const char *descr[3] = 
	{
		"Length of ellipsoid in x direction", "Length of ellipsoid in y direction", "Length of ellipsoid in z direction"
	};
	if (num < 3)
		return descr[num];
	else
		return NULL;
}

REGISTER_TARGET(AniEll3,3,false,-1,9,"Three touching identical anisotropic ellipsoids")
void Target_AniEll3::SayHello(FILE *stream)
{
	fprintf(stream, "The target %s parameters are:\n", shortDescr.c_str());
	fprintf(stream, "length of one ellipsoid in x, y, z directions are:\n");
	fprintf(stream, "%lf %lf %lf\n", shpar[0], shpar[1], shpar[2]);
}

void Target_AniEll3::Composer(int index, int item)
{
	AnisotropicComposer(index, item);
}

void Target_AniEll3::PrepareIaniso(void)
{
	ianiso = TargetIsAnisotropic;
}

const char *TargetVerboseDescriptor_Ellipso3(int num)
{
	static const char *descr[3] = 
	{
		"x length/d for one ellipsoid", "y length/d for one ellipsoid", "z length/d for one ellipsoid"
	};
	if (num < 3)
		return descr[num];
	else
		return NULL;
}

REGISTER_TARGET(Ellipso3,3,false,-1,3,"Three touching ellipsoids; second ellipsoid displaced from firstst in +x-direction, third displaced from second in +x-direction")
void Target_Ellipso3::SayHello(FILE *stream)
{
	fprintf(stream, "The target %s parameters are:\n", shortDescr.c_str());
	fprintf(stream, "x,y,z length/d for one ellipsoid:\n");
	fprintf(stream, "%lf %lf %lf\n", shpar[0], shpar[1], shpar[2]);
}

void Target_Ellipso3::Composer(int index, int item)
{
	IsotropicComposer(index, item);
}
