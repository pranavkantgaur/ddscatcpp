#include "StdAfx.h"

#include "Definitions.h"
#include "DDscatMain.h"
#include "AbstractTarget.h"

void Evale(AbstractTarget *curTarget, Vect3<Complex> &cxe00, Vect3<real> &akd, Complex *cxe)
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
This code is covered by the GNU General Public License.
** */
//
// Evaluate electric field vector at each dipole location.
// If NAT=NAT0, then evaluate E only at occupied sites.
// If NAT>NAT0, then evaluate E at all sites.

	int curTaregtNat = curTarget->Nat();
	int curTaregtNat0 = curTarget->Nat0();
	double *akdData = akd.Data();
	double *dxData = curTarget->Dx().Data();
	double *x0Data = curTarget->X0().Data();

	int ia, m;
	if (curTargetNat == curTargetNat0)
	{
		for(ia=0; ia<curTargetNat0; ++ia)
		{
			real x = zero_;
			for(m=0; m<3; ++m)
			{
				x += (akdData[m] * dxData[m] * ((real)(curTarget->Ixyz().Value(ia, m)) + x0Data[m]));
			}
			Complex cxfac = Complex(Cos(x), Sin(x));
			for(m=0; m<3; ++m)
			{
				cxe[ia + curTargetNat0 * m] = cxe00.Data(m) * cxfac;
			}
		}
	}
	else
	{
		ia = 0;
		for(int iz=0; iz<curTarget->Nz(); ++iz)
		{
			real x1 = akdData[2] * dxData[2] * (x0Data[2] + iz);
			for(int iy=0; iy<curTarget->Ny(); ++iy)
			{
				real x2 = x1 + akdData[1] * dxData[1] * (x0Data[1] + iy);
				for(int ix=0; ix<curTarget->Nx(); ++ix)
				{
					real x = x2 + akdData[0] * dxData[0] * (x0Data[0] + ix);
					Complex cxfac = Complex(Cos(x), Sin(x));
					for(m=0; m<3; ++m)
					{
						cxe[ia + curTargetNat0 * m] = cxe00.Data(m) * cxfac;
					}
					++ia;
				}
			}
		}
	}
}
