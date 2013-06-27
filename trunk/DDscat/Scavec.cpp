#include "StdAfx.h"

#include "Definitions.h"
#include "Complex.h"
#include "Vect3.h"

void Scavec(int nscat, real *thetan, real *phin, Vect3<real> *ensc, Vect3<real> *em1, Vect3<real> *em2)
{
/* **
Given:
      MXSCA=dimension information for ENSC,PHIN,THETAN
      NSCAT=number of scattering directions
      THETAN(1-NSCAT)=scattering angles theta
      PHIN(1-NSCAT)=scattering angles phi
      CXE01(1-3)=Complex polarization vector 1 (phi=0 direction is defined by x,y plane (Lab Frame), where incident radiation propagates along the x axis.

Returns:
      ENSC(1-3,1-NSCAT)=scattering vectors in Lab Frame
      EM1(1-3,1-NSCAT)=scattered pol vectors parallel to scat. plane in Lab Frame
      EM2(1-3,1-NSCAT)=scattered pol vectors perp. to scat. plane in Lab Frame
 It is assumed that incident propagation vector is in x-direction in Lab Frame

 History:
 96.11.06 (BTD): Changed definition of scattering angle phi
                 Previously, phi was measured from plane containing
                 incident k vector (i.e., Lab x-axis) and Re(CXE01)
                 Henceforth, phi is measured from Lab x,y plane.
 10.01.30 (BTD): cosmetic changes
 end history
Copyright (C) 1993,1996,2010 B.T. Draine and P.J. Flatau
This code is covered by the GNU General Public License.
** */

	for(int i=0; i<nscat; ++i)
	{
		real cosphi = Cos(phin[i]);
		real sinphi = Sin(phin[i]);
		real costhe = Cos(thetan[i]);
		real sinthe = Sin(thetan[i]);
		ensc[i].Set(costhe, sinthe*cosphi, sinthe*sinphi);
		em1[i].Set(-sinthe, costhe*cosphi, costhe*sinphi);
		em2[i].Set((real)0., -sinphi, cosphi);
	}
}