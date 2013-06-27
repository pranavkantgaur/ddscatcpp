#include "StdAfx.h"

#include "Definitions.h"
#include "DDscatMain.h"
#include "DDscatParameters.h"
#include "DielectricManager.h"
#include "LoadableTarget.h"
#include "Matrix.h"

void Alphadiag(AbstractTarget *currentTarget, Vect3<real> &akr, Matrix *theTensor, Vect3<Complex> &cxe0r, DielectricManager *dielec, int myid)
{
/* **
 Given:
		AKR(1-3)=(kx,ky,kz)*d, where d=effective lattice spacing
		BETADF(1-NAT), PHIDF(1-NAT), THETADF(1-NAT): orientation angles beta,phi,theta (radians)
		specifying orientation of "Dielectric Frame" (DF relative to the "Target Frame" (TF)

       CALPHA = polarizability prescription
              = 'LATTDR' for LDR of Draine & Goodman (1993)
              = 'GKDLDR' for LDR of Gutkowicz-Krusin & Draine (2004)
              = 'FLTRCD' for filtered coupled dipole approach
                             of Gay-Balmaz & Martin (2002) and Yurkin, Min & Hoekstra (2010)
              [= 'SCLDR' not supported in present version] 
       CMETHD = determines 3-d FFT routine used by ESELF
       CSHAPE = descriptor of target shape
       CXE0R(1-3) = polarization vector in lattice coordinates (assumed to be normalized)
       CXEPS(1-NCOMP)=distinct values of dielectric constant
       CXSC(1-NAT,1-3,1-3) = Complex scratch space for SCLDR calculatio
       CXSCR1(1-NAT,1-3)   = Complex scratch space for SCLDR calculatio
       CXZC = Complex scratch space needed by ESELF
       CXZW = Complex scratch space needed by ESELF
       DX(1-3)=(dx/d,dy/d,dz/d), where dx,dy,dz=lattice spacings in x,y,z directions, and d=(dx*dy*dz)**(1/3)
       IBETH= MYID+IBETH1 if first time through combined BETA/THETA orientation loop
       IBETH1=starting value of IBETH (see above)
       ICOMP(1-NAT3)=composition identifier for each lattice site and direction (ICOMP=0 if lattice site is unoccupied)
                storage scheme: 1x,2x,...,NATx,1y,2y,...,NATy,1z,2z,...,NATz
       IOCC(1-NAT) == 0 if site unoccupied
                   == 1 if site occupied
       IPHI = IPHI1 if first time through PHI orientation loop
       IPHI1= starting value of IPHI (see above)
       JORTH= 1 if first incident polarization state
       MXCOMP = dimensioning information
       MXN3 = dimensioning information
       MYID = parallel process identifier (=0 if only 1 process)
       NAT  = number of sites in extended target (incl. vacuum sites)
       NAT3 = 3*number of sites in extended target (incl. vacuum sites)
       NCOMP= number of different dielectric tensor elements in target
       NX   = x-dimension of extended target
       NY   = y-dimension of extended target
       NZ   = z-dimension of extended target
       SHPAR(1-10)=target shape parameters

 Returns:
       CXALPH(J,1-3)=(alpha_11,alpha_22,alpha_33)/d^3 for dipole J=1-NA where alpha_ij=Complex polarizability tensor.
       CXALOF(J,1-3)=(alpha_23,alpha_31,alpha_12)/d^3 for dipole J=1-NA

 If CALPHA = LATTDR:
    Compute dipole polarizability using "Lattice Dispersion Relation" of Draine & Goodman (1993,ApJ,March 10).  It is required that
    polarizability be such that an infinite lattice of such dipoles reproduce the continuum dispersion relation for radiation
    propagating with direction and polarization of radiation incident on the DDA target.

 If CALPHA = GKDLDR:
    Compute dipole polarizability using "Lattice Dispersion Relation" of Gutkowicz-Krusin and Draine (2004)..  It is required that
    polarizability be such that an infinite lattice of such dipoles reproduce the continuum dispersion relation for radiation
    propagating with direction and polarization of radiation incident on the DDA target.  This is the recommended option.  It is nearly
    but not exactly identical to LATTDR.

 If CALPHA = FLTRCD = "Filtered Discrete Dipole"
    Compute dipole polarizability for "Filtered Coupled Dipole" approach
    of Piller & Martin (1998) and Gay-Balmaz & Martin (2002)
    and recently discussed by Yurkin, Min, & Hoekstra (2010)

 Note: CXALPH = polarizability/d^3 
 In the event that ICOMP=0, then we set CXALPH=1.

B.T.Draine, Princeton Univ. Obs.
History records of Fortran versions removed.
Copyright (C) 1993,1996,1997,1998,1999,2003,2004,2006,2007,2008,2011
              B.T. Draine and P.J. Flatau
Copyright (C) 2012,2013, C++ versions, Choliy V.
This code is covered by the GNU General Public License.
** */
//
	const real zero_ = (real)0.;

	real ak2 = akr.ModSquared();
	real ak1 = Sqrt(ak2);
	real ak3 = ak1 * ak2;
	Complex cxrr = Complex((real)0., -ak3 / (real)1.5);
//
// Lattice dispersion relation (Draine & Goodman 1993)
	AlphaMethod alphaMethod = DDscatParameters::GetInstance()->Calpha();
	if (alphaMethod != AlphaMethod_End)
	{
		if (myid == 0)
		{
			char Buffer[256];
			real emkd = Sqrt(dielec->GetCxeps(0).mod()) * ak1;			// EMKD = |m|*k_0*d , where |m|=refractive index, k_0 = wave vector in vacuo, d = lattice spacing
			sprintf(Buffer, "%s Lattice dispersion relation for |m|k_0d=%6.4lf", AlphaEnumerator(alphaMethod), emkd);
			Wrimsg("Alpha ", Buffer);
		}
	}
//
	switch(alphaMethod)
	{
		case AlphaMethod_LATTDR:
		{
//
// Compute sum (a_j*e_j)^2 , where a_j=unit propagation vector e_j=unit polarization vector
			real sum = zero_;
			for(int l=0; l<3; ++l)
			{
				real tmp = Sqrt(cxe0r.data[l].re * cxe0r.data[l].re + cxe0r.data[l].im * cxe0r.data[l].im);
				sum += (akr.data[l] * tmp) * (akr.data[l] * tmp);
			}
			sum /= ak2;
			real b1 = -(real)(1.8915316 * ak2);
			real b2 =  (real)((0.1648469 - 1.7700004 * sum) * ak2);
			theTensor->AbstractAlphadiag(currentTarget, dielec->GetCxeps(), b1, b2, NULL, cxrr);
		}
			break;
//
// Lattice dispersion relation: modified (Gutkowicz-Krusin & Draine 2004)
		case AlphaMethod_GKDLDR:
		{
//
// B1 = (c_1/pi)*ak2 = (-5.9424219/pi)*ak2 = -1.8915316*ak2
// B2 = (c_2/pi)*ak2 = (0.5178819/pi)*ak2 = 0.1648469*ak2
// B3 = -[(3c_2+c_3)/pi]*ak2 = -[(3*0.5178819+4.0069747)/pi]*ak2 = -1.7700004*ak2
// B3L = B3*A(I)**2 where a_i = unit vector in direction of propagation
			real b31x[3];
			real b1 = -(real)(1.8915316 * ak2);
			real b2 =  (real)(0.1648469 * ak2);
			real b3 = -(real)(1.7700004 * ak2);
			for(int l=0; l<3; ++l)
			{
				b31x[l] = b3 * akr.data[l] * akr.data[l] / ak2;
			}
			theTensor->AbstractAlphadiag(currentTarget, dielec->GetCxeps(), b1, b2, b31x, cxrr);
		}
			break;

		case AlphaMethod_FLTRCD:
		{
// x = (kd)
// B1 = (4/3)*x^2 
// B2 = (2/3*pi)*ln[(pi-x)/(pi+x)]*x^3
//
// prescription for alpha from Yurkin, Min & Hoekstra (2010):
// *** Radiative-reaction correction
//     Note that we are applying it differently from other authors
//     (e.g., Piller & Martin 1998, Gay-Balmaz & Martin 2002, 
//     Yurkin, Min & Hoekstra 2010) who would have
//                 CXTERM=CXTERM/(1._WP+B1*CXTERM+CXTERM*CXRR)
//     whereas we write
//                  CXTERM=[CXTERM/(1._WP+B1*CXTERM)]/[1._WP+CXTERM*CXRR/(1._WP+B1*CXTERM)]
//     although to leading order (x^3) they are the same
			real b1 = ((real)4. * ak2 + ((real)2. / Pi) * Log((Pi - ak1) / (Pi + ak1)) * ak3) / (real)3.;
			theTensor->AbstractAlphadiag(currentTarget, dielec->GetCxeps(), b1, zero_, NULL, cxrr);
		}
			break;

		default:
			Wrimsg("Alpha ", "Error: invalid option for subroutine ALPHA");
			exit(0);
			break;
	}
//
// Now enter module to handle possible microcrystal rotation at each
// lattice site.  BETADF,PHIDF,THETADF = 3 rotation angles specifying
// orientation of "Dielectric Frame" (in which dielectric tensor is
// diagonal) relative to Target Frame.

	if (!currentTarget->IsLoadableTarget())
		return;

	theTensor->LoadableAlphadiag((LoadableTarget *)currentTarget);
}
