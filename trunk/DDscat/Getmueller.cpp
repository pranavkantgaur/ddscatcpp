#include "StdAfx.h"

#include "Definitions.h"
#include "DDscatMain.h"
#include "DDscatParameters.h"

void Getmueller(DDscatParameters *param, int ibeta, int iphi, int itheta, PeriodicBoundaryFlag jpbc, real *orderm, real *ordern, real ak1, 
	Vect3<real> *aks_tf, Vect3<real> *ensc_lf, Vect3<real> *ensc_tf, real pyddx, real pzddx, real *sm, real *smori, real *s1111, real *s2121, 
	Complex *cx1121, FourArray *cxfData, SumPackage &sumPackage, Vect3<real> *em1_lf, Vect3<real> *em2_lf, Vect3<real> *em1_tf, Vect3<real> *em2_tf)
{
/* **
 Given:
   IBETA =index specifying rotation of target around axis a1
   IORTH =number of incident polarization states being calculated for
   IPHI  =index specifying phi value for orientation of target axis a1
   ITHETA=index specifying theta value for orientation of target axis a1
   JPBC  = 0 if PBC not used
           1 if PBC used in y direction only
           2 if PBC used in z direction only
           3 if PBC used in both y and z directions
   MXBETA=dimensioning information
   MXSCA =dimensioning information
   MXPHI =dimensioning information
   MXTHET=dimensioning information
   NSCAT =number of scattering directions for evaluation of Mueller matrix
   CMDTRQ='DOTORQ' or 'NOTORQ' depending on whether or not torque calculation is to be carried out
   ENSC_LF(1-3,1-NSCAT)=normalized scattering direction vector in Lab Frame for scattering directions 1-NSCAT
   ENSCR(1-3,1-NSCAT)=normalized scattering direction vector in Target Frame for scattering directions 1-NSCA
   EM1_LF(1-3,1-NSCAT)=unit scattering polarization vector 1 in Lab Frame
   EM2_LF(1-3,1-NSCAT)=unit scattering polarization vector 2 in Lab Frame
   EM1_TF(1-3,1-NSCAT)=unit scattering polarization vector 1 in Target Frame
   EM2_TF(1-3,1-NSCAT)=unit scattering polarization vector 2 in Target Frame
   PHIN(1-NSCAT)=values of PHI for scattering directions 1-NSCAT in Lab Frame [only used if JPBC=0: scattering by finite target]
   AKS_TF(1-3,1-NSCAT)=(kx,ky,kz)*d for NSCAT scattering directions in Target Frame
   ORDERM(1-NSCAT)=scattering order M for case of 1-d or 2-d target
   ORDERN(1-NSCAT)=scattering order N for case of 2-d target
   S1111(1-NSCAT)=weighted sum of |f_11|^2 over previous orientations
   S2121(1-NSCAT)=weighted sum of |f_21|^2 over previous orientations
   CXE01_LF(1-3)    =incident pol state 1 in Lab Frame
   CXE02_LF(1-3)    =                   2
   CXF11(1-NSCAT)=f_11 for current orientation
   CXF12(1-NSCAT)=f_12 for current orientation
   CXF21(1-NSCAT)=f_21 for current orientation
   CXF22(1-NSCAT)=f_22 for current orientation
   QABS(2)       =Q_abs for two incident polarizations
   QABSUM(2)     =weighted sum of Qabs over previous orientations
   QBKSCA(2)     =Q_bk for two incident polarizations
   QBKSUM(2)     =weighted sum of Q_bk over previous orientations
   QEXSUM(2)     =weighted sum of Q_ext over previous orientations
   QEXT(2)       =Q_ext for two incident polarizations
   QPHA(2)       =Q_pha for two incident polarizations
   QPHSUM(2)     =weighted sum of Q_ph over previous orientations
   QSCAG(3,2)    =Q_pr vector for two incident polarizations
   QSCAT(2)      =Q_sca for two incident polarizations
   QSCGSUM(3,2)  =weighted sum of Q_pr vector over prev. orients.
   QSCSUM(2)     =weighted sum of Q_sca over previous orientations
   QTRQAB(3,2)   =absorptive part of Q_gamma for two inc. pols.
   QTRQABSUM(3,2)=weighted sum of Q_gamma over previous orientations
   QTRQSC(3,2)   =scattering part of Q_gamma for two inc. pols.
   QTRQSCSUM(3,2)=weighted sum of Q_gamma over previous orientations
   WGTA(ITHETA,IPHI)=weight factor for theta,phi orientation
   WGTB(IBETA)      =weight factor for beta orientation
   SMORI(1-NSCAT,4,4)=weighted sum of Mueller matrix over previous orientations
 If IORTH=1, returns:
   S1111(1-NSCAT)=updated weighted sum of |f_11|^2 over NSCAT scattering directions
   S2121(1-NSCAT)=updated weighted sum of |f_21|^2
   CX1121(1-NSCAT)=updated weighted summ of f_11*conjg(f_21)

If IORTH=2, returns:
   SM(1-NSCAT,4,4)=Mueller matrix for current orientation, and NSCAT scattering directions
   CXS1(1-NSCAT)=amplitude scattering matrix element S_1 for current orientation
   CXS2(1-NSCAT)=amplitude scattering matrix element S_2 for current orientation
   CXS3(1-NSCAT)=amplitude scattering matrix element S_3 for current orientation
   CXS4(1-NSCAT)=amplitude scattering matrix element S_4 for current orientations
   QABSUM(2)     =updated weighted sum of Q_abs
   QBKSUM(2)     =updated weighted sum of Q_bk
   QEXSUM(2)     =updated weighted sum of Q_ext
   QSCGSUM(3,2)  =updated weighted sum of Q_pr vector
   QSCSUM(2)     =updated weighted sum of Q_sca
   QTRQABSUM(3,2)=updated weighted sum of absorptive contribution to Q_gamma
   QTRQSCSUM(3,2)=updated weighted sum of scattering contribution to Q_gamma

   SMORI(1-NSCAT,4,4)=updated weighted sum of Mueller matrix elements

History:
Fortran versions history removed.

Copyright (C) 1996,1997,1998,2000,2003,2006,2007,2008,2011,2012
              B.T. Draine and P.J. Flatau
This code is covered by the GNU General Public License.
** */

/* **
!*********  Sum scattering properties over orientations ****************
Correspondence to notation of Bohren & Huffman (1983):
Our f_jk correspond to notation of Draine (1988).
For the special case where incident polarization modes 1 and 2 are
arallel and perpendicular to the scattering plane, we have the
relationship between the f_jk and the elements S_j of the amplitude
scattering matrix as defined by Bohren & Huffman (1983) is as follows:
    S_1 = -i*f_22
    S_2 = -i*f_11
    S_3 =  i*f_12
    S_4 =  i*f_21
In this case (incident pol. states parallel,perpendicular to the
scattering plane), the elements of the 4x4 Mueller matrix (see Bohren
Huffman 1983) are related to our f_jk as follows:
    S_11 = [|f_11|^2 + |f_22|^2 + |f_12|^2 + |f_21|^2]/2
    S_12 = [|f_11|^2 + |f_21|^2 - |f_22|^2 - |f_12|^2]/2
    S_13 = -Re[f_11*conjg(f_12) + f_22*conjg(f_21)]
    S_14 =  Im[f_22*conjg(f_21) - f_11*conjg(f_12)]
    S_21 = [|f_11|^2 + |f_12|^2 - |f_22|^2 - |f_21|^2]/2
    S_22 = [|f_11|^2 + |f_22|^2 - |f_21|^2 - |f_12|^2]/2
    S_23 =  Re[f_22*conjg(f_21) - f_11*conjg(f_12)]
    S_24 = -Im[f_11*conjg(f_12) + f_22*conjg(f_21)]
    S_31 = -Re[f_11*conjg(f_21) + f_22*conjg(f_22)]
    S_32 =  Re[f_22*conjg(f_12) - f_22*conjg(f_21)]
    S_33 =  Re[f_22*conjg(f_11) + f_12*conjg(f_21)]
    S_34 =  Im[f_11*conjg(f_22) + f_21*conjg(f_22)]
    S_41 = -Im[f_21*conjg(f_11) + f_22*conjg(f_12)]
    S_42 =  Im[f_22*conjg(f_12) - f_22*conjg(f_11)]
    S_43 =  Im[f_22*conjg(f_11) - f_12*conjg(f_21)]
    S_44 =  Re[f_22*conjg(f_11) - f_12*conjg(f_21)]

Notation internal to this program:
ij.ne.kl: CX_ijkl = \f_ij* \times f_kl  (summed over orientations)
           S_ijij = \f_ij* \times f_ij  (summed over orientations)
4 pure real elements: S1111,S1212,S2121,S2222
8 Complex elements:   CX1112,CX1121,CX1122,CX1221,CX1222,CX2122
4 redundant elts.:    CX1211 = Conjg(CX1112)
                      CX2111 = Conjg(CX1121)
                      CX2112 = Conjg(CX1221)
                      CX2211 = Conjg(CX1122)

In this routine we will carry out fully general calculation of the
4x4 Mueller matrix elements, valid for arbitrary incident polarization
states j=1,2 for which we have the scattering matrix f_ij

First, we compute the amplitude scattering matrix elements S_1,S_2,
S_3,S_4 as defined by Bohren & Huffman:
** */

// Notation:
//        khat_0     = unit vector along incident direction
//        ehat_01    = CXE01 = complex incident polarization vector 1
//        ehat_02    = CXE02 = complex incident polarization vector 2
//        ehat_0perp = incident pol vector perp to scattering plane scattered pol vector perp to scattering plane em2
//        ehat_0para = incident pol vector para to scattering plane
//                   = khat_0 cross ehat_0perp   [conventional def.]
//                   = xhat_LF cross ehat_0perp
//                   = xhat_LF cross [yhat_LF (yhat_LF dot ehat_0perp) + zhat_LF (zhat_LF dot ehat_0perp)]
//                   = - yhat_LF (ehat_0perp dot zhat_LF) + zhat_LF (ehat_0perp dot yhat_LF) - yhat_LF em2_LF(3) + zhat_LF em2_LF(2)
//        em2        = scattered pol vector perp to scattering plane
//        cxa        = conjg(ehat_01) dot yhat_LF
//        cxb        = conjg(ehat_01) dot zhat_LF
//        cxc        = conjg(ehat_02) dot yhat_LF
//        cxd        = conjg(ehat_02) dot zhat_LF
//
// WG=weighting factor for this orientation

	const Complex conei((real)0., (real)1.);
	const real half_ = (real)0.5;
	Complex cxa, cxb, cxc, cxd, cxfac0, cxfac1, cxfac2;
	real fac;
    int nd, k;

	real wg = param->GetOriData()->GetWgtaValue(itheta, iphi) * param->GetOriData()->GetWgtb(ibeta);
	if (param->Iorth() == 2)
	{
		FourArray *myData = new FourArray;
		myData->Alloc(param->Nscat());
		Complex *myS1 = myData->Cxs1();
		Complex *myS2 = myData->Cxs2();
		Complex *myS3 = myData->Cxs3();
		Complex *myS4 = myData->Cxs4();
//
// Compute Complex scattering amplitudes S_1,S_2,S_3,S_4 for this particular target orientation and NSCAT scattering directions:
		switch(jpbc)				// ChB: TODO the code snippets are identical for every jpbc
		{
		case PeriodicNo:
			cxa = param->Cxe01_lf().data[1].conjg();
            cxb = param->Cxe01_lf().data[2].conjg();
            cxc = param->Cxe02_lf().data[1].conjg();
            cxd = param->Cxe02_lf().data[2].conjg();
			for(nd=0; nd<param->Nscat(); ++nd)
			{
				real cosphi =  em2_lf[nd].data[2];
				real sinphi = -em2_lf[nd].data[1];
				Complex cxaa = cxa * cosphi + cxb * sinphi;
				Complex cxbb = cxb * cosphi - cxa * sinphi;
				Complex cxcc = cxc * cosphi + cxd * sinphi;
				Complex cxdd = cxd * cosphi - cxc * sinphi;
				myS1[nd] = -conei * (cxfData->Cx21(nd) * cxbb + cxfData->Cx22(nd) * cxdd);
				myS2[nd] = -conei * (cxfData->Cx11(nd) * cxaa + cxfData->Cx12(nd) * cxcc);
				myS3[nd] =  conei * (cxfData->Cx11(nd) * cxbb + cxfData->Cx12(nd) * cxdd);
				myS4[nd] =  conei * (cxfData->Cx21(nd) * cxaa + cxfData->Cx22(nd) * cxcc);
			}
			break;

		case PeriodicY:
		case PeriodicZ:
// JPBC = 1 or 2:
// ENSC(1-3,ND) = components of scattering unit vector in Lab Frame
// CXA,CXB = y,z components of incident polarization state 1 in Lab Frame
// CXC,CXD = y,z                                           2 in Lab Frame
			cxa = param->Cxe01_lf().data[1].conjg();
            cxb = param->Cxe01_lf().data[2].conjg();
            cxc = param->Cxe02_lf().data[1].conjg();
            cxd = param->Cxe02_lf().data[2].conjg();
			for(nd=0; nd<param->Nscat(); ++nd)
			{
				real cosphi =  em2_lf[nd].data[2];
				real sinphi = -em2_lf[nd].data[1];
				Complex cxaa = cxa * cosphi + cxb * sinphi;
				Complex cxbb = cxb * cosphi - cxa * sinphi;
				Complex cxcc = cxc * cosphi + cxd * sinphi;
				Complex cxdd = cxd * cosphi - cxc * sinphi;
				myS1[nd] = -conei * (cxfData->Cx21(nd) * cxbb + cxfData->Cx22(nd) * cxdd);
				myS2[nd] = -conei * (cxfData->Cx11(nd) * cxaa + cxfData->Cx12(nd) * cxcc);
				myS3[nd] =  conei * (cxfData->Cx11(nd) * cxbb + cxfData->Cx12(nd) * cxdd);
				myS4[nd] =  conei * (cxfData->Cx21(nd) * cxaa + cxfData->Cx22(nd) * cxcc);
			}
			break;

		case PeriodicBoth:
//fprintf(stderr, "Getmueller %p\n", cxfData);
//
// JPBC = 3
//     Target is periodic in y and z directions in Lab Frame Allowed scattering directions are identified by scattering order (M,N).  
//     For each allowed (M,N) there is one forward (transmitted) direction and one backward (reflected) direction.
//     Following convention set in subroutine PBCSCAVEC for JPBC=3, it is assumed that NSCAT is even, with the first NSCAT/2 directions 
//     corresponding to transmission, and the remaining NSCAT/2 directions corresponding to reflection
//
// ENSC(1-3,ND) = components of scattering unit vector in Lab Frame
// EM1(1-3,ND) = components of scattering polarization 1 in Lab Frame
// EM2(1-3,ND) = components of scattering polarization 2 in Lab Frame
			cxa = param->Cxe01_lf().data[1].conjg();
			cxb = param->Cxe01_lf().data[2].conjg();
			cxc = param->Cxe02_lf().data[1].conjg();
			cxd = param->Cxe02_lf().data[2].conjg();
			for(nd=0; nd<param->Nscat(); ++nd)
			{
				real cosphi =  em2_lf[nd].data[2];
				real sinphi = -em2_lf[nd].data[1];
				Complex cxaa = cxa * cosphi + cxb * sinphi;
				Complex cxbb = cxb * cosphi - cxa * sinphi;
				Complex cxcc = cxc * cosphi + cxd * sinphi;
				Complex cxdd = cxd * cosphi - cxc * sinphi;
				myS1[nd] = -conei * (cxfData->Cx21(nd) * cxbb + cxfData->Cx22(nd) * cxdd);
				myS2[nd] = -conei * (cxfData->Cx11(nd) * cxaa + cxfData->Cx12(nd) * cxcc);
				myS3[nd] =  conei * (cxfData->Cx11(nd) * cxbb + cxfData->Cx12(nd) * cxdd);
				myS4[nd] =  conei * (cxfData->Cx21(nd) * cxaa + cxfData->Cx22(nd) * cxcc);

//fprintf(stderr, "%lf %lf %lf %lf %lf %lf %lf %lf\n", cxfData->Cx11(nd).re, cxfData->Cx11(nd).im, cxfData->Cx12(nd).re, cxfData->Cx12(nd).im, 
//		 cxfData->Cx21(nd).re, cxfData->Cx21(nd).im, cxfData->Cx22(nd).re, cxfData->Cx22(nd).im);
//fprintf(stderr, "%lf %lf %lf %lf %lf %lf %lf %lf\n", myS1[nd].re, myS1[nd].im, myS2[nd].re, myS2[nd].im, myS3[nd].re, myS3[nd].im, myS4[nd].re, myS4[nd].im);

			}
			break;

		default:
			break;
		}
//
// if JPBC=1,2,3 (target periodic in y, z, or y and z directions):
// Check for special case: JPBC=3 (2-d target) and forward scattering wit
// ORDERM=ORDERN=0.  For this case we need to add incident wave to
// S_1 and S_2.
// Note that S_1 and S_2 are defined such that for M=N=0 we change
// iS_1 -> iS_1 +1  and iS_2 -> iS_2 +1
//  S_1 ->  S_1 -i       S_2 ->  S_2 -i
// note that CXS1 and CXS2 have yet to be multiplied by factor
// 2*pi/(ak1*aksr(1,nd)*pyddx*pzddx)
		if (jpbc == PeriodicBoth)
		{
			for(nd=0; nd<param->Nscat() / 2; ++nd)
			{
				if ((nint_(orderm[nd]) == 0) && (nint_(ordern[nd]) == 0))
				{
// Zero-degree forward-scattering need to add incident wave to radiated wave
// 2012.04.26 (BTD) change
// !                  CXFAC0=AKS_TF(1,ND)*AK1*PYDDX*PZDDX/(2._WP*PI)
					cxfac0 = Complex(Fabs(aks_tf[nd].data[0]) * ak1 * pyddx * pzddx / TwoPi, (real)0.);
// ! 2012.04.25 (BTD): for reasons not understood, incident wave added for S2 has different sign from S1
// !                   ... does this indicate a sign mistake in calculation of either S1 or S2?
// !                   ... is it possible that there is an error in eq. 68 of Draine & Flatau 2008?. Something like this could
// 							arise from a 180deg shift in directions between
// !                       the incident and scattered "parallel" basis vectors, or between the incident and scattered "perpendicular" basis vectors (should check for this in output files!)
// !                   This appears to indicate that for JPBC=3 the
// !                   sign of either CXS1 or CXS2 is incorrect
// !                   The Mueller matrix elements S11, S22, S12, and S21 appear
// !                   to be correct, but they would not be affected by
// !                   a sign error in either CXS1 or CXS2
// !                   Elements S13,S14,S23,S24,S31,S32,S33,S34,S41,S42,S43,S44
// !                   would be sensitive to a sign error, and it would be
// !                   a good idea to find a way to test them.
					myS1[nd] -= cxfac0;
					myS2[nd] += cxfac0;

//fprintf(stderr, "%d %lf %lf\n", nd, cxfac0.re, cxfac0.im);

				}
			}
		}
//
// Calculation of scattering amplitudes CXS1, CXS2, CXS3, CXS4 is complete. 
// Now compute Mueller matrix elements for this particular target orientation:

		int index;
		for(nd=0; nd<param->Nscat(); ++nd)
		{
			index = 16 * nd;
			sm[index +  0] = half_ * (myS1[nd] * myS1[nd].conjg() + myS2[nd] * myS2[nd].conjg() + myS3[nd] * myS3[nd].conjg() + myS4[nd] * myS4[nd].conjg()).re;
			sm[index +  1] = half_ * (myS2[nd] * myS2[nd].conjg() - myS1[nd] * myS1[nd].conjg() + myS4[nd] * myS4[nd].conjg() - myS3[nd] * myS3[nd].conjg()).re;
			sm[index +  2] =         (myS2[nd] * myS3[nd].conjg() + myS1[nd] * myS4[nd]).re;
			sm[index +  3] =         (myS2[nd] * myS3[nd].conjg() - myS1[nd] * myS4[nd]).im;
			sm[index +  4] = half_ * (myS2[nd] * myS2[nd].conjg() - myS1[nd] * myS1[nd].conjg() - myS4[nd] * myS4[nd].conjg() + myS3[nd] * myS3[nd].conjg()).re;
			sm[index +  5] = half_ * (myS2[nd] * myS2[nd].conjg() + myS1[nd] * myS1[nd].conjg() - myS4[nd] * myS4[nd].conjg() - myS3[nd] * myS3[nd].conjg()).re;
			sm[index +  6] =         (myS2[nd] * myS3[nd].conjg() - myS1[nd] * myS4[nd].conjg()).re;
			sm[index +  7] =         (myS2[nd] * myS3[nd].conjg() + myS1[nd] * myS4[nd].conjg()).im;
			sm[index +  8] =         (myS2[nd] * myS4[nd].conjg() + myS1[nd] * myS3[nd].conjg()).re;
			sm[index +  9] =         (myS2[nd] * myS4[nd].conjg() - myS1[nd] * myS3[nd].conjg()).re;
			sm[index + 10] =         (myS1[nd] * myS2[nd].conjg() + myS3[nd] * myS4[nd].conjg()).re;
			sm[index + 11] =         (myS2[nd] * myS1[nd].conjg() + myS4[nd] * myS3[nd].conjg()).im;
			sm[index + 12] =         (myS4[nd] * myS2[nd].conjg() + myS1[nd] * myS3[nd].conjg()).im;
			sm[index + 13] =         (myS4[nd] * myS2[nd].conjg() - myS1[nd] * myS3[nd].conjg()).im;
			sm[index + 14] =         (myS1[nd] * myS2[nd].conjg() - myS3[nd] * myS4[nd].conjg()).im;
			sm[index + 15] =         (myS1[nd] * myS2[nd].conjg() - myS3[nd] * myS4[nd].conjg()).re;
		}
//
		real fac0 = (real)0.;
		if (jpbc != PeriodicNo)
		{
			switch(jpbc)
			{
			case PeriodicY:							// Prepare to calculate S^{(1d)} for target with 1-d periodicity in y
				fac0 = TwoPi / (ak1 * pyddx * pyddx);
				break;

			case PeriodicZ:							// Prepare to calculate S^{(1d)} for target with 1-d periodicity in z
				fac0 = TwoPi / (ak1 * pzddx * pzddx);
				break;

			case PeriodicBoth:						// Prepare to calculate S^{(2d)} for target with periodicity in y and z
				fac0 = TwoPi / (ak1 * pyddx * pzddx);
				fac0 = fac0 * fac0;
				break;

			default:
				break;
			}
//
// k*sin(alpha_s)=component of k_s perpendicular to target
			for(nd=0; nd<param->Nscat(); ++nd)
			{
				switch(jpbc)
				{
				case PeriodicY:
					fac = fac0 / Sqrt(aks_tf[nd].data[0] * aks_tf[nd].data[0] + aks_tf[nd].data[2] * aks_tf[nd].data[2]);
					break;

				case PeriodicZ:
					fac = fac0 / Sqrt(aks_tf[nd].data[0] * aks_tf[nd].data[0] + aks_tf[nd].data[1] * aks_tf[nd].data[1]);
					break;

				case PeriodicBoth:
					fac = fac0 / (aks_tf[nd].data[0] * aks_tf[nd].data[0]);
					break;

				default:
					Errmsg("Fatal", "Getmueller", " invalid jpbc");
					break;
				}
               
				index = 16 * nd;
				for(k=0; k<16; ++k)
				{
					sm[index + k] *= fac;
				}
			}
		}
//
// Augment Mueller matrix for orientational average
		index = 16 * param->Nscat();
		for(nd=0; nd<index; ++nd)
		{
			smori[nd] += sm[nd] * wg;
		}
		delete myData;
	}
	else
	{
// 
// When IORTH=1, cannot compute full Mueller matrix. Compute three scattering properties:
		for(nd=0; nd<param->Nscat(); ++nd)
		{
			s1111[nd]  += (cxfData->Cx11(nd).conjg() * cxfData->Cx11(nd)).re * wg;
			s2121[nd]  += (cxfData->Cx21(nd).conjg() * cxfData->Cx21(nd)).re * wg;
			cx1121[nd] += (cxfData->Cx11(nd).conjg() * cxfData->Cx21(nd)) * wg;
		}
	}
// 
// Now augment sums of QABS,QBKSCA,QEXT,QSCA,G*QSCA over incident directions and polarizations.
	for(int jo=0; jo<param->Iorth(); ++jo)
	{
		sumPackage.Qext().Qsum1()[jo]  += sumPackage.Qext().Q()[jo] * wg;
		sumPackage.Qabs().Qsum1()[jo]  += sumPackage.Qabs().Q()[jo] * wg;
		sumPackage.Qbksca().Qsum1()[jo] += sumPackage.Qbksca().Q()[jo] * wg;
		sumPackage.Qpha().Qsum1()[jo] += sumPackage.Qpha().Q()[jo] * wg;
		sumPackage.Qsca().Qsum1()[jo] += sumPackage.Qsca().Q()[jo] * wg;
		sumPackage.Qscag2().Qsum1()[jo] += sumPackage.Qscag2().Q()[jo] * wg;
        sumPackage.Qscag().Qsum1()[jo] += sumPackage.Qscag().Q()[jo] * wg;
		if (DDscatParameters::GetInstance()->Cmdtrq() == TorqMethod_DOTORQ)
		{
			sumPackage.Qtrqab().Qsum1()[jo] += sumPackage.Qtrqab().Q()[jo] * wg;
			sumPackage.Qtrqsc().Qsum1()[jo] += sumPackage.Qtrqsc().Q()[jo] * wg;
		}
	}
}
