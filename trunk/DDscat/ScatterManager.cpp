#include "StdAfx.h"

#include "ScatterManager.h"
#include "Timeit.h"
#include "DDscatCommons.h"
#include "Functions.h"
#include "DDscatMain.h"
#include "Vect3.h"

/* **
class Ttt
{
protected:
	vector<real> vect;
	vector<string> vvvv;

public:
	Ttt() { Clear(); Fix("start", 0, 0); }
	~Ttt() { Fix("finis", 0, 0); Print(); }

	void Clear()
	{
		vect.clear();
		vvvv.clear();
	}
	void Print()
	{
		for(unsigned int i=0; i<vect.size(); ++i)
		{
			cerr << "Time #" << i << " = " << vect.at(i) << "  " << vvvv.at(i);
			if (i < vect.size()-1)
				cerr << "   delta = " << vect.at(i+1) - vect.at(i);
			cerr << endl;
		}
	}
	void Fix(char *label, int n, int m)
	{
		char Buf[64];
		vect.push_back(Cpu_time());
		sprintf(Buf, "%s %d %d", label, n, m);
		vvvv.push_back(Buf);
	}
};
** */

ScatterManager::ScatterManager(void)
{
	currentTarget = NULL;
	Init();
}

ScatterManager::ScatterManager(AbstractTarget *target)
{
	currentTarget = target;
	Init();
}	

ScatterManager::~ScatterManager(void)
{
	CleanDelete2(scrrs1);
	CleanDelete2(scrrs2);
	CleanDelete2(cscr1);
	CleanDelete2(cscr2);
	CleanDelete2(cscr3);
	CleanDelete2(cscr4);
}

void ScatterManager::Init()
{
	if (currentTarget)
	{
		int nat = currentTarget->Nat();
		scrrs1 = new real[3 * nat];
		scrrs2 = new real[nat];
		cscr1 = new Complex[nat];
		cscr2 = new Complex[nat];
		cscr3 = new Complex[nat];
		cscr4 = new Complex[3 * nat];
	}
	else
	{
		scrrs1 = scrrs2 = NULL;
		cscr1 = cscr2 = cscr3 = cscr4 = NULL;
	}
	torqMethod = TorqMethod_NOTORQ;
}

void ScatterManager::Scat(Vect3<real> &ak_tf, Vect3<real> *aks_tf, Vect3<real> *em1_tf, Vect3<real> *em2_tf, real e02, real etasca,  
	real &cbksca, real &csca, Vect3<real> &cscag, real &cscag2, Vect3<real> &ctrqin, Vect3<real> &ctrqsc, DipoleData *theDipoleData,
	Vect3<Complex> &cxe01_tf, Complex *cxf1l, Complex *cxf2l, int myid, PeriodicBoundaryFlag jpbc, int &navg, int ndir)
{
/* **
SCAT computes energy radiated by array of oscillating dipoles and corresponding scattering properties if oscillations are in response to incident E field.

Note: *** All vectors are given in the same frame *** (e.g., the Target Frame)

Given:
     AK_TF(1-3) = (kx,ky,kz)*d = where (kx,ky,kz)=incident k vector in TF d = (dx*dy*dz)**(1/3) = effective lattice spacing
     AKS_TF(3,MXSCA)=scattered k vectors in TF
     DX(1-3) = (dx,dy,dz)/d where dx,dy,dz=x,y,z lattice spacing
     CMDTRQ = 'DOTORQ' to calculate torques = 'NOTORQ' to skip calculation of torques
     EM1_TF(3,MXSCA)=pol.vector 1 in TF for each scattering direction
     EM2_TF(3,MXSCA)=pol.vector 2 in TF for each scattering direction
     E02 = |E_0|^2 , where E_0 = incident Complex E field
     NAT = number of dipoles
     NAT3 = 3*NAT
     IXYZ(NAT,1-3) = [x-X0(1)]/d, [y-X0(2)]/d, [z-X0(3)]/d (integers) for dipoles 1-NAT
     X0(1-3) = location/d in TF of lattice site with IXYZ=(0,0,0)
     CXE_TF(1-NAT,1-3) = components in TF of E field at each dipole at t=0
     CXE01_TF(1-3) = (Complex) reference polarization vector in TF
     CXP_TF(1-NAT,1-3) = components in TF of polarization of each dipole at t=0
     NDIR = number of directions at which scattering matrix elements are to be computed
     JPBC = 0 for finite target
            1 for target periodic in y_TF direction
            2 for target periodic in z_TF direction
            3 for target periodic in y_TF and z_TF directions
     MYID = id of this thread 

 Returns:
     CBKSCA   =differential scattering cross section for theta=pi (lattice units)
     CSCA     = scattering cross section (lattice units)
     CSCAG(1) = CSCA*<cos(theta)> , where theta=scattering angle
     CSCAG(2) = CSCA*<sin(theta)cos(phi)>
     CSCAG(3) = CSCA*<sin(theta)sin(phi)>
     CSCAG2   = CSCA*<cos^2(theta)>
     NAVG     = number of directions used for calculating CSCAG(1-3) and CSCAG2
    CTRQIN(1-3)=cross section for torque on grain due to incident E and B fields, for torque in directions x,y,z assuming incident
               radiation is along x direction, and "reference pol state" is in y direction.
   CTRQSC(1-3)=cross section for torque on grain due to scattered E and B fields,
               for torque in directions x,y,z assuming incident radiation is along x direction,
               and "reference pol state" is in y direction Here "torque cross section" is defined so that
               torque(1-3)=(momentum flux)*(1/k)*C_torque(1-3) Total torque cross section=CTRQIN+CTRQSC
   CXF1L(NDIR)=f_{1L} for direction NDIR
   CXF2L(NDIR)=f_{2L} for direction NDIR
               where f_{1L} connects incident polarization L to outgoing theta polarization (E in scattering plane)
               f_{2L} connects incident polarization L to outgoing phi polarization (E perpendicular to scattering plane)
               This is same f_{ml} notation as used by Draine (1988; Ap.J., 333, 848).

 Notes:
     AKS0(1-3)=k_s*d=scattered k vector (lattice units)
     AKSN(1-3)=k_s/|k_s|=unit vector in scattering direction
 Scratch vectors introduced to permit vectorization:
     SCRRS1 = real vector of length.GE.3*NAT
     SCRRS2 = real vector of length.GE.NAT
     CXSCR1 = Complex scratch array of length 3*MXNAT
     CXSCR2 = Complex scratch array of length 3*MXNAT
     CXSCR3 = Complex scratch array of length 3*MXNAT
     CXSCR4 = Complex scratch array of length 3*MXNAT

 B.T.Draine, Princeton Univ. Obs., 87.01.04
 History:
 Fortran versions history removed

 end history
 Copyright (C) 1993,1994,1995,1997,1998,2003,2007,2008,2011,2012
               B.T. Draine and P.J. Flatau
 This code is covered by the GNU General Public License.
** */

//	Ttt *ttt = new Ttt;

	const real zero_ = (real)0.;
	const real half_ = (real)0.5;
	const real threeFourth = (real)0.75;

	Complex *cxe_tf = theDipoleData->Cxe_tf();
	Complex *cxp_tf = theDipoleData->Cxxi(); 
//fprintf(stderr, " --- DDD\n");
//theDipoleData->Debug(stderr);

	int curTargetNat0 = currentTarget->Nat0();

	real ak2 = ak_tf.ModSquared();
	real akk = Sqrt(ak2);
	real ak3 = akk * ak2;

// =======================================================================
//
//         Calculate scattering and possibly torque for finite target
//
//    Method used to choose scattering directions for calculation of radiation force and torque:
//
//    Define function
//                        A*(theta/theta0)
//       s(theta)=theta + ----------------
//                        1 + theta/theta0
//
//    We choose angles theta to be uniformly-distributed in s. 
//    With A > 0 this gives increased resolution at small scattering angles, particularly for theta < theta0
//
//    We take
//       A=1
//       theta0=2*pi/(1.+x)
//
//    Parameter etasca determines the numbers of angles, chosen so that the largest interval is
//                                 pi/2
//    [Delta theta]_max = etasca * ----
//                                 3+x
//
//    A reasonable choice is etasca = 1
//
//    With etasca=1:
//         x=0 --> [Delta theta]_max = 30 deg
//           1                         22.5 deg
//           10                         6.92
//           15                         5 deg
//
//    With this requirement for [Delta theta]_max, the number of theta values (including 0 and pi) becomes
//
//                 2*(3+x)       1 + A/(pi+theta0)
//    NTHETA = 1 + ------- * --------------------------
//                 etasca    1 + A*theta0/(pi+theta0)^2
//
//    Thus, with A=1 and etasca=1
//    x=0: theta0 = 2*pi   NTHETA=7.19818 -> 7
//      1            pi           9.48970 -> 9
//      2          2*pi/3         12.0646 -> 12
//     10          2*pi/11        32.6897 -> 33
//
// etasca = 1 appears to work well.
// Use smaller value of ETASCA for improved accuracy
//
// Incident k vector AK_TF defines one axis (for spherical integration)
// Use (real part of) reference polarization vector CXE01_TF to define second coordinate axis.
// Construct third coordinate axis by taking cross product of first two axes.
// AKS1, AKS2 = vectors of length equal to AK_TF lying along second and third coordinate axes.

	const Complex conei((real)0., (real)1.);
	Complex cxes[3], cxes_l[3], cxbs[3];
	Vect3<real> aks0, aks1, aks2, aksn;
	real xyzcm[3], ctrqsctf[3];
	int j, jj, k, icosth, iphi, nphi;

//	ttt->Fix("a", 0, 0);
	navg = 0;
	if (jpbc == PeriodicNo)
	{
		aks1.Set(cxe01_tf.data[0].re, cxe01_tf.data[1].re, cxe01_tf.data[2].re);
		real einc = aks1.Mod();								// EINC = magnitude of (real part of) reference electric field
		if (einc <= zero_)
			Errmsg("Fatal", "Scat", "cxe01_tf is pure imaginary!");
		real rrr = akk / einc;
		aks1 *= rrr;
		aks2 = (ak_tf * aks1) / akk;
//
		real afac = (real)1.;
		real x = akk * Pow(threeFourth * curTargetNat0 / Pi, (real)(1./3.));
		real theta0 = TwoPi / ((real)1. + x);
		int ntheta = 1 + nint_(((real)2. * ((real)3. + x) / etasca) * ((real)1. + afac/(Pi + theta0)) / ((real)1. + afac*theta0/((Pi + theta0)*(Pi + theta0))));
//
// Prepare for first interval in theta
		real theta = zero_;
		real thetau = zero_;

		csca = zero_;
		cscag.Clear();
		ctrqsctf[0] = ctrqsctf[1] = ctrqsctf[2] = zero_;
		cscag2 = zero_;

//		ttt->Fix("b", 0, 0);
		if (torqMethod == TorqMethod_DOTORQ)
		{
// Additional quantities required for torque computations
// XYZCM(1-3) = coordinates of grain centroid
// SCRRS1(J,K)= r_j = coordinates of dipole J relative to centroid.
// CXSCR3(J) = r_j dot p_j
// CXSCR4(J,1-3) = p_j cross r_j

#ifdef openmp
!$omp parallel 
#endif

			real temp[3];
			xyzcm[0] = zero_;
			xyzcm[1] = zero_;
			xyzcm[2] = zero_;

#ifdef openmp
!$omp single
#endif

#ifdef openmp
!$OMP END SINGLE
!$OMP DO                  &
!$OMP    PRIVATE(J)       &
!$OMP    REDUCTION(+:XYZCM)
#endif
			for(j=jj=0; j<curTargetNat0; ++j, jj+=3)
			{
				currentTarget->DipoleVector(j, temp);
				xyzcm[0] += temp[0];
				xyzcm[1] += temp[1];
				xyzcm[2] += temp[2];
				memcpy(scrrs1+jj, temp, 3*sizeof(real));
			}
#ifdef openmp
!$omp end do
#endif
			xyzcm[0] /= (real)curTargetNat0;
			xyzcm[1] /= (real)curTargetNat0;
			xyzcm[2] /= (real)curTargetNat0;
//
#ifdef openmp
!$OMP DO          &
!$OMP&   PRIVATE(J)
#endif
			for(j=jj=0; j<curTargetNat0; ++j, jj+=3)
			{
				for(k=0; k<3; ++k)
				{
					scrrs1[jj + k] -= xyzcm[k];
				}
			}
#ifdef openmp
!$omp end do
#endif
#ifdef openmp
!$OMP DO            &
!$OMP&   PRIVATE(J,K)
#endif
			for(j=jj=0; j<curTargetNat0; ++j, jj+=3)
			{
				cscr3[j] = Scalar(cxp_tf+jj, scrrs1+jj);	// (cxp_tf[jj] * scrrs1[jj] + cxp_tf[jj+1] * scrrs1[jj+1] + cxp_tf[jj+2] * scrrs1[jj+2]);
			}
#ifdef openmp
!$OMP END DO
!$OMP DO          &
!$OMP&   PRIVATE(J)
#endif
			for(j=jj=0; j<curTargetNat0; ++j, jj+=3)
			{
				Vector(cxp_tf+jj, scrrs1+jj, cscr4+jj);
			}
#ifdef openmp
!$omp end do
!$omp end parallel
#endif
		}
		else
		{
#ifdef openmp
!$omp parallel
#endif
//#ifdef openmp
//!$OMP DO          &
//!$OMP&   PRIVATE(J)
//#endif
			for(j=jj=0; j<curTargetNat0; ++j, jj+=3)
			{
				currentTarget->DipoleVector(j, scrrs1 + jj);
			}
//#ifdef openmp
//!$omp end do
//#endif

#ifdef openmp
!$omp end parallel
#endif
		}																	//  end IF(CMDTRQ=='DOTORQ')
//		ttt->Fix("c", 0, 0);
//
// Proceed to sum over scattering angles:
		Complex cxscl1, cxscl2, cxscl3;
		real es2;
		for(icosth=0; icosth<ntheta; ++icosth)
		{
//			ttt->Fix("d", icosth, 0);
// Set the scattering angle THETA and calculate THETAU = scattering angle for next value of ICOSTH
			real thetal = theta;
			theta = thetau;
			if (icosth < ntheta - 1)
			{
// THETAU will be the *next* value of theta: SI = value of S for *next* value of ICOSTH
				real si = (icosth + 1) * (Pi + afac * Pi / (theta0 + Pi)) / (ntheta - 1);
				real term = si - afac - theta0;
				thetau = half_ * (term + Sqrt(term*term + (real)4. * theta0*si));
			}
			else
			{
				thetau = Pi;
			}

			real domega = TwoPi * (Cos(half_ * (thetal + theta)) - Cos(half_ * (theta + thetau)));

			real costh = Cos(theta);
			real sinth = Sin(theta);

			if (icosth == 0 || icosth == ntheta - 1)
			{
				nphi = 1;
			}
            else
			{
				nphi = nint_(FourPi * sinth / (thetau - thetal));
				nphi = max(nphi, 3);
			}
//
// Evaluate weight factor W
			real w = domega / (real)nphi;
			real dphi = TwoPi / (real)nphi;
			navg += nphi;
			for(iphi=0; iphi<nphi; ++iphi)
			{
//				ttt->Fix("e", icosth, iphi);
//
// Offset phi values by 0.5*DPHI each time
				real phi = dphi * ((real)(iphi + 1) - half_ * ((icosth + 1) % 2));
				real sinphi = Sin(phi);
				real cosphi = Cos(phi);
//
// Compute scattered k vector AKS0(1-3)
// Initialize CXES(1-3) = scattered electric field
				aks0 = ak_tf * costh + (aks1 * cosphi + aks2 * sinphi) * sinth;
				memset(cxes, 0, 3*sizeof(Complex));
//				cxes[0] = cxes[1] = cxes[2] = czero;
//
// Compute normalized scattering vector = nhat
				aksn = aks0 / akk;
//
//Now compute CXES(1-3) = sum_j[p_j-nhat(nhat dot p_j)]exp[-ik_s dot x_j
//      =(r/k*k)*exp(-ikr+iwt)*E(k_s,r,t)
//      where E(k_s,r,t) is Complex electric field vector propagating
//      with wave vector k_s, at distance r from origin
// Notation:
//    If CMDTRQ='DOTORQ', then r_j = position relative to centroid XYZCM
//    If CMDTRQ='NOTORQ', then r_j = original lattice coords, and XYZCM=0
// Note that since we are not concerned with overall phase shifts for
// evaluation of forces and torques, we do NOT need to worry about change
// in phase resulting from using two different definitions of r_j when
// computing phase factors exp[-i k_s dot r_j]

#ifdef openmp
!$OMP PARALLEL                  & 
!$OMP&   PRIVATE(CXES_L,CXSCL2_L)
#endif
				memset(cxes_l, 0, 3*sizeof(Complex));
				//cxes_l[0] = cxes_l[1] = cxes_l[2] = czero;							// Art
				const real *aksnData = aksn.data;
#ifdef openmp
!$OMP DO          &
!$OMP& PRIVATE(J,K)
#endif
//				ttt->Fix("f", icosth, iphi);
				for(j=jj=0; j<curTargetNat0; ++j, jj+=3)
				{
//
// CXSCR1(J) = nhat dot p(j)
					cscr1[j] = aksn.Scalar(cxp_tf + jj);
//
// SCRRS2(J) = nhat dot r(j)
					scrrs2[j] = aksn.Scalar(scrrs1 + jj);
//
// CXSCR2(J) = exp(-i*k_s dot r_j) factor for atom j
					cscr2[j] = Complex(Cos(akk * scrrs2[j]), -Sin(akk * scrrs2[j]));
//
// CXES(1-3) = sum_j[p_j-nhat(nhat dot p_j)]exp[-ik_s dot x_j]
					cxes_l[0] += (cxp_tf[jj  ] - cscr1[j] * aksnData[0]) * cscr2[j];
					cxes_l[1] += (cxp_tf[jj+1] - cscr1[j] * aksnData[1]) * cscr2[j];
					cxes_l[2] += (cxp_tf[jj+2] - cscr1[j] * aksnData[2]) * cscr2[j];
				}
				cxes[0] += cxes_l[0];
				cxes[1] += cxes_l[1];
				cxes[2] += cxes_l[2];
//
//				ttt->Fix("g", icosth, iphi);
				if (torqMethod == TorqMethod_DOTORQ)
				{
// Additional terms required to compute torque on grain
// Have already computed
// CXES(1-3)=sum_j[p_j - k_s*(k_s dot p_j)/|k|^2]exp[-ik_s dot x_j]
//             where p_j = polarization of dipole j, x_j = position of dipole j rel. to CM, k_s = scattered k vector
// CXBS = -nhat cross CXES  = -(r/k*k)*exp(-ikr+iwt)*B(k_s,r,t)
					aksn.Vector(cxes, cxbs);
//
// CXSCL1 = sum_j p_j dot (r_j cross nhat)*exp(-i*k_s dot r_j)
//        = nhat dot sum_j (p_j cross r_j)*exp(-i*k_s dot r_j)
// recall that have previously evaluated
// CXSCR2(J) = exp(-i*k_s dot r_j) factor for atom j
// CXSCR4(J,K)=p_j cross r_j
#ifdef openmp
!$omp single
#endif
					cxscl1.clear();
#ifdef openmp
!$omp end single
#endif

//!Art                  DO K=1,3
//!Art                     CXTRM(K)=CXZERO
//!Art                  ENDDO

#ifdef openmp
!$omp single
#endif
					Complex cxtrm_k[3];
					cxtrm_k[0].clear();
					cxtrm_k[1].clear();
					cxtrm_k[2].clear();
#ifdef openmp
!$omp end single
!$OMP DO                    &     
!$OMP&   PRIVATE(J)         &
!$OMP&   REDUCTION(+:CXTRM_K)
#endif
					for(j=jj=0; j<curTargetNat0; ++j, jj+=3)
					{
// !Art                       CXTRM(K)=CXTRM(K)+CXSCR2(J)*CXSCR4(J,K)
						cxtrm_k[0] += (cscr2[j] * cscr4[jj  ]);
						cxtrm_k[1] += (cscr2[j] * cscr4[jj+1]);
						cxtrm_k[2] += (cscr2[j] * cscr4[jj+2]);
					}
#ifdef openmp
!$omp end do
!$omp atomic
#endif
// !Art                    CXSCL1=CXSCL1+AKSN(K)*CXTRM(K)
					cxscl1 += aksn.Scalar(cxtrm_k);

//! CXSCL2 = sum_j[(r_j dot p_j) - (nhat dot p_j)*((nhat dot r_j) + 2i/k_s)] * exp(-i*k_s dot r_j)
//! Have previously computed
//! CXSCR1(J) = nhat dot p(j)
//! CXSCR2(J) = exp(-i*k_s dot r_j) factor for atom j
//! CXSCR3(J) = r_j dot p_j
//! SCRRS2(J) = nhat dot r(j)

#ifdef openmp
!$omp single
#endif
					cxscl2.clear();
					cxscl3 = Complex(zero_, (real)2. / akk);
#ifdef openmp
!$OMP END SINGLE
!$OMP DO                    &
!&OMP    PRIVATE(J)         &       
!$OMP&   REDUCTION (+:CXSCL2)
#endif
					for(j=0; j<curTargetNat0; ++j)
					{
						cxscl2 += (cscr2[j] * (cscr3[j] - cscr1[j] * (cxscl3 + scrrs2[j])));
					}
#ifdef openmp
!$omp end do
#endif

// Note following notational correspondence:
//    present    Draine & Weingartner (1996)
//    -------    ---------------------------
//     CXES   =        V_E
//     CXBS   =        V_B
//    CXSCL1  =        S_B
//    CXSCL2  =        S_E
//
				}
#ifdef openmp
!$omp end parallel
#endif
//
// Have completed computation of vector CXES(1-3) and terms to calculate
// angular momentum flux for scattering direction (THETA,PHI).
// CTRQSCTF(K) is "torque cross section" due to "scattered radiation"
// resulting in torque in lattice direction k (i.e., in "target
// frame".  After CTRQSCTF has been evaluated, we will compute CTRQSC in
// the "scattering frame", where radiation is incident along the x direct
// and "reference polarization state" is along the y direction.
				es2 = zero_;
				for(k=0; k<3; ++k)
				{
					es2 += (cxes[k] * cxes[k].conjg()).re;
				}
				real tmp = w * es2;
				csca += tmp;
				cscag += Vect3 <real>(costh, sinth * cosphi, sinth * sinphi) * tmp;
				cscag2 += (costh * costh * tmp);
				if (torqMethod == TorqMethod_DOTORQ)
				{
					for(k=0; k<3; ++k)
					{
						ctrqsctf[k] -= (w * (cxbs[k].conjg() * cxscl2 + cxes[k].conjg() * cxscl1).re);
					}
				}
			}															// End DO IPHI=1,NPHI
//
// Save backscattering cross section:
			if (icosth == ntheta - 1)
				cbksca = ak2 * ak2 * es2 / e02;
		}																// End DO ICOSTH=1,NTHETA
//
		csca = ak2 * ak2 * csca / e02;
// 
// Convert CSCAG,CSCAG2, and CTRQSCTF to cross sections measured in lattice units:
		cscag *= (ak2 * ak2 / e02);
		cscag2 *= (ak2 * ak2 / e02);
//
//		ttt->Fix("k", 0, 0);
		if (torqMethod == TorqMethod_DOTORQ) 
		{
//
// Compute the two contributions to the torque on the target:
// (1) torque due to E and B fields of incident radiation (CTRQIN)
// (2) torque due to E and B fields of scattered radiation (CTRQSC)
// These are first computed in the target frame, then transformed to
// the lab frame (frame where incident radiation propagates along
// the x axis, and "reference pol. state 1" is along the y axis).
			ctrqsctf[0] *= (ak2 * ak2 / e02);
			ctrqsctf[1] *= (ak2 * ak2 / e02);
			ctrqsctf[2] *= (ak2 * ak2 / e02);

// CTRQSCTF(1-3) is the vector torque cross section/|k| measured in the
// target frame.  We now compute CTRQSC(1-3)=the vector torque cross
// section in the "scattering frame", where the radiation is incident alo
// the x axis, and the "reference polarization state" (at t=0) is
// linearly polarized along the y axis.
// Note that the following computation introduces an additional factor of |k|.
			ctrqsc.Set(ak_tf.Scalar(ctrqsctf), aks1.Scalar(ctrqsctf), aks2.Scalar(ctrqsctf));
//
// Now compute vector torque cross section due to incident radiation:
			ctrqin.Clear();
//
// Redefine scratch vector CXSCR1: CXSCR1(J)=conjg(p_j) dot E_0j at t=0
#ifdef openmp
!$OMP PARALLEL 
!$OMP DO            & 
!$OMP&   PRIVATE(J,K)
#endif
			for(j=jj=0; j<curTargetNat0; ++j, jj+=3)
			{
				cscr1[j] = (cxp_tf[jj].conjg() * cxe_tf[jj]) + (cxp_tf[jj+1].conjg() * cxe_tf[jj+1]) + (cxp_tf[jj+2].conjg() * cxe_tf[jj+2]);
			}
#ifdef openmp
!$omp end do
!$omp parallel
#endif
#ifdef openmp
!$OMP SINGLE
#endif

// ! 12.04.28 (BTD) I cannot see what is accomplished by this
// !                isolated OMP SINGLE section
			cxscl1.clear();
			cxscl2.clear();
			cxscl3.clear();
#ifdef openmp
!$omp end single
#endif
//
// Compute CXSCL1,CXSCL2,CXSCL3=components of the torque due to
// incident E and B fields in the target frame.
// 95.06.22 (BTD): added factor CXI in r cross k term appearing in evaluations of CXSCL1,CXSCL2,CXSCL2
// 95.07.21 (BTD): corrected sign mistake (+CXI -> -CXI)
// 95.07.24 (BTD): changed def. of CXSCR1 to agree with paper, changed calc. of CXSCLn to more closely mirror paper, and changed AKS0 to AK
#ifdef openmp
!$OMP PARALLEL 
!$OMP DO                                 &
!$OMP&   PRIVATE(J)                      &
!$OMP&   REDUCTION(+:CXSCL1,CXSCL2,CXSCL3)
#endif
			const real *aktfData = ak_tf.data;
			for(j=jj=0; j<curTargetNat0; ++j, jj+=3)
			{
				cxscl1 += (cxp_tf[jj+1].conjg() * cxe_tf[jj+2] - cxp_tf[jj+2].conjg() * cxe_tf[jj+1] - conei * (aktfData[1] * scrrs1[jj+2] - aktfData[2] * scrrs1[jj+1]) * cscr1[j]);
				cxscl2 += (cxp_tf[jj+2].conjg() * cxe_tf[jj  ] - cxp_tf[jj  ].conjg() * cxe_tf[jj+2] - conei * (aktfData[2] * scrrs1[jj  ] - aktfData[0] * scrrs1[jj+2]) * cscr1[j]);
				cxscl3 += (cxp_tf[jj  ].conjg() * cxe_tf[jj+1] - cxp_tf[jj+1].conjg() * cxe_tf[jj  ] - conei * (aktfData[0] * scrrs1[jj+1] - aktfData[1] * scrrs1[jj  ]) * cscr1[j]);
			}
#ifdef openmp
!$omp end do
!$omp end parallel
#endif
//
//Compute CTRQIN(1-3)=2*components of torque due to incident E and B fields (in lab frame), expressed as a cross section:
			ctrqin.Set(  aktfData[0] * cxscl1.re +  aktfData[1] * cxscl2.re +  aktfData[2] * cxscl3.re,
						aks1.data[0] * cxscl1.re + aks1.data[1] * cxscl2.re + aks1.data[2] * cxscl3.re,
						aks2.data[0] * cxscl1.re + aks2.data[1] * cxscl2.re + aks2.data[2] * cxscl3.re);
			ctrqin *= (FourPi / e02);
		}
//
//
		char cmsgnm[256];
		sprintf(cmsgnm, "%8d  scattering directions used to calculate <cos>, etc.", navg);
		Wrimsg("Scat", cmsgnm);

		if (ndir <= 0) 
			return;
	}															// ! end IF(JPBC==0)
//
//============== Following code is used to calculate ====================
//               CXF1L(J=1-NDIR) and CXF2L(J=1-NDIR)
//               for JPBC=0,1,2,3 (i.e., all cases)
//
// Calculate the 2 matrix elements FM1 for selected directions
// F1L(NDIR) is to the polarization state e1 (E in scattering plane)
// F2L(NDIR) is to the polarization state e2 (E perp. to scatt. plane)
//
// IMPORTANT NOTE: If the incident radiation is not linearly polarized, with zero imaginary component to the polarization vector, then the
// F1L and F2L computed here are NOT the "usual" f_ml, which describe scattering from linearly polarized states l to linearly polarized states m.
// If one wishes to reconstruct the "usual" f_ml even when computing with elliptically polarized light, then one must add additional
// code to subroutine SCAT to reconstruct the f_ml from the F1L and F2L computed here.
//
//	ttt->Fix("l", 0, 0);
//    fprintf(stderr, "%d\n", curTargetNat0);
//	fprintf(stderr, "Scat x0: %lf %lf %lf\n", currentTarget->X0().data[0], currentTarget->X0().data[1], currentTarget->X0().data[2]);
//  fprintf(stderr, "%lf %lf %lf\n", currentTarget->Dx().data[0], currentTarget->Dx().data[1], currentTarget->Dx().data[2]);
//	fprintf(stderr, "IXYZ\n");
//    for(j=0; j<32; ++j)
//	{
//		fprintf(stderr, "%d %d %d %d %lf %lf %lf %lf %lf %lf\n", 
//			j, currentTarget->Ixyz().Value(j, 0), currentTarget->Ixyz().Value(j, 1), currentTarget->Ixyz().Value(j, 2), 
//			cxp_tf[3*j].re, cxp_tf[3*j].im, cxp_tf[3*j+1].re, cxp_tf[3*j+1].im, cxp_tf[3*j+2].re, cxp_tf[3*j+2].im);
//	}

	Complex cxf1l_l, cxf2l_l;
	real term = ak3 / Sqrt(e02);
//fprintf(stderr, "Term = %lf\n", term);
//currentTarget->Debug(NULL);
	for(int nd=0; nd<ndir; ++nd)
	{
#ifdef openmp
!$OMP parallel 
!$omp single
#endif
		cxf1l_l.clear();
		cxf2l_l.clear();
#ifdef openmp
!$OMP END SINGLE
!$OMP DO                              &
!$OMP&   PRIVATE(K)                   &
!$OMP&   REDUCTION(+: CXF1L_L, CXF2L_L)
#endif
		for(j=jj=0; j<curTargetNat0; ++j, jj+=3)
		{
			real tmp = currentTarget->DipoleScalar(aks_tf[nd], j);	
			cscr1[j] = Complex(Cos(tmp), -Sin(tmp));
			cscr2[j] = em1_tf[nd].Scalar(cxp_tf + jj);
			cscr3[j] = em2_tf[nd].Scalar(cxp_tf + jj);
			cxf1l_l += cscr2[j] * cscr1[j];
			cxf2l_l += cscr3[j] * cscr1[j];
//int *dq = currentTarget->Ixyz().Row(j);
//fprintf(stderr, "Sca: %d %d %d %d %lf %lf %lf %lf %lf %lf\n", j, dq[0], dq[1], dq[2], cscr1[j].re, cscr1[j].im, cscr2[j].re, cscr2[j].im, cscr3[j].re, cscr3[j].im);
//fprintf(stderr, "Sca: %d %lf %lf %lf %lf %lf %lf\n", j, (*(cxp_tf + jj)).re, (*(cxp_tf + jj)).im, (*(cxp_tf + jj + 1)).re, (*(cxp_tf + jj + 1)).im, (*(cxp_tf + jj + 2)).re, (*(cxp_tf + jj + 2)).im);
		}
//fprintf(stderr, "%d %lf %lf %lf\n", nd, aks_tf[nd].data[0], aks_tf[nd].data[1], aks_tf[nd].data[2]);
//fprintf(stderr, "%d %lf %lf %lf\n", nd, em1_tf[nd].data[0], em1_tf[nd].data[1], em1_tf[nd].data[2]);
//fprintf(stderr, "%d %lf %lf %lf\n", nd, em2_tf[nd].data[0], em2_tf[nd].data[1], em2_tf[nd].data[2]);
//fprintf(stderr, "%d %lf %lf %lf %lf\n", nd, cxf1l_l.re, cxf1l_l.im, cxf2l_l.re, cxf2l_l.im);
#ifdef openmp
!$omp end do
!$omp end parallel
#endif

// !Art!$omp end parallel
//
//   Units: Dipole polarizations are in units of [E]*d**3, where d=lattice spacing, and [E]=dimensions of E field
//          Therefore, up to this point CXF1L, CXF2L are in units of [E]*d**3
//          Now multiply by AK3/SQRT(E02) (units of 1/([E]d**3)) to get dimensionless result.
		cxf1l[nd] = cxf1l_l * term;
		cxf2l[nd] = cxf2l_l * term;
//fprintf(stderr, "Cxf1l %d %lf %lf\n", nd, cxf1l[nd].re, cxf1l[nd].im);
//fprintf(stderr, "Cxf2l %d %lf %lf\n", nd, cxf2l[nd].re, cxf2l[nd].im);
	}
//	ttt->Fix("z", 0, 0);  // !--- end DO ND=1,NDIR
//	delete ttt;
}
