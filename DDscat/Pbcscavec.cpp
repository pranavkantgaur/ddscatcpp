#include "StdAfx.h"

#include "Definitions.h"
#include "DDscatMain.h"

void Pbcscavec(PeriodicBoundaryFlag jpbc, int nscat, real pyddx, real pzddx, Vect3<real> &a1, Vect3<real> &a2, real theta, real beta, 
	Vect3<real> &xlr, Vect3<real> &ylr, Vect3<real> &zlr, real ak1, Vect3<real> &en0r, Vect3<Complex> &cxe01r, 
	real *orderm, real *ordern, real *thetan, real *phin, Vect3<real> *aksr, Vect3<real> *em1r, Vect3<real> *em2r, Vect3<real> *ensc, Vect3<real> *em1, Vect3<real> *em2)
{
/* **
subroutine PBCSCAVEC
given:

    JPBC = 1 for target periodic in y_TF only (PZD=0)
           2                        z_TF only (PYD=0)
           3                        both y_TF and z_TF

    MXSCA = dimensioning information
    PYDDX = periodicity in y direction/d
    PZDDX = periodicity in z direction/d
    A1(3) = target axis a1 in TF (TF = Target Frame)    [**NOT USED**]
    A2(3) = target axis a2 in TF                        [**NOT USED**]
    THETA = rotation angle (rad) of a1 relative to x_LF
    BETA  = rotation angle (rad) of target around axis a1
    XLR(3)= unit vector xlab in TF = direction of incident radiation
    YLR(3)= unit vector ylab in TF
    ZLR(3)= unit vector zlab in TF

    AK1   = k*d in vacuo
    NSCAT = number of scattering directions
    EN0R(1-3) = unit vector in incident wave direction in TF
    CXE01R(1-3) = (Complex) incident pol.vector 1 in TF
    ORDERM(1-NSCAT)=diffraction order in y_TF direction if JPBC=1 or 3
                                      in z_TF direction if JPBC=2
    ORDERN(1-NSCAT)=rotation angle (rad) around scattering cone if JPBC=1 or 2
                   =diffraction order in z_TF direction if JPBC=3

returns:

    AKSR(3,1-NSCAT) = scattering vectors in TF
    EM1R(3,1-NSCAT) = scattered polarization vector parallel to scattering plane in TF
    EM2R(3,1-NSCAT) = scattered polarization vector perpendicular to scattering plane in TF
    ENSC(3,1-NSCAT) = scattering unit vectors in Lab Frame
    EM1(3,1-NSCAT)  = scattered polarization vector parallel to scattering plane in LF
    EM2(3,1-NSCAT)  = scattered polarization vector perpendicular to scattering plane in LF

if JPBC=1 or 2:
    THETAN(1-NSCAT) = alpha_s (rad) in TF
    PHIN(1-NSCAT)   = zeta_s (rad) in TF

where the scattering directions are determined as follows:

if JPBC=1:
   input ORDERM(J) = scattering order for y (int)
         ORDERN(J) = azimuthal angle zeta (radians) around y_TF 
                     where ORDERM=0 and zeta=0 for forward scattering
   computes THETAN(J) = angle between k_s and y_TF
            PHIN(J)   = azimuthal angle zeta (radians)

if JPBC=2:
   input ORDERM(J) = scattering order for z (int)
         ORDERN(J) = azimuthal angle zeta (radians) around z_TF 
                     where ORDERM=0 and zeta=0 for forward scattering
   computes THETAN(J) = angle between k_s and z_TF
            PHIN(J)   = azimuthal angle zeta (radians)

if JPBC=3:
   input ORDERM(J) = scattering order for y (int)
         ORDERN(J) = scattering order for z (int)
                     where ORDERM=0 and ORDERM=0 for forward scattering

    NSCAT is even: each scattering order (M,N) appears twice:
                   first for transmission, later for reflection
                   (this is established in REAPAR)
      J = 1         -> NSCAT/2 correspond to transmission
      J = NSCAT/2+1 -> NSCAT   correspond to reflection

    computes THETAN(J) = angle (rad) between k_s and x_TF
             PHIN(J)   = azimuthal angle (around x_TF)
                         k_sx=k_s*cos(thetan)
                         k_sy=k_s*sin(thetan)*cos(phin)
                         k_sz=k_s*sin(thetan)*sin(phin)

    in special case of zeroth-order transmission, the angle
    THETAN(J)=0, and the angle PHIN(J) is taken to be the same as for (0,0) reflection.
    In the special case of (0,0) normal incidence, where k_sy=k_sz=0, we set PHIN=0 

Fortran history records removed.

Copyright (C) 2004,2007,2008,2009,2010 B.T. Draine and P.J. Flatau
Copyright (C) 2012, C++ version, Choliy V.

This code is covered by the GNU General Public License.
** */

// for use in this module, define
// EM01R = unit polarization vector for incident polarization state 1 in TF
// AK2 = |kd|^2

	int j, k;
	Vect3<real> em01r;
	for(k=0; k<3; ++k)
	{
		em01r.data[k] = cxe01r.data[k].re;
	}
	real sum = em01r.Mod();
	em01r /= sum;
	real ak2 = ak1 * ak1;

	real akperp0, akperp;
	switch(jpbc)
	{
	case PeriodicY:
// JPBC=1 : target periodic in y_TF direction
//       -> input ORDERM = diffraction order in y_TF direction
//                ORDERN = azimuthal angle around y_TF axis
//                         with ORDERN=0 for forward scattering
		akperp0 = ak1 * Sqrt(en0r.data[0] * en0r.data[0] + en0r.data[2] * en0r.data[2]);
		for(j=0; j<nscat; ++j)
		{
			aksr[j].data[1] = ak1 * en0r.data[1] + orderm[j] * TwoPi / pyddx;
			if (Fabs(aksr[j].data[1]) <= ak1)
			{
				phin[j] = ordern[j];
				thetan[j] = Acos(aksr[j].data[1] / ak1);
				akperp = Sqrt(ak2 - aksr[j].data[1] * aksr[j].data[1]);
				aksr[j].data[2] = (akperp/akperp0) * ak1 * (en0r.data[2] * Cos(phin[j]) - en0r.data[0] * Sin(phin[j]));
				aksr[j].data[0] = (akperp/akperp0) * ak1 * (en0r.data[0] * Cos(phin[j]) + en0r.data[2] * Sin(phin[j]));
			}
			else
			{
				Errmsg("Fatal", "Pbcscavec", "invalid diffraction order");
			}
		}
		break;

	case PeriodicZ:
// JPBC=2 : target periodic in z direction
//       -> input ORDERM = diffraction order in z_TF direction
//                ORDERN = azimuthal angle around z_TF axis
//                         with ORDERN=0 for forward scattering
		akperp0 = ak1 * Sqrt(en0r.data[0]*en0r.data[0] + en0r.data[1]*en0r.data[1]);
		for(j=0; j<nscat; ++j)
		{
			aksr[j].data[2] = ak1 * en0r.data[2] + orderm[j] * TwoPi / pzddx;
			if (Fabs(aksr[j].data[2]) <= ak1)
			{
				phin[j] = ordern[j];
				thetan[j] = Acos(aksr[j].data[2]/ak1);
				akperp = Sqrt(ak2 - aksr[j].data[2] * aksr[j].data[2]);
				aksr[j].data[0] = (akperp/akperp0) * ak1 * (en0r.data[0] * Cos(phin[j]) - en0r.data[1] * Sin(phin[j]));
				aksr[j].data[1] = (akperp/akperp0) * ak1 * (en0r.data[1] * Cos(phin[j]) + en0r.data[0] * Sin(phin[j]));
			}
			else
			{
				Errmsg("Fatal", "Pbcscavec", "invalid diffraction order");
			}
		}
		break;

	case PeriodicBoth:
// JPBC=3 : infinite and periodic in y and z directions
//      -> input ORDERM = diffraction order in y_TF direction
//               ORDERN = diffraction order in z_TF direction
//         NSCAT is assumed to be an even number (this is set by subroutine REAPAR)
//         every (M,N) case appears twice: 
//            first for transmission (directions J = 1 -> NSCAT/2)
//            then for reflection (directions J = NSCAT/2+1 -> NSCAT)
//         EN0R(1-3) = components of unit vector for incident direction in TF
//         calculate
//             AKSR(1-3,J) = scattered vector in TF for J=1,NSCAT
//             THETAN(J)   = angle between k_s and x_TF
//             PHIN(J)     = azimuthal angle in TF
//             where  k_sx = k_s*cos(thetan)
//                    k_sy = k_s*sin(thetan)*cos(phin)
//                    k_sz = k_s*sin(thetan)*sin(phin)
//             Note: for special case k_sy=k_sz=0, we set phin=0
		real akperp2;
		int jt;
		for(jt=0; jt<nscat/2; ++jt)
		{
			int jr = nscat/2 + jt;
			aksr[jt].data[1] = ak1*en0r.data[1] + orderm[jt] * TwoPi / pyddx;
			aksr[jt].data[2] = ak1*en0r.data[2] + ordern[jt] * TwoPi / pzddx;
			aksr[jr].data[1] = aksr[jt].data[1];
			aksr[jr].data[2] = aksr[jt].data[2];
			akperp2 = aksr[jr].data[1] * aksr[jr].data[1] + aksr[jr].data[2] * aksr[jr].data[2];
			if (akperp2 <= ak2)
			{
				aksr[jt].data[0] = Sqrt(ak2-akperp2) * sign_(en0r.data[0]);
				aksr[jr].data[0] = -aksr[jt].data[0];
				thetan[jt] = Acos(aksr[jt].data[0] / ak1);
				thetan[jr] = Acos(aksr[jr].data[0] / ak1);
				akperp = Sqrt(akperp2);
				if (akperp > (real)0.)
				{
					real cosphi = aksr[jt].data[1] / akperp;
					real sinphi = aksr[jt].data[2] / akperp;
					if (sinphi < (real)0.)
						phin[jt] = TwoPi - Acos(cosphi);
					else
						phin[jt] = Acos(cosphi);
				}
				else
                    phin[jt] = (real)0.;
				phin[jr] = phin[jt];
			}
			else
				Errmsg("Fatal", "Pbcscavec", "invalid diffraction order");
		}								// ! end loop over JT
		break;

	default:
		Errmsg("Fatal", "Pbcscavec", " invalid value of JPBC");
		break;
	}
// scattering directions AKSR are now defined.
// For each transmission direction that is not either parallel or 
// antiparallel to the incident beam there is a scattering plane defined
// by vectors k_0 and k_s.
// Construct unit vector for scattered polarization
// parallel to plane (EM1R) and perpendicular to plane (EM2R)
// We adopt convention of Bohren & Huffman 1984 sec 3.2 for
// unit vectors parallel and perpendicular to the scattering plane:
//   em1 = ehat_{para,s} is in direction of increasing scattering angle theta
//   em2 = ehat_{perp,s} = ehat_{para,s} cross khat_s

// For special case of |k_s cross k_0| = 0  
//   (i.e., k_s dot k_0 = +/- |k_0|^2 , or k_s = +/- k_0: 
//   forward scattering or 180deg backscattering):
//   Let k_sperp = component of k_s perpendicular to target normal
//   If |k_sperp| = 0, then incident and scattered rays are normal
//                     to target.  For this special case, let
//                     
//                     ehat_spar =  ehat_01
//                     ehat_sperp = ehat_spar cross khat_s
//                               
//   If |k_s cross k_0} > 0, then
//                                k_s cross k_0 
//                     e_sperp = ---------------
//                               |k_s cross k_0|
//
// For all cases, we take:
//                     e_iperp = e_sperp
//                     e_iparr = (k_0 cross e_iperp) / |k_0|
//                     e_sparr = (k_s cross e_sperp) / |k_s|

// Here we evaluate em1r(1-3,J) = (x,y,z) components of e_sparr in TF
//                  em2r(1-3,J) = (x,y,z) components of e_sperp in TF

	for(j=0; j<nscat; ++j)
	{
// determine which case we have
		real costheta = (real)0.;
		for(k=0; k<3; ++k)
		{
			costheta += aksr[j].data[k] * en0r.data[k];
		}
		costheta /= ak1;

		if (Fabs(costheta) <= (real)0.9999)
		{
// scattering angle is neither zero nor pi
// scattering plane is therefore defined by k_0 and k_s
// obtain EM2R(1-3,J) = components of em2 = (ks cross k0)/|ks cross k0|= in TF
			em2r[j] = (aksr[j] * en0r);
			em2r[j].Ortify();
		}
		else								// Scattering angle is either zero or pi
		{
			switch(jpbc)
			{
			case PeriodicY:
			case PeriodicZ:
// JPBC=1 or 2 scattering angle = 0  -> ORDERM=0 and PHIN=0
//             scattering angle = pi -> ORDERM=0 and PHIN=pi and alpha=pi/2
// Check to see whether scattering angle theta = 0 or pi
				if (costheta > 0)		// forward scattering: em2 = [chat - nhat (nhat dot chat)]/sin(alpha)
				{
					real sinalphai, cosalphai;
					if (jpbc == PeriodicY)
					{
						cosalphai = en0r.data[1];
						sinalphai = Sqrt((real)1. - cosalphai * cosalphai);
						em2r[j].Set(en0r.data[0] * cosalphai / sinalphai, (en0r.data[1] * cosalphai - (real)1.) / sinalphai, en0r.data[2] * cosalphai / sinalphai);
					}
					else
					{
						cosalphai = en0r.data[2];
						sinalphai = Sqrt((real)1. - cosalphai * cosalphai);
						em2r[j].data[0] =  en0r.data[0] * cosalphai / sinalphai;
						em2r[j].data[1] =  en0r.data[1] * cosalphai / sinalphai;
						em2r[j].data[2] = (en0r.data[2] - (real)1.) / sinalphai;
					}
				}
				else					// backscattering: em2 = chat
				{
					if (jpbc == PeriodicY)
                        em2r[j].Set((real)0., -(real)1., (real)0.);
					else
						em2r[j].Set((real)0., (real)0., -(real)1.);
				}
				break;

			case PeriodicBoth:						// Check to see if k_0 is normal to target plane.
				if ((en0r.data[1] * en0r.data[1] + en0r.data[2] * en0r.data[2]) < 0.0001)
				{
// incident and scattered ray are both normal to the target plane:
					for(k=0; k<3; ++k)
					{
						em1r[j].data[k] = cxe01r.data[k].re;
					}
					em1r[j].Ortify();
					em2r[j] = em1r[j] * aksr[j];
					em2r[j] /= ak1;
				}
				else
				{
// incident ray is not normal to the target plane.  Therefore we can
// use the reflected ray to construct a scattering plane, and use this
// to define the perp polarization direction.
// Note that AKSR for the reflected ray differs from AKSR for the
// transmitted ray only by change in sign of AKSR(1,J)
					em2r[j].data[0] =  aksr[j].data[1] * en0r.data[2] - aksr[j].data[2] * en0r.data[1];
					em2r[j].data[1] =  aksr[j].data[2] * en0r.data[0] + aksr[j].data[0] * en0r.data[2];
					em2r[j].data[2] = -aksr[j].data[0] * en0r.data[1] - aksr[j].data[1] * en0r.data[0];
					em2r[j].Ortify();
				}
	            break;

			default:
				break;
			}
		}					// end IF(ABS(SUM).LE.0.9999_WP) ...
// Now use EM2R to obtain EM1R from e_par = khat cross e_perp
		em1r[j] = aksr[j] * em2r[j];
		em1r[j] /= ak1;
// Now evaluate scattering unit vector in Lab Frame
		ensc[j].Set(aksr[j].Scalar(xlr), aksr[j].Scalar(ylr), aksr[j].Scalar(zlr));
		ensc[j] /= ak1;
// Evaluate scattering polarization vectors in Lab Frame
		em1[j].Set(em1r[j].Scalar(xlr), em1r[j].Scalar(ylr), em1r[j].Scalar(zlr));
		em2[j].Set(em2r[j].Scalar(xlr), em2r[j].Scalar(ylr), em2r[j].Scalar(zlr));
	}
}
