#include "StdAfx.h"

#include "Definitions.h"
#include "DDscatMain.h"
#include "DDscatParameters.h"
#include "Matc3.h"

void Rotate(Vect3<real> &a1, Vect3<real> &a2, real ak1, real beta, real theta, real phi, 
	Vect3<real> &en0r, Vect3<Complex> &cxe01r, Vect3<Complex> &cxe02r, int nscat, Vect3<real> *ensc, 
	Vect3<real> *em1, Vect3<real> *em2, Vect3<real> *aksr, Vect3<real> *em1r, Vect3<real> *em2r)
{
/* **
Given:
       A1(1-3)=axis 1 of (original, unrotated) target
       A2(1-3)=axis 2 of (original, unrotated) target
       AK1=magnitude of k vector
       CXE01(1-3)=original incident polarization vector 1
       CXE02(1-3)=original incident polarization vector 2
       BETA=angle (radians) through which target is to be rotated around target axis A1
       THETA=angle (radians)which target axis A1 is to make with Lab x-axis
       PHI=angle (radians) between Lab x,A1 plane and Lab x,y plane
       ENSC(1-3,1-NSCAT)=scattering directions in Lab Frame
       EM1(1-3,1-NSCAT)=scattering polarization vector 1 in Lab Frame
       EM2(1-3,1-NSCAT)=scattering polarization vector 2 in Lab Frame

Returns:
       EN0R(1-3)=incident propagation vector in Target Frame
       CXE01R(1-3)=incident polariz. vector 1 in Target Frame
       CXE02R(1-3)=incident polariz. vector 2 in Target Frame
       AKSR(1-3,1-NSCAT)=scattering k vectors in Target Frame
       EM1R(1-3,1-NSCAT)=scattering pol. vector 1 in Target Frame
       EM2R(1-3,1-NSCAT)=scattering pol. vector 2 in Target Frame

Purpose:
       In the Lab Frame, we hold the propagation direction (x-axis)
       and incident polarization vectors (CXE01 and CXE02) fixed,
       and rotate target to an orientation specified by the three
       angles BETA,THETA,PHI.
       However, computations in the main program DDSCAT are carried out
       in the Target Frame, in which the target is held fixed and the
       propagation vectors and polarization vectors are rotated.
       Hence, this routine computes the incident propagation vector,
       scattering propagation vectors, and associated polarization
       vectors in the Target Frame.

B. T. Draine, Princeton Univ. Obs., 89.11.20

Fortran history records removed.

Copyright (C) 1996,2003,2005 B.T. Draine and P.J. Flatau
Copyright (C) 2012, C++ version, Choliy V.

This code is covered by the GNU General Public License.
** */

//  Determine initial orientation THETA0,PHI0,BETA0 of target:  theory:
//
//      theta = arccos(a_1x)
//
//  if sin(theta) .ne. 0 , then
//
//      cos(phi)=a_1y/sin(theta)      sin(phi)=a_1z/sin(theta)
//      if(a_1z/sin(theta) > 0) then phi = arccos(a_1y/sin(theta))
//      if(a_1z/sin(theta) < 0) then phi = 2*pi-arccos(a_1y/sin(theta))
//
//      cos(beta)=-a_2x/sin(theta)    sin(beta)=a_3x/sin(theta)
//      if(a_3x/sin(theta) > 0) then beta = arccos(-a_2x/sin(theta))
//      if(a_3x/sin(theta) < 0) then beta = 2*pi-arccos(-a_2x/sin(theta))
//
//  note: a_3x = a_1y*a_2z - a_1z*a_2y
//
//  if sin(theta) = 0 , then
//
//      phi = 0
//
//      cos(beta)=a_2y                sin(beta)=a_2z
//      if(a_2z) > 0)      then      beta = arccos(a_2y)
//      if(a_2z) < 0)      then      beta = 2*pi-arccos(a_2y)

	const real onex_ = (real)1.;
	const real zero_ = (real)0.;
	real phi0, beta0;
	real theta0 = Acos(a1.data[0]);
	real sinthe = Sin(theta0);
	if (sinthe != zero_)
	{
// guard against roundoff errors:
		real term = a1.data[1] / sinthe;
		if (term >  onex_) term =  onex_;
		if (term < -onex_) term = -onex_;
//
		if (a1.data[2] >= zero_)
			phi0 = Acos(term);
		else
			phi0 = TwoPi - Acos(term);
// guard against roundoff errors:
		term = -a2.data[0] / sinthe;
		if (term >  onex_) term =  onex_;
		if (term < -onex_) term = -onex_;
//
		if (a1.data[1]*a2.data[2] - a1.data[2]*a2.data[1] >= zero_)
			beta0 = Acos(term);
        else
			beta0 = TwoPi - Acos(term);
	}
	else
	{
		phi0 = zero_;
		if (a2.data[2] >= zero_)
			beta0 = Acos(a2.data[1]);
		else
			beta0 = TwoPi - Acos(a2.data[1]);
	}
//
// First rotate the target through angle PHI-PHI0 around x axis (i.e., rotate the lab through angle PHI0-PHI around x axis)
// obtain rotation matrix:
	Vect3<real> vec;
	Matc3<real> rm1, rm2, rm3;
	vec.Set(onex_, zero_, zero_);
	Rot2(vec, phi0-phi, rm1);
//
// Now rotate the target through angle THETA-THETA0 around axis VEC = xaxis cross a1
// (i.e., rotate the lab through angle THETA0-THETA around xaxis cross a1 )
// compute rotation axis
	vec.Set(zero_, -a1.data[2], a1.data[1]);
//
// For special case where a1=xaxis, we want rotation axis=a1 cross a2
	if ((vec.data[1] * vec.data[1] + vec.data[2] * vec.data[2]) < (real)1.e-10)
	{
		vec.data[1] = a1.data[2]*a2.data[0] - a1.data[0]*a2.data[2];
		vec.data[2] = a1.data[0]*a2.data[1] - a1.data[1]*a2.data[0];
	}

	Rot2(vec, theta0-theta, rm2);
	rm3 = rm2 * rm1;
//
// Now rotate the target through angle BETA-BETA0 around a1 (i.e., rotate the lab through angle BETA0-BETA around a1)
	Rot2(a1, beta0-beta, rm1);
	rm2 = rm1 * rm3;
//
// RM2 is now the full rotation matrix
//
// Now rotate all required vectors:
// Incident propagation vector:
	vec.Set(onex_, zero_, zero_);
	en0r = rm2 * vec;
//
// Incident polarization vectors:
	Prod3c(rm2, DDscatParameters::GetInstance()->Cxe01_lf(), cxe01r);
	Prod3c(rm2, DDscatParameters::GetInstance()->Cxe02_lf(), cxe02r);
//
// Scattering vectors:
	Prod3v(rm2, ensc, aksr, nscat);
	for(int i=0; i<nscat; ++i)
	{
		aksr[i] *= ak1;
	}
//
// Scattering polarization vectors:
	Prod3v(rm2, em1, em1r, nscat);
	Prod3v(rm2, em2, em2r, nscat);
}
