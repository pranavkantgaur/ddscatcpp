#include "StdAfx.h"

#include "Definitions.h"
#include "Vect3.h"
#include "Matc3.h"
#include "Complex.h"

void Rot2(const Vect3<real> &a, real theta, Matc3<real> &rm)
{
/* ** 
This subroutine finds the (3 by 3) rotation matrix for the given axis (vector) - a and the rotation angle - theta.

On input:
a(3)    --- axis of rotation  (vector, any length)
theta   --- angle of rotation (in radians) around vector a

On output:
rm(3,3) --- rotation matrix such that new vector v2 is obtained from the old vector v1 by the transformation: 
            v2 = rm \dot v1
(If more rotations are needed: generate matrices rm with the help of this subroutine and (matrix) multiply them before performing v2 = rm \dot v1.)

Reference:
Leubner, C., 1977,  Coordinate-free  rotation operator,
Am. J. Phys. 47, 727---729.

History:
Written by PJF, May 1989, for use in the rotation-averaged DDA calculations.

Copyright (C) 1993, B.T. Draine and P.J. Flatau
Copyright (C) 2013, C++ version, Choliy V.
This code is covered by the GNU General Public License.
** */

	real ct = Cos(theta);
	real st = Sin(theta);
	real anorm = a.Mod();
	Vect3<real> ahat(a);
	ahat /= anorm;
//
// cos(\theta) {\bf 1}  term:
	rm.Set(ct, -st*ahat.data[2], st*ahat.data[1], st*ahat.data[2], ct, -st*ahat.data[0], -st*ahat.data[1], st*ahat.data[0], ct);
//
// aa-dyadic (outer-product):
	for(int i=0; i<3; ++i)
	{
		for(int j=0; j<3; ++j)
		{
			rm.Data(i,j) += ((real)1. - ct) * ahat.data[j] * ahat.data[i];
		}
	}
}

void Prod3c(Matc3<real> &rm, Vect3<Complex> &cvin, Vect3<Complex> &cvout)  
{
/* **
Given:
      RM=3x3 rotation matrix
      CVIN=Complex 3-vector
Computes:
      CVOUT=RM \dot CVIN

Copyright (C) 1993, B.T. Draine and P.J. Flatau
Copyright (C) 2013, C++ version, Choliy V.
This code is covered by the GNU General Public License.
** */

	for(int i=0; i<3; ++i)
	{
		cvout.data[i].clear();
		for(int j=0; j<3; ++j)
		{
			cvout.data[i] += cvin.data[j] * rm.Data(i, j);
		}
	}
}

void Prod3v(const Matc3<real> &rm, Vect3<real> *vin, Vect3<real> *vout, int nv)
{
/* **
Given:
      RM = 3x3 rotation matrix
      VIN = NV 3-vectors
Computes:
      VOUT=RM \dot VIN for each of the NV 3-vectors

Copyright (C) 1993,1996 B.T. Draine and P.J. Flatau
Copyright (C) 2013, C++ version, Choliy V.

Fortran versions history records removed.
This code is covered by the GNU General Public License.
** */
	for(int n=0; n<nv; ++n)
	{
		vout[n] = rm * vin[n];
	}
}