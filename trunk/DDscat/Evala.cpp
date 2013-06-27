#include "StdAfx.h"

#include "Definitions.h"
#include "DDscatMain.h"
#include "Matrix.h"

void Evala(Matrix *theMatrix, Matrix *theTensor, int nat)
{
/* **
Purpose: 
	to evaluate 3x3 ***diagonal*** blocks of matrix A for use in ziscrete-dipole scattering calculation.

Given:   CXALPH(J,1-3)=(alpha_11,alpha_22,alpha_33) where alpha = Complex polarizability tensor of dipole J (J=1-NAT, first 3*NAT elements of CXALPH used)
          CXALOF(J,1-3)=(alpha_23,alpha_31,alpha_12) for dipole J (J=1-NAT, first 3*NAT elements of CXALOF used)
          NAT=number of dipoles
          MXN3,MXNAT=dimensioning information

 Returns: CXADIAG(1-3*NAT)=diagonal elements of matrix A
          CXAOFF(J,1-3)=(a_23,a_31,a_12) where a_ij=3x3 matrix equal to inverse of 3x3 matrix alpha_ij for dipole J (first 3*NAT elements of CXAOFF are used)

 Note that it is assumed that vectors E and P will have data
 in order E_1x,E_2x,...,E_Nx,E_1y,E_2y,...,E_Ny,E_1z,E_2z,...,E_Nz

 B. T. Draine, Princeton Univ. Obs., 87.01.07
 History:
 88.05.05 (BTD): Modifications...
 90.11.06 (BTD): Modified to pass array dimensions
 97.12.25 (BTD): Added CXAOFF,CXALOF to argument list
                 Modified to allow nondiagonal polarizabilities alpha_ij
 98.01.01 (BTD): Correct inconsistency in assumed data ordering.

 End history
 Copyright (C) 1993,1997,1998 B.T. Draine and P.J. Flatau
 This code is covered by the GNU General Public License.

Given symmetric matrix m with six independent elements a1, a2, a3, b1,
b2, b3 (see below), we seek inverse matrix M with six independent
elements A1, A2, A3, B1, B2, B3 such that

      (a1 b3 b2)   (A1 B3 B2)   (1 0 0)
      (b3 a2 b1) x (B3 A2 B1) = (0 1 0)
      (b2 b1 a3)   (B2 B1 A3)   (0 0 1)

 solution:
                     a2*a3 - b1*b1
 A1 = ----------------------------------------------
      a1*a2*a3+2*b1*b2*b3-a1*b1**2-a2*b2**2-a3*b3**2

                     a3*a1 - b2*b2
 A2 = ----------------------------------------------
      a1*a2*a3+2*b1*b2*b3-a1*b1**2-a2*b2**2-a3*b3**2

                     a1*a2 - b3*b3
 A3 = ----------------------------------------------
      a1*a2*a3+2*b1*b2*b3-a1*b1**2-a2*b2**2-a3*b3**2

                     b2*b3 - a1*b1
 B1 = ----------------------------------------------
      a1*a2*a3+2*b1*b2*b3-a1*b1**2-a2*b2**2-a3*b3**2

                     b3*b1 - a2*b2
 B2 = ----------------------------------------------
      a1*a2*a3+2*b1*b2*b3-a1*b1**2-a2*b2**2-a3*b3**2

                     b1*b2 - a3*b3
 B3 = ----------------------------------------------
      a1*a2*a3+2*b1*b2*b3-a1*b1**2-a2*b2**2-a3*b3**2
** */

	Complex *cxadia = theMatrix->Diagonal();
	Complex *cxaoff = theMatrix->OffDiagonal();
	Complex *cxalph = theTensor->Diagonal();
	Complex *cxalof = theTensor->OffDiagonal();
	for(int j=0; j<nat; ++j)					// TODO: are they the vectorial products?
	{
		int index = 3*j;
//fprintf(stderr, "%d %lf %lf %lf %lf %lf %lf\n", index, cxalph[index].re, cxalph[index].im, cxalph[index+1].re, cxalph[index+1].im, cxalph[index+2].re, cxalph[index+2].im);
//fprintf(stderr, "%d %lf %lf %lf %lf %lf %lf\n", index, cxalof[index].re, cxalof[index].im, cxalof[index+1].re, cxalof[index+1].im, cxalof[index+2].re, cxalof[index+2].im);
		Complex cx1 = cxalof[index  ] * cxalof[index  ];
		Complex cx2 = cxalof[index+1] * cxalof[index+1];
		Complex cx3 = cxalof[index+2] * cxalof[index+2];
		Complex denom = cxalph[index]*cxalph[index+1]*cxalph[index+2] + cxalof[index]*cxalof[index+1]*cxalof[index+2]*(real)2. - cxalph[index]*cx1 - cxalph[index+1]*cx2 - cxalph[index+2]*cx3;
		cxadia[index  ] = (cxalph[index+1] * cxalph[index+2] - cx1) / denom;
		cxadia[index+1] = (cxalph[index+2] * cxalph[index  ] - cx2) / denom;
		cxadia[index+2] = (cxalph[index  ] * cxalph[index+1] - cx3) / denom;
		cxaoff[index  ] = (cxalof[index+1] * cxalof[index+2] - cxalph[index  ] * cxalof[index  ]) / denom;
		cxaoff[index+1] = (cxalof[index+2] * cxalof[index  ] - cxalph[index+1] * cxalof[index+1]) / denom;
		cxaoff[index+2] = (cxalof[index  ] * cxalof[index+1] - cxalph[index+2] * cxalof[index+2]) / denom;
//fprintf(stderr, "%d %lf\n", index, denom);
//fprintf(stderr, "%d %lf %lf %lf %lf %lf %lf\n", index, cxadia[index].re, cxadia[index].im, cxadia[index+1].re, cxadia[index+1].im, cxadia[index+2].re, cxadia[index+2].im);
//fprintf(stderr, "%d %lf %lf %lf %lf %lf %lf\n", index, cxaoff[index].re, cxaoff[index].im, cxaoff[index+1].re, cxaoff[index+1].im, cxaoff[index+2].re, cxaoff[index+2].im);
	}
}
