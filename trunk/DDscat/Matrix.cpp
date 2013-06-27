#include "StdAfx.h"

#include "Matrix.h"

Matrix::Matrix(void)
{
	cxadia = cxaoff = NULL;
	curSize = 0;
}

Matrix::~Matrix(void)
{
	Delete();
}

void Matrix::Allocate(unsigned int newSize)
{
	if (newSize != curSize)
	{
		Delete();
	}
	cxadia = new Complex[newSize];
	cxaoff = new Complex[newSize];
	curSize = newSize;
}

void Matrix::Delete()
{
	CleanDelete2(cxadia);
	CleanDelete2(cxaoff);
}

int Matrix::WriteDiagonal(int file)
{
	return write(file, cxadia, curSize * sizeof(Complex));
}

int Matrix::WriteOffDiagonal(int file)
{
	return write(file, cxaoff, curSize * sizeof(Complex));
}

int Matrix::ReadDiagonal(int file)
{
	return read(file, cxadia, curSize * sizeof(Complex));
}

int Matrix::ReadOffDiagonal(int file)
{
	return read(file, cxaoff, curSize * sizeof(Complex));
}

void Matrix::Debug(FILE *file, int from, int toto)
{
	if ((from == 0) && (toto == 0))
		toto = curSize;
	for(int i=from; i<toto; ++i)
	{
		fprintf(file, "%3d %20.8e %20.8e %20.8e %20.8e\n", i, cxadia[i].re, cxadia[i].im, cxaoff[i].re, cxaoff[i].im);
	}
}

void Matrix::Evala(Matrix *theTensor, int nat)
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

	Complex *cxalph = theTensor->Diagonal();
	Complex *cxalof = theTensor->OffDiagonal();
	for(int j=0; j<nat; ++j)					// TODO: are they the vectorial products?
	{
		int index = 3*j;
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
	}
}

void Matrix::Evalq(DipoleData *theDipoleData, Vect3<real> &ak, int nat3, real e02, real &cabs, real &cext, real &cpha, int imethd)
{
/* **
Given: 
	CXADIA(J,1-3)=(a_11,a_22,a_33) for dipoles J=1-NAT, where symmetric 3x3 matrix a_ij is inverse of Complex polarizability tensor alpha_ij for dipole J
	CXAOFF(J,1-3)=(a_23,a_31,a_12) for dipoles J=1-NAT
	AK(1-3) = (k_x,k_y,k_z)*d d = (d_x*d_y*d_z)**(1/3) = effective lattice spacing
	NAT3 = 3*number of dipoles
	E02 = |E_0|^2 , where E_0 = incident Complex E field.
	CXE(1-NAT3) = components of E field at each dipole, in order E_1x,E_2x,...,E_NATx,E_1y,E_2y,...,E_NATy,E_1z,E_2z,...,E_NATz
	CXP(1-NAT3) = components of polarization vector at each dipole, in order P_1x,P_2x,...,P_NATx,P_1y,P_2y,...,P_NATy,P_1z,P_2z,...,P_NATz
	IMETHD = 0 or 1

Finds:
        CEXT = extinction cross section /d**2 
  and, if IMETHD=1, also computes
        CPHA = phase-lag cross section /d**2
        CABS = absorption cross section /d**2

Note: present definition of CPHA is 1/2 of Martin's definition

B.T.Draine, Princeton Univ. Obs., 87/1/4

History:
Fortran version history removed.

Copyright (C) 1993,1997,1998,2008 B.T. Draine and P.J. Flatau
This code is covered by the GNU General Public License.
** */

//
// Compute magnitude of kd

	Complex *cxe = theDipoleData->Cxe_tf();
	Complex *cxp = theDipoleData->Cxxi();

	int j1;
	int nat = nat3 / 3;
	cext = cabs = cpha = (real)0.;
	real ak2 = ak.ModSquared();
	real ak1 = Sqrt(ak2);
	real ak3 = ak1 * ak2;
//
// Initialization complete.
	if (imethd == 0)
	{
		for(j1=0; j1<nat3; ++j1)
		{
			cext = cext + cxp[j1].im * cxe[j1].re - cxp[j1].re * cxe[j1].im;
		}
		cext = FourPi * ak1 * cext / e02;
	}
	else							// imethd == 1
	{
// C_abs=  (4*pi*k/|E_0|^2)*sum_J { Im[P_J*conjg(a_J*P_J)] - (2/3)*k^3*|P_J|^2 } = (4*pi*k/|E_0|^2)*sum_J {-Im[conjg(P_J)*a_J*P_J] - Im[i*(2/3)*k^3*P_J*conjg(P_J)]}
//      = -(4*pi*k/|E_0|^2)*Im{ sum_J [ conjg(P_J)*(a_J*P_J + i*(2/3)*k^3*P_J) ] }
		Complex cxa;
		Complex rabs((real)0., ak3 / 1.5);
		for(j1=0; j1<nat; ++j1)
		{
			int index = 3*j1;
			Complex dcxa = cxp[index  ].conjg() * ((cxadia[index  ] + rabs) * cxp[index  ] + cxaoff[index+1] * cxp[index+2] + cxaoff[index+2] * cxp[index+1]) +
				           cxp[index+1].conjg() * ((cxadia[index+1] + rabs) * cxp[index+1] + cxaoff[index+2] * cxp[index  ] + cxaoff[index  ] * cxp[index+2]) + 
				           cxp[index+2].conjg() * ((cxadia[index+2] + rabs) * cxp[index+2] + cxaoff[index  ] * cxp[index+1] + cxaoff[index+1] * cxp[index  ]);
			cxa += dcxa;
		}
		cabs = -FourPi * ak1 * cxa.im / e02;
//
		cxa.clear();
		for(j1=0; j1<nat3; ++j1)
		{
			cxa += cxp[j1] * cxe[j1].conjg();
		}
		cext = FourPi * ak1 * cxa.im / e02;
		cpha =  TwoPi * ak1 * cxa.re / e02;
	}
}

void Matrix::LoadableAlphadiag(LoadableTarget *loadableTarget)
{
	if (loadableTarget->GetDFdata() == NULL)
		return;

	Complex s1, s2, s3, s4, s5, s6;
	real r[3][3];

	int curTargetNat = loadableTarget->Nat();
	short *curTargetIcompData = loadableTarget->Icomp().GetData();
	Complex *cxalph = Diagonal();
	Complex *cxalof = OffDiagonal();
	AbstractDFData *dfData = loadableTarget->GetDFdata();
	for(int ia=0; ia<curTargetNat; ++ia)
	{
		int index = 3*ia;
		if (curTargetIcompData[index] > 0)
		{
			if (loadableTarget->IsDipolePhysicallyAnisotropic(ia))
//			if ((curTargetIcompData[index] != curTargetIcompData[index + 1]) || 
//				(curTargetIcompData[index] != curTargetIcompData[index + 2]))
			{
				if (dfData->AreAnglesZero(ia))
					continue;

				real costh = Cos(dfData->GetThetadf(ia));
				real cosph = Cos(dfData->GetPhidf(ia));
				real cosbe = Cos(dfData->GetBetadf(ia));
				if (costh * cosbe * cosph < (real)1.)
				{
//
// if nonzero rotation, recalculate CXALPH and CXALOF
					real sinth = Sin(dfData->GetThetadf(ia));
					real sinph = Sin(dfData->GetPhidf(ia));
					real sinbe = Sin(dfData->GetBetadf(ia));
//
// Define R = rotation matrix // ! Define RI = inverse of R
					r[0][0] =  costh;
					r[0][1] =  sinth*cosph;
					r[0][2] =  sinth*sinph;
					r[1][0] = -sinth*cosbe;
					r[1][1] =  costh*cosbe*cosph - sinbe*sinph;
					r[1][2] =  costh*cosbe*sinph + sinbe*cosph;
					r[2][0] =  sinth*sinbe;
					r[2][1] = -costh*sinbe*cosph - cosbe*sinph;
					r[2][2] = -costh*sinbe*sinph + cosbe*cosph;
//
// Calculate diagonal elements:
					s1.clear();
					s2.clear();
					s3.clear();
					s4.clear();
					s5.clear();
					s6.clear();
					for(int kl=0; kl<3; ++kl)
					{
						s1 += (cxalph[index + kl] * r[0][kl] * r[0][kl]);
						s2 += (cxalph[index + kl] * r[1][kl] * r[1][kl]);
						s3 += (cxalph[index + kl] * r[2][kl] * r[2][kl]);
						s4 += (cxalph[index + kl] * r[2][kl] * r[1][kl]);
						s5 += (cxalph[index + kl] * r[0][kl] * r[2][kl]);
						s6 += (cxalph[index + kl] * r[1][kl] * r[0][kl]);
					}
					cxalph[index    ] = s1;
					cxalph[index + 1] = s2;
					cxalph[index + 2] = s3;
					cxalof[index    ] = s4;
					cxalof[index + 1] = s5;
					cxalof[index + 2] = s6;
				}
			}
		}
	}
}

void Matrix::AbstractAlphadiag(AbstractTarget *currentTarget, Complex *dielec, real b1, real b2, real *b3, const Complex &cxrr)
{
	const real onex_ = (real)1.;
	int curTargetNat = currentTarget->Nat();
	short *curTargetIcompData = currentTarget->Icomp().GetData();
	int nncomp = currentTarget->Ncomp();
	Vect3<Complex> *cax = new Vect3<Complex>[nncomp];
	Complex *cxalph = Diagonal();
	Complex *cxalof = OffDiagonal();
//
	int ia, l;
	for(ia=0; ia<nncomp; ++ia)
	{
		Complex temp = dielec[ia];
		for(l=0; l<3; ++l)
		{
			real b33 = (b3) ? b3[l] : (real)0.;
			Complex cxterm = (temp - onex_) / (temp + (real)2.) * (real)(0.75 / Pi);			// First compute Clausius-Mossotti polarizability:
			cxterm = cxterm / (cxterm * (temp * (b2 + b33) + b1) + onex_);						// Now apply corretion term
			cax[ia].data[l] = cxterm / (cxterm * cxrr + onex_);
		}
	}
//
	for(ia=0; ia<curTargetNat; ++ia)
	{
		for(l=0; l<3; ++l)
		{
			int index = 3*ia + l;
			short ic = curTargetIcompData[index];
			if (ic > nncomp) ic = nncomp;
			if (ic > 0)
				cxalph[index] = cax[ic-1].data[l];
			else 
				cxalph[index].unityRe();														// To avoid divisions by zero, etc., set CXALPH=1 for vacuum sites.
			cxalof[index].clear();																// Set off-diagonal terms to zero
		}
	}
	delete [] cax;
}