#include "StdAfx.h"
#include "ReadnfHolder.h"

ReadnfHolder::ReadnfHolder(void)
{
	cxeinc = NULL;
	cxbinc = NULL;
	cxbsca = NULL;
	cxeps = cxpol = cxesca = cxadia = NULL;
	icomp = NULL;
}

ReadnfHolder::~ReadnfHolder(void)
{
	CleanDelete(icomp);
	CleanDelete2(cxeinc);
	CleanDelete2(cxpol);
	CleanDelete2(cxeps);
	CleanDelete2(cxesca);
	CleanDelete2(cxadia);
	CleanDelete2(cxbinc);
	CleanDelete2(cxbsca);
}

bool ReadnfHolder::Readnf(const char *cflename)
{
/* **
subroutine READNF
purpose: to read file with stored polarization and EM field
         produced by DDSCAT v7.3.0
         to support postprocessing
given:
   CFLENAME = name of file with stored polarization and EM field
   IDVOUT   = unit number for output (e.g., 7)

returns

via arguments:
   CSTAMP     = CHARACTER*26 string with DDSCAT version used to create file
                CFLENAME
                e.g., 'DDSCAT 7.3.0 [13.03.18]'
   VERSNUM    = integer string with version number
                e.g., 730 for version 7.3.0
   NRFLDB     = 0 if scattered magnetic field was not stored in file CFLENAME
              = 1 if scattered magnetic field was stored in file CFLENAME
   AEFF       = effective radius of target (physical units)
   DPHYS      = d = dipole spacing (physical units)
   NX,NY,NZ   = dimensions of computational volume
   NAT0       = number of occupied sites (dipoles)
   X0(3)      = x/d,y/d,z/d for site with indices (I,J,K)=(0,0,0)
   XMIN,XMAX  = x_min,x_max (physical units) for computational volume (in TF)
   YMIN,YMAX  = y_min,y_max (physical units) for computational volume (in TF)
   ZMIN,ZMAX  = z_min,z_max (physical units) for computational volume (in TF)
   NAMBIENT   = refractive index of ambient medium (real)
   WAVE       = wavelength **in vacuo** (physical units)
   AK_TF(3)   = (k_x,k_y,k_z)*d in the Target Frame, where (k_x,k_y,k_z) = 2*pi*NAMBIENT/WAVE = wavevector in the medium
   CXE0_TF(3) = complex (E_x,E_y,E_z)_TF at (x,y,z)_TF=0 and t=0 for the incident plane wave
   CXB0_TF(3) = complex (B_x,B_y,B_z)_TF at (x,y,z)_TF=0 and t=0 for the incident plane wave
   NCOMP      = number of compositions

via module READNF_ECOM:
   CXADIA(J,3)= complex diagonal elements of the "A matrix"
              = d^3/diagonal elements of complex polarizability tensor
                at lattice sites J=1-NX*NY*NZ
   CXEINC(J,3)= complex "macroscopic" (E_x,E_y,E_z)_TF of incident wave at lattice  site J
   CXEPS(IC)  = complex dielectric function for composition IC=1-NCOMP
   CXESCA(J,3)= complex "macroscopic" (E_x,E_y,E_z)_TF of scattered wave at lattice site J
   CXPOL(J,3) = complex polarization (P_x,P_y,P_z)_TF of dipole J
   ICOMP(J,3) = integer*2 composition identifier for lattice site J
              = 0 for ambient medium

via module READNF_BCOM (only if NRFLDB=1):
   CXBINC(J,3)= complex (B_x,B_y,B_z)_TF of incident wave at lattice site J
   CXBSCA(J,3)= complex (B_x,B_y,B_z)_TF of scattered wave at lattice site J

NB: macroscopic and microscopic E fields are related by
                  3
   E_macro = ----------- * E_micro
             (epsilon+2)

The dipoles respond to E_micro : P = alpha*E_micro

B.T. Draine, Princeton University, 2013.03.20

Fortran history records removed.
** */

//              >>>>> Important Note! <<<<<
// The structure of the READ statements below *must* conform to the
// structure of the corresponding WRITE statements in nearfield.f90 
// Any changes must be made in both modules.

#ifdef _WIN32
	int mode = O_RDONLY|O_BINARY;
#else	
	int mode = O_RDONLY;
#endif	
	int iobinFile = open(cflename, mode);
	if (iobinFile == -1)
	{
		printf("Cannot open file %s\n", cflename);
		return false;
	}
	else
	{
		char Buffer[256];
		int versnum;
		int i = 0;
		while(1)
		{
			read(iobinFile, Buffer+i, sizeof(char));
			if (!Buffer[i])
				break;
			++i;
		}
		read(iobinFile, &versnum, sizeof(int));
		fprintf(stdout, ">READNF data from %s\n", Buffer);
		if (versnum != 730)
		{
			fprintf(stderr, "File = %s  was written by version = %d\n", Buffer, versnum);
			fprintf(stderr, "File is incompatible with present version of subroutine Readnf\n");
			close(iobinFile);
			return false;
		}
		else
		{
			int nrword_nf;
			read(iobinFile, &nrword_nf, sizeof(int));
			if (nrword_nf != sizeof(real))
			{
				fprintf(stderr, "READNF fatal error:\n");
				fprintf(stderr, "  Word length = %d in file %s\n", nrword_nf, Buffer);
				fprintf(stderr, "  Word length = %d in subroutine Readnf\n", sizeof(real));
				close(iobinFile);
				return false;
			}
			read(iobinFile, &nrfldb, sizeof(int));
			read(iobinFile, &nxyz, sizeof(int));
			read(iobinFile, &nat0, sizeof(int));
			read(iobinFile, &nat3, sizeof(int));
			read(iobinFile, &ncomp, sizeof(int));
			read(iobinFile, &nx, sizeof(int));
			read(iobinFile, &ny, sizeof(int));
			read(iobinFile, &nz, sizeof(int));
			x0.Read(iobinFile);
			read(iobinFile, &aeff, sizeof(real));
			read(iobinFile, &nambient, sizeof(real));
			read(iobinFile, &wave, sizeof(real));
			ak_tf.Read(iobinFile);
			cxe0_tf.Read(iobinFile);
			cxb0_tf.Read(iobinFile);
			AllocateEArrays();
			read(iobinFile, &cxeps, ncomp * sizeof(Complex));
			for(int j=0; j<nxyz; ++j)							// TOREFACT
				icomp->ReadLine(iobinFile, j);
			read(iobinFile, cxpol, nat3 * sizeof(Complex));
			read(iobinFile, cxesca, nat3 * sizeof(Complex));
			read(iobinFile, cxadia, nat3 * sizeof(Complex));
			if (nrfldb == 1)
			{
				AllocateBArrays();
				read(iobinFile, cxbsca, nat3 * sizeof(Complex));
			}
			close(iobinFile);
		}
//
// compute phase for JX=JY=JZ=0
		real phi0 = ak_tf.Scalar(x0);
		for(int jx=0; jx<nx; ++jx)
		{
			real phix = jx * ak_tf.data[0] + phi0;
			for(int jy=0; jy<ny; ++jy)
			{
				real phiy = phix + jy * ak_tf.data[1];
				for(int jz=0; jz<nz; ++jz)
				{
					real arg = phiy + jz*ak_tf.data[2];
					Complex cxphas = Complex(Cos(arg), Sin(arg));
					int index = (jx * ny + jy) * nz + jz;
// compute E_macro contributed by incident wave
					for(int k=0; k<3; ++k)
					{
						short ic = icomp->Value(index, k);
						Complex cxfac((real)1., (real)0.);
						if (ic > 0) 
							cxfac = Complex((real)3., (real)0.) / (cxeps[ic-1] + (real)2.);
						cxeinc[3*index + k] = cxfac * cxe0_tf.data[k] * cxphas;
					}
					if(nrfldb == 1)
					{
						for(int k=0; k<3; ++k)
						{
							cxbinc[3*index + k] = cxb0_tf.data[k] * cxphas;
						}
					}
				}
			}
		}
		dphys = aeff * Pow(((real)4. * Pi / nat0 / (real)3.), (real)1./(real)3.);
//		int nxy = nx * ny;
//		int nzy = nz * ny;
		xmin = (real)((x0.data[0] + (real)1. - 0.5001) * dphys); 
		xmax = (real)((x0.data[0] +       nx + 0.5001) * dphys);
		ymin = (real)((x0.data[1] + (real)1. - 0.5001) * dphys);
		ymax = (real)((x0.data[1] +       ny + 0.5001) * dphys);
		zmin = (real)((x0.data[2] + (real)1. - 0.5001) * dphys);
		zmax = (real)((x0.data[2] +       nz + 0.5001) * dphys);
//
		printf(">READNF lambda= %lf physical units\n", wave);
		printf(">READNF %11.4e = AEFF (vol. equiv. radius, phys. units)\n", aeff);
		printf(">READNF %11d  = NAT0 (number of physical dipoles in target)\n", nat0);
        printf(">READNF %11.4e = d = interdipole separation (phys. units)\n", dphys);
        printf(">READNF %11.4e %11.4lf = target xmin,xmax\n", xmin, xmax);
		printf(">READNF %11.4e %11.4lf = target ymin,ymax\n", ymin, ymax);
		printf(">READNF %11.4e %11.4lf = target zmin,zmax\n", zmin, zmax);
// check solution
// this check neglects the off-diagonal elements of A
// this will not be a valid assumption for anisotropic materials
// that do not have optical axes aligned with TF xyz axes
		real einc2 = (real)0.;
		for(int k=0; k<3; ++k)
		{
			einc2 += (cxe0_tf.data[k] * cxe0_tf.data[k].conjg()).re;
		}
//
		real sumerr2 = (real)0.;
		int j1 = 0;
		for(int jx=0; jx<nx; ++jx)
		{
			for(int jy=0; jy<ny; ++jy)
			{
				for(int jz=0; jz<nz; ++jz)
				{
					int index = (jx * ny + jy) * nz + jz;
					if (icomp->Value(index, 0) > 0)
					{
						for(int k=0; k<3; ++k)
						{
// CXEINC and CXESCA are macroscopic E fields.
// DDA equation is for microscopic E field
// E_micro = E_macro*(epsilon+2)/3
							short ic = icomp->Value(index, k);
							Complex cxerr = cxpol[3*index + k] * cxadia[3*index + k] - (cxeinc[3*index + k] + cxesca[3*index + k]) * (cxeps[ic-1] + (real)2.) / (real)3.;
							sumerr2 += (cxerr * cxerr.conjg()).re;
						}
// count number of occupied sites as a sanity check 
						++j1;
					}
				}
			}
		}
// *** sanity check
		if (j1 != nat0)
		{
			fprintf(stderr, "Readnf sanity failure: inconsistent j1= %d and nat0 = %d\n", j1, nat0);
			return false;	
		}
// ***
		sumerr2 /= (einc2 * nat0);
// write information to unit IDVOUT (usually a log file)
		fprintf(stdout, ">READNF %11.4lf = normalized error |P/alpha-E|^2/|E_inc|^2\n", sumerr2);
		fprintf(stdout, ">READNF %11.4lf = AEFF (vol. equiv. radius, phys. units)\n", aeff);
		fprintf(stdout, ">READNF %d	= NAT0 (number of physical dipoles in target)\n", nat0);
		fprintf(stdout, ">READNF %11.4lf = d = interdipole separation (phys. units)\n", dphys);
		fprintf(stdout, ">READNF %11.4lf = wavelength in vacuo (phys. units)\n", wave);
		fprintf(stdout, ">READNF %11.4lf = wavelength in ambient medium (phys. units)\n", wave/nambient);

		return true;
	}
}

void ReadnfHolder::AllocateEArrays()
{
	cxeps = new Complex[ncomp];
	cxpol = new Complex[nat3];
	cxesca = new Complex[nat3];
	cxadia = new Complex[nat3];
	cxeinc = new Complex[nat3];

	icomp = new Array2F<short>;			icomp->Dimension(nxyz, 3);
}

void ReadnfHolder::AllocateBArrays()
{
	cxbinc = new Complex[nat3];
	cxbsca = new Complex[nat3];
}

bool ReadnfHolder::IsPointInside(real xa, real xb, real ya, real yb, real za, real zb)
{
	return (xa >= xmin) && (xb <= xmax) && (ya >= ymin) && (yb <= ymax) && (za >= zmin) && (zb <= zmax);
}

void ReadnfHolder::Output(FILE *outputFile, const Vect3<real> &s_inc, real snorm)
{
	fprintf(outputFile, "%10.4e = a_eff (radius of equal-volume sphere)\n", aeff);
	fprintf(outputFile, "%10.4e = d = dipole spacing\n", dphys);
	fprintf(outputFile, "%10d = N = number of dipoles in target or TUC\n", nat0);
	fprintf(outputFile, "%8d%8d%8d = NX,NY,NZ = extent of computational volume\n", nx, ny, nz);
	fprintf(outputFile, "%10.5e = wavelength in vacuo\n", wave);
	fprintf(outputFile, "%10.5lf = refractive index of ambient medium\n", nambient);
	ak_tf.Fprintf(outputFile, "%10.5lf", NULL, " = k_{inc,TF} * d\n");
	fprintf(outputFile, "%10.5lf %10.5lf = Re(E_inc,x) Im(E_inc,x) at x_TF=0,y_TF=0,z_TF=0\n", cxe0_tf.data[0].re, cxe0_tf.data[0].im);
	fprintf(outputFile, "%10.5lf %10.5lf = Re(E_inc,y) Im(E_inc,y) \"\n", cxe0_tf.data[1].re,  cxe0_tf.data[1].im);
	fprintf(outputFile, "%10.5lf %10.5lf = Re(E_inc,z) Im(E_inc,z) \"\n", cxe0_tf.data[2].re,  cxe0_tf.data[2].im);
	fprintf(outputFile, "%10.5lf %10.5lf = Re(B_inc,x) Im(B_inc,x) at x_TF=0,y_TF=0,z_TF=0\n", cxb0_tf.data[0].re, cxb0_tf.data[0].im);
	fprintf(outputFile, "%10.5lf %10.5lf = Re(B_inc,y) Im(B_inc,y) \"\n", cxb0_tf.data[1].re,  cxb0_tf.data[1].im);
	fprintf(outputFile, "%10.5lf %10.5lf = Re(B_inc,z) Im(B_inc,z) \"\n", cxb0_tf.data[2].re,  cxb0_tf.data[2].im);
	fprintf(outputFile, "%10.5lf %10.5lf %10.5lf = 8*pi*<S_inc>/c (S_inc=incident Poynting vector)\n", s_inc.data[0], s_inc.data[1], s_inc.data[2]);
	fprintf(outputFile, "%10.5lf = 8*pi*<S_inc>/c (<|S_inc|> = Snorm\n", snorm);
	if (nrfldb == 0)
	{
		fprintf(outputFile, "   x_TF       y_TF       z_TF      Re(E_x)   Im(E_x)   Re(E_y)   Im(E_y)   Re(E_z)   Im(E_z)\n");
	}
	else
	{
		fprintf(outputFile, "                                 ----------------------- E field -------------------------------------------------- B field ------------------------------ -normalized Poynting vector-\n");
		fprintf(outputFile, "   x_TF       y_TF       z_TF      Re(E_x)   Im(E_x)   Re(E_y)   Im(E_y)   Re(E_z)   Im(E_z)   Re(B_x)   Im(B_x)   Re(B_y)   Im(B_y)   Re(B_z)   Im(B_z) -- <(Sx,Sy,Sz)>/|<S_inc>| --\n");
	}
}

void ReadnfHolder::Evaluate(FILE *outputFile, real xa, real xb, real ya, real yb, real za, real zb, int nab, Array3F <Vect3<real> > &s)
{
	const real onex_ = (real)1.;
	Complex cxe_inc[3], cxe_sca[3], cxp[3], cxe[3], cxb_inc[3], cxb_sca[3], cxb[3], svec[3];
	int jj[8], j;
	real w[8], xtf[3];
	int nzy = ny * nz;

	for(int ja=0; ja<=nab; ++ja)
	{
		real zeta = (real)ja / (real)nab;
		xtf[0] = xa + (xb - xa) * zeta;
		xtf[1] = ya + (yb - ya) * zeta;
		xtf[2] = za + (zb - za) * zeta;
//
// XTF is assumed to be in physical units
// XTF/DPHYS is in dipole units
// I + X0 = XTF/DPHYS
// I = XTF/DPHYS - X0

		int ix1 = (int)(xtf[0] / dphys - x0.data[0]);
		int iy1 = (int)(xtf[1] / dphys - x0.data[1]);
		int iz1 = (int)(xtf[2] / dphys - x0.data[2]);
		if (ix1 < 1) ix1 = 1;
		if (ix1 >= nx) ix1 = nx - 1;
		if (iy1 < 1) iy1 = 1;
		if (iy1 >= ny) iy1 = ny - 1;
		if (iz1 < 1) iz1 = 1;
		if (iz1 >= nz) iz1 = nz - 1;
//
// -------- determine indices ----------------- zyx
		jj[0] = (iz1 - 1) + nz * (iy1 - 1) + nzy * (ix1 - 1);		// 000
		jj[1] = jj[0] + nzy;										// 001
		jj[2] = jj[0] + nz;											// 010
		jj[3] = jj[2] + nzy;										// 011
		jj[4] = jj[0] + 1;											// 100
		jj[5] = jj[4] + nzy;										// 101
		jj[6] = jj[4] + nz;											// 110
		jj[7] = jj[6] + nzy;										// 111
		real wx = xtf[0] / dphys - x0.data[0] - ix1;
		real wy = xtf[1] / dphys - x0.data[1] - iy1;
		real wz = xtf[2] / dphys - x0.data[2] - iz1;
//
// ------------ determine weights ---------- zyx
		w[0] = (onex_ - wx) * (onex_ - wy) * (onex_ - wz);			// 000
		w[1] = wx * (onex_ - wy) * (onex_ - wz);					// 001
		w[2] = (onex_ - wx) * wy * (onex_ - wz);					// 010
		w[3] = wx * wy * (onex_ - wz);								// 011
		w[4] = (onex_ - wx) * (onex_ - wy) * wz;					// 100
		w[5] = wx * (onex_ - wy) * wz;								// 101
		w[6] = (onex_ - wx) * wy * wz;								// 110
		w[7] = wx * wy * wz;										// 111
//
// -------- evaluate weighted averages -----------------------------------
		for(j=0; j<3; ++j)
		{
			cxe_inc[j].clear();
			cxe_sca[j].clear();
			cxp[j].clear();
			for(int ji=0; ji<8; ++ji)
			{
				int index = 3 * jj[ji];
				cxe_inc[j] += (cxeinc[index + j] * w[ji]);
				cxe_sca[j] += (cxesca[index + j] * w[ji]);
				cxp[j]     += ( cxpol[index + j] * w[ji]);
			}
		}
		for(j=0; j<3; ++j)
		{
			cxe[j] = cxe_inc[j] + cxe_sca[j];
		}
		switch(nrfldb)
		{
		case 0:
//
// write out total E field at point on track
			fprintf(outputFile, "%11.3e%11.3e%11.3e%10.5lf%10.5lf%10.5lf%10.5lf%10.5lf%10.5lf\n", 
				xtf[0], xtf[1], xtf[2], cxe[0].re, cxe[0].im, cxe[1].re, cxe[1].im, cxe[2].re, cxe[2].im);
				break;

		case 1:
			for(j=0; j<3; ++j)
			{
				cxb_inc[j].clear();
				cxb_sca[j].clear();
				for(int ji=0; ji<8; ++ji)
				{
					int index = 3 * jj[ji];
					cxb_inc[j] += (cxbinc[index + j] * w[ji]);
					cxb_sca[j] += (cxbsca[index + j] * w[ji]);
				}
			}
// calculate time-averaged Poynting vector at each point, normalized by
// magnitude of time-averaged Poynting vector of incident plane wave
			for(j=0; j<3; ++j)
			{
				cxb[j] = cxb_inc[j] + cxb_sca[j];
				svec[j].clear();
				for(int ji=0; ji<8; ++ji)
				{
					svec[j] += s.Value(ix1, iy1, iz1).data[ji] * w[ji];
				}
			}

			fprintf(outputFile, "%11.3e%11.3e%11.3e%10.5lf%10.5lf%10.5lf%10.5lf%10.5lf%10.5lf%10.5lf%10.5lf%10.5lf%10.5lf%10.5lf%10.5lf%10.5lf%10.5lf%10.5lf%10.5lf%10.5lf%10.5lf\n", 
				xtf[0], xtf[1], xtf[2], cxe[0].re, cxe[0].im, cxe[1].re, cxe[1].im, cxe[2].re, cxe[2].im, cxb[0].re, cxb[0].im, cxb[1].re, cxb[1].im, cxb[2].re, cxb[2].im, 
				svec[0].re, svec[0].im, svec[1].re, svec[1].im, svec[2].re, svec[2].im);
				break;

			default:
				fprintf(stderr, "DDpostprocess fatal error: NRFLDB = %d\n", nrfldb);
				break;
		}
	}
}
